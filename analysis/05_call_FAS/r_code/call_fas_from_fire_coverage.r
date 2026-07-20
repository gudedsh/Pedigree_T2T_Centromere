#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(zoo)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx == length(args)) stop("Missing value for ", flag)
  args[idx + 1]
}

coverage_file <- get_arg("--coverage")
cen_file <- get_arg("--cen")
out_file <- get_arg("--out", "FAS_results.tsv")
bin_out <- get_arg("--bin_out", "FAS_bin_summary.tsv")
plotdir <- get_arg("--plotdir", "FAS_QC_plots")

stat_out <-  file.path(plotdir,"FAS_cell_difference_stats.tsv")
stat_plot <- file.path(plotdir,"FAS_cell_difference_summary.pdf")
pdf_file <-file.path(plotdir,"all_chr_FAS_QC.pdf")

bin_size <- as.integer(get_arg("--bin_size", "2000"))
smooth_k <- as.integer(get_arg("--smooth_k", "10"))
cdr_expansion <- as.integer(get_arg("--cdr_expansion", "0"))
out_CDR_cutoff <- as.numeric(get_arg("--out_CDR_cutoff", "0.10"))
in_CDR_cutoff <- as.numeric(get_arg("--in_CDR_cutoff", "0.05"))

parameter_summary <- paste0("bin",bin_size,"_smooth",smooth_k,"_outCDRCutoff",out_CDR_cutoff,"_inCDRCutoff",in_CDR_cutoff,"_cdrExp",cdr_expansion)


if (is.null(coverage_file) || is.null(cen_file)) {
  stop("Usage: Rscript call_fas_from_fire_coverage.r --coverage activeHOR_coverage.tsv --cen active_hor_cdr.bed")
}
if (!file.exists(coverage_file)) stop("Coverage file not found: ", coverage_file)
if (!file.exists(cen_file)) stop("Active-HOR/CDR annotation not found: ", cen_file)
if (is.na(bin_size) || bin_size < 1) stop("--bin_size must be a positive integer")
if (is.na(smooth_k) || smooth_k < 1) stop("--smooth_k must be a positive integer")
if (is.na(cdr_expansion)) stop("--cdr_expansion must be an integer")
if (!is.finite(out_CDR_cutoff) || out_CDR_cutoff < 0 || out_CDR_cutoff >= 0.5) {
  stop("--out_CDR_cutoff must be in [0, 0.5)")
}
if (!is.finite(in_CDR_cutoff) || in_CDR_cutoff < 0 || in_CDR_cutoff >= 0.5) {
  stop("--in_CDR_cutoff must be in [0, 0.5)")
}

dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)


## ---------- read coverage ----------
dt <- fread(coverage_file)
req <- c("sample","chr","start","end","fire","linker","nucleosome")
if (!all(req %in% names(dt))) stop("Coverage file must contain: ", paste(req, collapse=", "))

dt[, `:=`(
  sample = as.character(sample),
  chr = as.character(chr),
  start = as.integer(start),
  end = as.integer(end),
  fire = as.numeric(fire),
  linker = as.numeric(linker),
  nucleosome = as.numeric(nucleosome)
)]
dt <- dt[end > start]
dt[, total := fire + linker + nucleosome]
dt <- dt[total > 0]
if (nrow(dt) == 0) stop("No valid positive-coverage records remain after input filtering")

## ---------- read active HOR / CDR ----------
cen_lines <- readLines(cen_file, warn = FALSE)
cen_lines <- cen_lines[!grepl("^[[:space:]]*(#|$)", cen_lines)]
if (length(cen_lines) == 0) stop("Active-HOR/CDR annotation contains no intervals")
cen <- fread(text = paste(cen_lines, collapse = "\n"), header = FALSE)
if (ncol(cen) < 5) stop("CEN file must have: chr active_start active_end cdr_start cdr_end")
setnames(cen, 1:5, c("chr","active_start","active_end","cdr_start","cdr_end"))

cen <- cen[, .(
  chr = as.character(chr),
  active_start = as.integer(active_start),
  active_end = as.integer(active_end),
  cdr_start = as.integer(cdr_start),
  cdr_end = as.integer(cdr_end)
)]

cen[, `:=`(
  cdr_start_exp = pmax(active_start, cdr_start - cdr_expansion),
  cdr_end_exp = pmin(active_end, cdr_end + cdr_expansion),
  region_id = paste0(chr, ":", active_start, "-", active_end)
)]
if (cen[is.na(active_start) | is.na(active_end) | is.na(cdr_start) | is.na(cdr_end) |
        active_start < 0 | active_end <= active_start | cdr_start < active_start |
        cdr_end > active_end | cdr_end <= cdr_start, .N] > 0) {
  stop("Active-HOR/CDR annotation contains invalid coordinates")
}
if (cen[cdr_end_exp <= cdr_start_exp, .N] > 0) {
  stop("CDR expansion produces an invalid interval")
}

## ---------- make bins ----------
bins <- rbindlist(lapply(seq_len(nrow(cen)), function(i) {
  starts <- seq(cen$active_start[i], cen$active_end[i] - 1, by = bin_size)
  data.table(
    chr = cen$chr[i],
    start = starts,
    end = pmin(starts + bin_size, cen$active_end[i]),
    cdr_start = cen$cdr_start_exp[i],
    cdr_end = cen$cdr_end_exp[i],
    region_id = cen$region_id[i]
  )
}))

bins[, `:=`(
  in_cdr = start < cdr_end & end > cdr_start,
  bin_width = end - start
)]

## ---------- length-weighted binning ----------
setkey(dt, chr, start, end)
setkey(bins, chr, start, end)

ov <- foverlaps(
  dt, bins,
  by.x = c("chr","start","end"),
  by.y = c("chr","start","end"),
  type = "any",
  nomatch = 0
)

if (nrow(ov) == 0) stop("No overlap between coverage and active HOR bins.")

ov[, `:=`(
  ov_start = pmax(start, i.start),
  ov_end = pmin(end, i.end)
)]
ov[, ov_len := ov_end - ov_start]
ov <- ov[ov_len > 0]

bin_dt <- ov[, .(
  covered_bp = sum(ov_len),
  fire_bp = sum(fire * ov_len, na.rm = TRUE),
  linker_bp = sum(linker * ov_len, na.rm = TRUE),
  nuc_bp = sum(nucleosome * ov_len, na.rm = TRUE)
), by = .(
  sample, chr, region_id, start, end,
  cdr_start, cdr_end, in_cdr, bin_width
)]

bin_dt[, total_bp := fire_bp + linker_bp + nuc_bp]
bin_dt <- bin_dt[total_bp > 0]

bin_dt[, `:=`(pfire=fire_bp/total_bp, plinker=linker_bp/total_bp, pnuc=nuc_bp/total_bp)]

setorder(bin_dt, chr, region_id, sample, start)

if (smooth_k > 1) {
  bin_dt[, pfire_smooth := zoo::rollmean(pfire, k=smooth_k, fill=NA, align="center"), by=.(chr,region_id,sample)]
} else {
  bin_dt[, pfire_smooth := pfire]
}

bin_dt <- bin_dt[is.finite(pfire_smooth)]

## ---------- functions ----------
winsorize_vec <- function(x, trim_fraction) {
  if (is.null(trim_fraction) || trim_fraction <= 0) return(x)
  fx <- x[is.finite(x)]
  if (length(fx) == 0) return(x)
  q <- quantile(fx, probs = c(trim_fraction, 1 - trim_fraction), na.rm = TRUE)
  pmin(pmax(x, q[1]), q[2])
}

get_cell_group <- function(x) {
  fifelse(grepl("PBMC", x, ignore.case = TRUE), "PBMC",
  fifelse(grepl("iPSC", x, ignore.case = TRUE), "iPSC",
  fifelse(grepl("NPC", x, ignore.case = TRUE), "NPC",
  fifelse(grepl("Monocyte", x, ignore.case = TRUE), "Monocyte",
  fifelse(grepl("Macrophage", x, ignore.case = TRUE), "Macrophage", "Other")))))
}

calc_one_fas <- function(x,out_cutoff,in_cutoff){
  bg <- copy(x[in_cdr==FALSE]); cdr <- copy(x[in_cdr==TRUE])
  empty <- data.table(FAS=NA_real_, log2_FAS=NA_real_, FIRE_enrichment=NA_real_, background_FIRE=NA_real_, CDR_FIRE=NA_real_, CDR_covered_bp=NA_real_, n_active_bins=nrow(x), n_bg_bins=nrow(bg), n_CDR_bins=nrow(cdr))
  if(nrow(bg)==0 || nrow(cdr)==0) return(empty)

  bg[, pfire_w := winsorize_vec(pfire_smooth,out_cutoff)]
  cdr[, pfire_w := winsorize_vec(pfire_smooth,in_cutoff)]

  background <- mean(bg$pfire_w,na.rm=TRUE); cdr_mean <- mean(cdr$pfire_w,na.rm=TRUE); cdr_bp <- sum(cdr$covered_bp,na.rm=TRUE)
  if(!is.finite(background) || !is.finite(cdr_mean) || !is.finite(cdr_bp) || cdr_bp<=0) return(empty)

  eps <- 1e-6
  cdr[, delta := (pfire_w-background)/(background+eps)]
  fas <- sum(cdr$delta*cdr$covered_bp,na.rm=TRUE)

  data.table(FAS=round(fas,4), log2_FAS=round(sign(fas)*log2(abs(fas)+1),4), FIRE_enrichment=round(cdr_mean/(background+eps),4), background_FIRE=round(background,4), CDR_FIRE=round(cdr_mean,4), CDR_covered_bp=round(cdr_bp,4), n_active_bins=nrow(x), n_bg_bins=nrow(bg), n_CDR_bins=nrow(cdr))
}

## ---------- call FAS ----------
res <- bin_dt[, calc_one_fas(.SD, out_CDR_cutoff, in_CDR_cutoff),
              by = .(sample, chr, region_id, cdr_start, cdr_end)]
if (nrow(res) == 0) stop("No FAS results were produced")

res[, `:=`(
  cell = get_cell_group(sample),
  asm = sub("_.*", "", sample),
  bin_size = bin_size,
  smooth_k = smooth_k,
  cdr_expansion = cdr_expansion,
  out_CDR_cutoff = out_CDR_cutoff,
  in_CDR_cutoff = in_CDR_cutoff
)]

setcolorder(res, c(
  "sample","asm","cell","chr","region_id","cdr_start","cdr_end",
  "FAS","log2_FAS","FIRE_enrichment",
  "background_FIRE","CDR_FIRE","CDR_covered_bp",
  "n_active_bins","n_bg_bins","n_CDR_bins",
  "bin_size","smooth_k","cdr_expansion","out_CDR_cutoff","in_CDR_cutoff"
), skip_absent = TRUE)

num_cols <- names(res)[sapply(res, is.numeric)]
res[, (num_cols) := lapply(.SD, function(x) round(x, 4)), .SDcols = num_cols]

fwrite(bin_dt, bin_out, sep = "\t")
fwrite(res, out_file, sep = "\t")

## ---------- statistics ----------
metrics <- c("log2_FAS", "FIRE_enrichment")

stat_long <- melt(
  res,
  id.vars = c("sample","asm","cell","chr","region_id"),
  measure.vars = metrics,
  variable.name = "metric",
  value.name = "value"
)[is.finite(value)]

stat_res <- rbindlist(lapply(metrics, function(met) {
  wide <- dcast(
    stat_long[metric == met],
    asm + chr + region_id ~ cell,
    value.var = "value"
  )
  if (!all(c("iPSC", "NPC") %in% names(wide))) return(NULL)
  wide <- wide[is.finite(iPSC) & is.finite(NPC)]
  if (nrow(wide) <= 1) return(NULL)

  wt <- wilcox.test(wide$NPC, wide$iPSC, paired = TRUE)
  wide[, `:=`(
    diff = NPC - iPSC,
    percent_change = (NPC - iPSC) / (abs(iPSC) + 1e-6) * 100
  )]

  data.table(
    metric = met,
    comparison = "NPC_vs_iPSC",
    n_pair = nrow(wide),
    mean_iPSC = round(mean(wide$iPSC), 4),
    mean_NPC = round(mean(wide$NPC), 4),
    mean_difference = round(mean(wide$diff), 4),
    median_difference = round(median(wide$diff), 4),
    mean_percent_change = round(mean(wide$percent_change), 4),
    median_percent_change = round(median(wide$percent_change), 4),
    p_value = wt$p.value
  )
}), fill = TRUE)

if (nrow(stat_res) > 0) {
  stat_res[, FDR := round(p.adjust(p_value, method = "BH"), 4)]
}
fwrite(stat_res, stat_out, sep = "\t")

## ---------- summary plot ----------
cell_cols <- c(
  PBMC = "#4E79A7", iPSC = "#E15759", NPC = "#59A14F",
  Monocyte = "#F28E2B", Macrophage = "#B07AA1", Other = "grey60"
)

p_stat <- ggplot(stat_long, aes(cell, value, fill = cell)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.65) +
  geom_jitter(aes(color = cell), width = 0.15, size = 1.1, alpha = 0.65) +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = cell_cols) +
  scale_color_manual(values = cell_cols) +
  labs(x = parameter_summary, y = "Value", title = "FAS comparison across cell types") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(face = "bold")
  )

ggsave(stat_plot, p_stat, width = 8, height = 4.5)

## ---------- per-chr QC plots ----------
plot_one_chr <- function(ch) {
  plot_dt <- bin_dt[chr == ch]
  if (nrow(plot_dt) == 0) return(NULL)

  plot_dt[, cell := get_cell_group(sample)]
  cen_row <- cen[chr == ch][1]

  label_dt <- res[chr == ch, .(
    sample,
    label = paste0("log2(FAS) = ", log2_FAS, "\nEnrichment = ", FIRE_enrichment)
  )]

  bg_dt <- res[chr == ch, .(
    sample,
    background_FIRE,
    cell = get_cell_group(sample)
  )]

  ggplot() +
    geom_rect(
      data = cen_row,
      aes(xmin = cdr_start_exp, xmax = cdr_end_exp, ymin = -Inf, ymax = Inf),
      fill = "grey85", alpha = 0.65, inherit.aes = FALSE
    ) +
    geom_line(
      data = plot_dt,
      aes(start, pfire_smooth, color = cell, group = sample),
      linewidth = 0.35, alpha = 0.85
    ) +
    geom_hline(
      data = bg_dt,
      aes(yintercept = background_FIRE, color = cell),
      linetype = "dashed", linewidth = 0.35
    ) +
    geom_text(
      data = label_dt,
      aes(x = Inf, y = Inf, label = label),
      hjust = 1.1, vjust = 1.4, size = 2.8,
      inherit.aes = FALSE
    ) +
    facet_wrap(~sample, ncol = 1) +
    scale_x_continuous(labels = label_number(scale = 1e-6, suffix = " Mb")) +
    scale_color_manual(values = cell_cols) +
    labs(
      x = "Position",
      y = "Smoothed FIRE fraction",
      title = paste0(ch, ": active HOR / CDR")
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      strip.text = element_text(size = 7.5, face = "bold"),
      plot.title = element_text(face = "bold")
    )
}

pdf(pdf_file, width = 8, height = 6)
for (ch in unique(cen$chr)) {
  p <- plot_one_chr(ch)
  if (!is.null(p)) print(p)
}
dev.off()

message("Done.")
message("FAS result: ", out_file)
message("Bin summary: ", bin_out)
message("Stats: ", stat_out)
message("QC plots: ", plotdir)

