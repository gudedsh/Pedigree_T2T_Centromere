#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(zoo)
  library(scales)
})

# -----------------------------
# Parse command-line arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx == length(args)) stop("Missing value for ", flag)
  args[idx + 1]
}

meth_file <- get_arg("--meth")
cen_file  <- get_arg("--cen")
out_file  <- get_arg("--out", "MARS_results.tsv")
plotdir   <- get_arg("--plotdir", "MARS_QC_plots")

smooth_k <- as.integer(get_arg("--smooth_k", "10"))
out_CDR_cutoff <- as.numeric(get_arg("--out_CDR_cutoff", "0.10"))
in_CDR_cutoff  <- as.numeric(get_arg("--in_CDR_cutoff", "0.05"))
bin_size <- as.integer(get_arg("--bin_size", "1000"))


if (is.null(meth_file) || is.null(cen_file)) {
  stop("Usage: Rscript call_mars.r --meth combined.tsv --cen active_hor_cdr.bed --out MARS_results.tsv --plotdir plots")
}

dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)

message("Reading methylation file: ", meth_file)
dt <- fread(meth_file)

required_cols <- c("chr", "start", "end", "nCG", "depth", "meth", "sample")

if (!all(required_cols %in% colnames(dt))) {
  stop("Methylation file must contain columns: ", paste(required_cols, collapse = ", "))
}

dt[, meth := as.numeric(meth)]
dt <- dt[is.finite(meth)]
setorder(dt, chr, sample, start)
dt
message("Smoothing methylation with window k = ", smooth_k)

if (smooth_k > 1) {
  dt[, meth_smooth := zoo::rollmean(meth, k = smooth_k,fill = NA,align = "center"), by = .(chr, sample)]
} else {
  dt[, meth_smooth := meth]
}

dt <- dt[is.finite(meth_smooth)]

message("Reading CDR annotation: ", cen_file)
cen <- fread(cen_file, header = FALSE)
cen
if (ncol(cen) < 5) {
  stop("CEN file must have at least 5 columns: chr active_hor_start active_hor_end cdr_start cdr_end")
}

setnames(cen, 1:5, c("chr", "active_start", "active_end", "cdr_start", "cdr_end"))
cen <- cen[, .(
  chr = as.character(chr),
  active_start = as.integer(active_start),
  active_end = as.integer(active_end),
  cdr_start = as.integer(cdr_start),
  cdr_end = as.integer(cdr_end)
)]

cen[, cdr_len := cdr_end - cdr_start]

winsorize_vec <- function(x, trim_fraction) {
  x2 <- x
  if (is.null(trim_fraction) || trim_fraction <= 0) {return(x2)}
  finite_x <- x[is.finite(x)]
  if (length(finite_x) == 0) {return(x2)}
  low <- quantile(finite_x, probs = trim_fraction, na.rm = TRUE)
  high <- quantile(finite_x, probs = 1 - trim_fraction, na.rm = TRUE)
  x2[is.finite(x2) & x2 < low] <- low
  x2[is.finite(x2) & x2 > high] <- high
  return(x2)
}


calc_one_mars <- function(x,cdr_start,cdr_end,active_start,active_end,out_cutoff,in_cutoff) {
  x <- x[start >= active_start & end <= active_end]
  empty_result <- function(
      background = NA_real_,
      bg_methylation = NA_real_,
      cdr_methylation = NA_real_,
      bg_methylation_raw = NA_real_,
      cdr_methylation_raw = NA_real_,
      n_active_bins = 0,
      n_bg_bins_raw = 0,
      n_bg_bins_used = 0,
      n_cdr_bins_raw = 0,
      n_cdr_bins_used = 0
  ) {
    data.table(
      background = background,
      bg_methylation = bg_methylation,
      cdr_methylation = cdr_methylation,
      bg_methylation_raw = bg_methylation_raw,
      cdr_methylation_raw = cdr_methylation_raw,
      methylation_difference = NA_real_,
      n_active_bins = n_active_bins,
      n_bg_bins_raw = n_bg_bins_raw,
      n_bg_bins_used = n_bg_bins_used,
      n_cdr_bins_raw = n_cdr_bins_raw,
      n_cdr_bins_used = n_cdr_bins_used,
      MARS = NA_real_,
      log2_MARS = NA_real_,
      mean_delta = NA_real_
    )
  }
  
  if (nrow(x) == 0) { return(empty_result())}
  
  bg <- x[end <= cdr_start | start >= cdr_end]
  cdr <- x[start >= cdr_start & end <= cdr_end]
  
  n_bg_raw <- nrow(bg)
  n_cdr_raw <- nrow(cdr)
  
  if (n_bg_raw == 0 || n_cdr_raw == 0) {
    return(empty_result(
      n_active_bins = nrow(x),
      n_bg_bins_raw = n_bg_raw,
      n_cdr_bins_raw = n_cdr_raw
    ))
  }
  
  bg_methylation_raw <- mean(bg$meth_smooth, na.rm = TRUE)
  cdr_methylation_raw <- mean(cdr$meth_smooth, na.rm = TRUE)
  
  bg_used <- copy(bg)
  bg_used[, meth_smooth_winsor := winsorize_vec(meth_smooth, out_cutoff)]
  
  cdr_used <- copy(cdr)
  cdr_used[, meth_smooth_winsor := winsorize_vec(meth_smooth, in_cutoff)]
  
  background <- mean(bg_used$meth_smooth_winsor, na.rm = TRUE)
  bg_methylation <- background
  cdr_methylation <- mean(cdr_used$meth_smooth_winsor, na.rm = TRUE)
  
  cdr_used[, delta := background - meth_smooth_winsor]
  
  if (
    nrow(cdr_used) == 0 ||
    !is.finite(background) ||
    !is.finite(cdr_methylation)
  ) {
    return(empty_result(
      background = background,
      bg_methylation = bg_methylation,
      cdr_methylation = cdr_methylation,
      bg_methylation_raw = bg_methylation_raw,
      cdr_methylation_raw = cdr_methylation_raw,
      n_active_bins = nrow(x),
      n_bg_bins_raw = n_bg_raw,
      n_bg_bins_used = nrow(bg_used),
      n_cdr_bins_raw = n_cdr_raw,
      n_cdr_bins_used = nrow(cdr_used)
    ))
  }
  
  mars <- round(sum(cdr_used$delta, na.rm = TRUE), 1)
  
  data.table(
    background = round(background, 4),
    bg_methylation = round(bg_methylation, 4),
    cdr_methylation = round(cdr_methylation, 4),
    bg_methylation_raw = round(bg_methylation_raw, 4),
    cdr_methylation_raw = round(cdr_methylation_raw, 4),
    methylation_difference = round(bg_methylation - cdr_methylation, 4),
    n_active_bins = nrow(x),
    n_bg_bins_raw = n_bg_raw,
    n_bg_bins_used = nrow(bg_used),
    n_cdr_bins_raw = n_cdr_raw,
    n_cdr_bins_used = nrow(cdr_used),
    MARS = mars,
    log2_MARS = round(log2(mars + 1e-6), 3),
    mean_delta = round(mean(cdr_used$delta, na.rm = TRUE), 3)
  )
}

get_cell_group <- function(x){
  fifelse(
    grepl("PBMC", x, ignore.case = TRUE), "PBMC",
  fifelse(
    grepl("iPSC", x, ignore.case = TRUE), "iPSC",
  fifelse(
    grepl("NPC", x, ignore.case = TRUE), "NPC",
  fifelse(
    grepl("Monocyte", x, ignore.case = TRUE), "Monocyte",
  fifelse(
    grepl("Macrophage", x, ignore.case = TRUE), "Macrophage",
    "Other"
  )))))
}

sample_order_by_celltype <- function(samples){
  dt <- data.table(sample = unique(samples))
  dt[, cell_group := get_cell_group(sample)]
  dt[, cell_rank := match(cell_group, c("PBMC","iPSC","NPC","Monocyte","Macrophage","Other"))]
  dt[order(cell_rank, sample)]$sample
}


plot_one_chr <- function(plot_dt, cen_row, outfile) {
  if (nrow(plot_dt) == 0) return(NULL)
  p <- ggplot() +
    geom_rect(data = cen_row,aes(xmin = cdr_start, xmax = cdr_end, ymin = -Inf, ymax = Inf),fill = "grey80",alpha = 0.6,inherit.aes = FALSE) +
    geom_line(data = plot_dt,aes(x = start, y = meth_smooth, color = sample, group = sample),linewidth = 0.35) +
    scale_x_continuous(labels = label_number(scale = 1e-6)) +
 #   coord_cartesian(ylim = c(0, 100)) +
    labs(x = "Genomic position (Mb)",y = "DNA methylation (%)",title = paste0(cen_row$chr, ": active HOR / CDR"),color = "Sample") +
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      legend.position = "right"
    )
  
  ggsave(outfile, p, width = 8, height = 4)
}

plot_one_chr_facet_mars <- function(plot_dt, cen_row, mars_res_chr, outfile, out_cutoff, in_cutoff) {
  if (nrow(plot_dt) == 0) return(NULL)
  plot_dt <- copy(plot_dt)
  plot_dt[, cell := get_cell_group(sample)]
  sample_order <- sample_order_by_celltype(plot_dt$sample)
  plot_dt[, sample := factor(sample, levels = sample_order)]
  delta_list <- list()
  for (sid in sample_order) {
    x <- plot_dt[sample == sid]
    bg <- x[end <= cen_row$cdr_start | start >= cen_row$cdr_end]
    cdr <- x[start >= cen_row$cdr_start & end <= cen_row$cdr_end]
    if (nrow(bg) == 0 || nrow(cdr) == 0) next
    bg_used <- copy(bg)
    bg_used[, meth_smooth_winsor := winsorize_vec(meth_smooth, out_cutoff)]
    background <- mean(bg_used$meth_smooth_winsor, na.rm = TRUE)
    cdr_used <- copy(cdr)
    cdr_used[, meth_smooth_winsor := winsorize_vec(meth_smooth, in_cutoff)]
    cdr_used[, background := background]
    cdr_used[, delta := background - meth_smooth_winsor]
    cdr_used[, cell := get_cell_group(as.character(sid))]
    delta_list[[as.character(sid)]] <- cdr_used
  }
  
  delta_dt <- rbindlist(delta_list, fill = TRUE)
  delta_dt[, sample := factor(sample, levels = sample_order)]
  
  label_dt <- mars_res_chr[, .(sample, label = paste0("log2(MARS) = ", round(log2_MARS, 3)))]
  label_dt[, sample := factor(sample, levels = sample_order)]
  
  bg_dt <- unique(delta_dt[, .(sample, background, cell)])
  bg_dt[, sample := factor(sample, levels = sample_order)]
  
  cell_cols2 <- c(PBMC = "#4E79A7",iPSC = "#E15759",NPC = "#59A14F",Monocyte = "#F28E2B",Macrophage = "#B07AA1",Other = "grey60")
  
  p <- ggplot() +
    geom_rect(data = cen_row,aes(xmin = cdr_start, xmax = cdr_end, ymin = -Inf, ymax = Inf),fill = "grey85",alpha = 0.65,inherit.aes = FALSE) +
    geom_point(data = plot_dt,aes(x = start, y = meth_smooth),size = 0.1,color = "darkgrey",alpha = 0.7) +
    geom_hline(data = bg_dt,aes(yintercept = background, color = cell),linetype = "dashed",linewidth = 0.35) +
    geom_segment(
      data = delta_dt,
      aes(x = start,xend = start,y = meth_smooth_winsor,yend = background,color = cell),linewidth = 0.25,alpha = 0.55) +
    geom_text(data = label_dt,aes(x = Inf, y = 25, label = label),hjust = 1.5,size = 3,inherit.aes = FALSE) +
    facet_wrap(~sample, ncol = 1) +
    scale_x_continuous(labels = label_number(scale = 1e-6)) +
    scale_color_manual(values = cell_cols2) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(x = "Genomic position (Mb)",y = "DNA methylation (%)",title = paste0(cen_row$chr, ": active HOR / CDR")) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      strip.text = element_text(size = 8, face = "bold"),
      plot.title = element_text(face = "bold")
    )
  
  ggsave(outfile,p,width = 8,height = max(4, 1.25 * uniqueN(plot_dt$sample)))
}

message("Calculating MARS...")

res_list <- list()

for (i in seq_len(nrow(cen))) {
  cr <- cen[i]
  chr_dt <- dt[chr == cr$chr & start >= cr$active_start & end <= cr$active_end]
  if (nrow(chr_dt) == 0) {
    warning("No methylation bins found for ", cr$chr)
    next
  }

  # Keep only bins detected in all samples for this active HOR
  nsamp <- uniqueN(chr_dt$sample)
  chr_dt_common <- chr_dt[, if (uniqueN(sample) == nsamp) .SD, by = .(chr, start, end)]

  for (sid in unique(chr_dt_common$sample)) {
    x <- chr_dt_common[sample == sid]
    one <- calc_one_mars(
      x = x,
      cdr_start = cr$cdr_start,
      cdr_end = cr$cdr_end,
      active_start = cr$active_start,
      active_end = cr$active_end,
      out_cutoff = out_CDR_cutoff,
      in_cutoff = in_CDR_cutoff
    )
    
    res_list[[length(res_list) + 1]] <- cbind(
      data.table(
        sample = sid,
        chr = cr$chr,
        active_start = cr$active_start,
        active_end = cr$active_end,
        cdr_start = cr$cdr_start,
        cdr_end = cr$cdr_end,
        cdr_len = cr$cdr_len,
        bin_size = bin_size,
        smooth_k = smooth_k,
        out_CDR_cutoff = out_CDR_cutoff,
        in_CDR_cutoff = in_CDR_cutoff,
        common_sample_n = nsamp
      ),
      one
    )
  }

  plot_file1 <- file.path(plotdir, paste0(cr$chr, "_MARS_QC1_all_samples.pdf"))
  plot_one_chr(chr_dt_common, cr, plot_file1)

  mars_res_chr <- rbindlist(res_list, fill = TRUE)[chr == cr$chr]
  plot_file2 <- file.path(plotdir, paste0(cr$chr, "_MARS_QC2_sample_facet.pdf"))
  
plot_one_chr_facet_mars(
  plot_dt = chr_dt_common,
  cen_row = cr,
  mars_res_chr = mars_res_chr,
  outfile = plot_file2,
  out_cutoff = out_CDR_cutoff,
  in_cutoff = in_CDR_cutoff
)

}

res <- rbindlist(res_list, fill = TRUE)

setcolorder(
  res,
  c(
    "sample",

    "cell",
    "ind",
    "clo",
    "tech",
    "hap",
    "asm",

    "chr",
    "active_start",
    "active_end",

    "cdr_start",
    "cdr_end",
    "cdr_len",

    "MARS",
    "log2_MARS",

    "mean_delta",
    "methylation_difference",

    "bg_methylation",
    "cdr_methylation",

    "bg_methylation_raw",
    "cdr_methylation_raw",

    "n_active_bins",
    "n_bg_bins_raw",
    "n_cdr_bins_raw",

    "bin_size",
    "smooth_k",

    "out_CDR_cutoff",
    "in_CDR_cutoff",

    "common_sample_n"
  ),
  skip_absent = TRUE
)

fwrite(res, out_file, sep = "\t")

message("Done.")
message("MARS result written to: ", out_file)
message("QC plots written to: ", plotdir)

