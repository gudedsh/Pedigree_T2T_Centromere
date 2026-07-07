#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(optparse)
  library(zoo)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sample"),
  make_option("--meth_file"),
  make_option("--active_hor"),
  make_option("--out_dir"),
  make_option("--bin_gap", type = "integer", default = 50000),
  make_option("--min_bins", type = "integer", default = 5),
  make_option("--cutoff_frac", type = "double", default = 0.15),
  make_option("--smooth_k", type = "integer", default = 10),
  make_option("--depth_cutoff", type = "double", default = 0)
)))

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# read methylation data
# -------------------------
dt <- fread(opt$meth_file)

dt <- dt[, .(
  chr = chr,
  start = as.integer(start),
  end = as.integer(end),
  sample = opt$sample,
  meth = as.numeric(meth),
  depth = as.numeric(depth)
)]

dt <- dt[depth >= opt$depth_cutoff & is.finite(meth)]
setorder(dt, chr, start)

# smooth only for QC visualization
dt[, meth_smooth := zoo::rollmean(meth, k = opt$smooth_k, fill = NA, align = "center"), by = chr]

# -------------------------
# read active HOR
# -------------------------
active <- fread(opt$active_hor, header = FALSE)
active <- active[, .(chr = V1, start = as.integer(V2), end = as.integer(V3))]
active[, len := end - start]

# -------------------------
# call CDR
# -------------------------
call_one_chr <- function(sub) {
  if (nrow(sub) == 0) return(NULL)

  med <- median(sub$meth, na.rm = TRUE)
  cutoff <- med * opt$cutoff_frac

  cdr_bin <- sub[meth < med - cutoff]
  if (nrow(cdr_bin) == 0) return(NULL)

  setorder(cdr_bin, start, end)
  cdr_bin[, gap := start - shift(end)]
  cdr_bin[, block := cumsum(is.na(gap) | gap > opt$bin_gap)]

  cdr_bin[, .(
    start = min(start),
    end = max(end),
    nbin = .N,
    mean_meth = mean(meth, na.rm = TRUE),
    median_active_hor_meth = med,
    cutoff = cutoff
  ), by = .(chr, block)][nbin >= opt$min_bins]
}

cdr <- rbindlist(
  lapply(unique(dt$chr), function(ch) call_one_chr(dt[chr == ch])),
  use.names = TRUE,
  fill = TRUE
)

if (is.null(cdr) || nrow(cdr) == 0) {
  cdr <- data.table(
    chr = character(),
    block = integer(),
    start = integer(),
    end = integer(),
    nbin = integer(),
    mean_meth = numeric(),
    median_active_hor_meth = numeric(),
    cutoff = numeric()
  )
}

cdr[, sample := opt$sample]

# -------------------------
# output CDR
# -------------------------
fwrite(
  cdr[, .(chr, start, end, sample, nbin, mean_meth, median_active_hor_meth, cutoff)],
  file.path(opt$out_dir, paste0(opt$sample, ".identified_CDR.bed")),
  sep = "\t"
)

fwrite(
  cdr[, .(
    n_cdr = .N,
    total_cdr_bp = sum(end - start),
    mean_cdr_bp = mean(end - start)
  ), by = sample],
  file.path(opt$out_dir, paste0(opt$sample, ".identified_CDR.summary.tsv")),
  sep = "\t"
)

# -------------------------
# QC plot: all chromosomes into one PDF
# -------------------------
pdf(file.path(opt$out_dir, paste0(opt$sample, ".CDR_QC.pdf")), width = 8, height = 4)

for (ch in unique(dt$chr)) {
  sub <- copy(dt[chr == ch])
  sub[, in_cdr := FALSE]

  cdr_ch <- cdr[chr == ch]

  if (nrow(cdr_ch) > 0) {
    for (i in seq_len(nrow(cdr_ch))) {
      sub[start >= cdr_ch$start[i] & end <= cdr_ch$end[i], in_cdr := TRUE]
    }
  }

  stat <- sub[, .(
    cdr_mean = mean(meth[in_cdr], na.rm = TRUE),
    noncdr_mean = mean(meth[!in_cdr], na.rm = TRUE)
  )]

  cdr_len_kb <- if (nrow(cdr_ch) > 0) sum(cdr_ch$end - cdr_ch$start) / 1000 else 0

  lab <- sprintf(
    "Non-CDR: %.1f; CDR: %.1f; Δ=%.1f; CDR len: %.0f kb",
    stat$noncdr_mean,
    stat$cdr_mean,
    stat$noncdr_mean - stat$cdr_mean,
    cdr_len_kb
  )

  p <- ggplot(sub, aes(start, meth)) +
    geom_line(linewidth = 0.25, alpha = 0.5) +
    geom_line(aes(y = meth_smooth), linewidth = 0.45, na.rm = TRUE) +
    geom_rect(
      data = cdr_ch,
      aes(xmin = start, xmax = end, ymin = 100, ymax = 103),
      inherit.aes = FALSE,
      fill = "grey30"
    ) +
    annotate(
      "text",
      x = min(sub$start),
      y = 5,
      label = lab,
      hjust = 0,
      size = 3
    ) +
    scale_x_continuous(labels = label_number(scale = 1e-6, suffix = " Mb")) +
    scale_y_continuous(limits = c(0, 103), breaks = c(0, 50, 100)) +
    labs(
      title = paste0(opt$sample, " | ", ch, " active HOR"),
      x = NULL,
      y = "Methylation (%)"
    ) +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank())

  print(p)
}

dev.off()
