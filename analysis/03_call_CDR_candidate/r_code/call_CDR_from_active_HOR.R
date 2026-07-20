#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option("--sample", type = "character", help = "Unique sample name"),
  make_option("--meth_file", type = "character", help = "Binned methylation TSV"),
  make_option("--active_hor", type = "character", help = "Active-HOR BED file"),
  make_option("--out_dir", type = "character", help = "Output directory"),
  make_option("--bin_gap", type = "integer", default = 50000L,
              help = "Maximum gap between candidate bins [default: %default]"),
  make_option("--min_bins", type = "integer", default = 5L,
              help = "Minimum bins per CDR [default: %default]"),
  make_option("--cutoff_frac", type = "double", default = 0.15,
              help = "Fraction below chromosome median [default: %default]"),
  make_option("--smooth_k", type = "integer", default = 10L,
              help = "Bins in centered QC rolling mean [default: %default]"),
  make_option("--depth_cutoff", type = "double", default = 0,
              help = "Minimum mean bin depth [default: %default]")
)

parser <- OptionParser(
  usage = "%prog --sample NAME --meth_file FILE --active_hor FILE --out_dir DIR [options]",
  description = "Call candidate centromeric dip regions from binned methylation data.",
  option_list = option_list
)
opt <- parse_args(parser)

abort <- function(...) {
  message("ERROR: ", sprintf(...))
  quit(save = "no", status = 1L)
}

required_options <- c("sample", "meth_file", "active_hor", "out_dir")
missing_options <- required_options[vapply(
  required_options,
  function(x) is.null(opt[[x]]) || !nzchar(opt[[x]]),
  logical(1)
)]
if (length(missing_options)) {
  print_help(parser)
  abort("Missing required option(s): %s", paste(missing_options, collapse = ", "))
}

if (!file.exists(opt$meth_file)) abort("Methylation file not found: %s", opt$meth_file)
if (!file.exists(opt$active_hor)) abort("Active-HOR file not found: %s", opt$active_hor)
if (is.na(opt$bin_gap) || opt$bin_gap <= 0) abort("--bin_gap must be a positive integer")
if (is.na(opt$min_bins) || opt$min_bins <= 0) abort("--min_bins must be a positive integer")
if (is.na(opt$smooth_k) || opt$smooth_k <= 0) abort("--smooth_k must be a positive integer")
if (!is.finite(opt$cutoff_frac) || opt$cutoff_frac < 0 || opt$cutoff_frac > 1) {
  abort("--cutoff_frac must be between 0 and 1")
}
if (!is.finite(opt$depth_cutoff) || opt$depth_cutoff < 0) {
  abort("--depth_cutoff must be non-negative")
}

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

# Read and validate binned methylation data.
dt <- tryCatch(
  fread(opt$meth_file),
  error = function(e) abort("Cannot read methylation file: %s", conditionMessage(e))
)
required_columns <- c("chr", "start", "end", "meth", "depth")
missing_columns <- setdiff(required_columns, names(dt))
if (length(missing_columns)) {
  abort("Methylation file is missing column(s): %s", paste(missing_columns, collapse = ", "))
}

dt <- dt[, .(
  chr = as.character(chr),
  start = suppressWarnings(as.integer(start)),
  end = suppressWarnings(as.integer(end)),
  meth = suppressWarnings(as.numeric(meth)),
  depth = suppressWarnings(as.numeric(depth))
)]

invalid_coordinates <- dt[is.na(chr) | !nzchar(chr) | is.na(start) | is.na(end) | start < 0 | end <= start]
if (nrow(invalid_coordinates)) {
  abort("Methylation file contains %d row(s) with invalid genomic coordinates", nrow(invalid_coordinates))
}

non_numeric <- dt[!is.finite(meth) | !is.finite(depth)]
if (nrow(non_numeric)) {
  message(sprintf("WARNING: removing %d row(s) with non-numeric methylation or depth", nrow(non_numeric)))
}
dt <- dt[is.finite(meth) & is.finite(depth) & depth >= opt$depth_cutoff]
if (!nrow(dt)) abort("No methylation bins remain after input and depth filtering")
setorder(dt, chr, start, end)

# Read the first three BED columns. Comment and blank lines are removed before
# parsing so files with a '#chr' header are handled correctly.
active_lines <- readLines(opt$active_hor, warn = FALSE)
active_lines <- active_lines[!grepl("^[[:space:]]*(#|$)", active_lines)]
if (!length(active_lines)) abort("Active-HOR file contains no intervals")

active <- tryCatch(
  fread(text = paste(active_lines, collapse = "\n"), header = FALSE, select = 1:3),
  error = function(e) abort("Cannot parse active-HOR file: %s", conditionMessage(e))
)
setnames(active, c("chr", "start", "end"))
active[, `:=`(
  chr = as.character(chr),
  start = suppressWarnings(as.integer(start)),
  end = suppressWarnings(as.integer(end))
)]
if (active[is.na(chr) | is.na(start) | is.na(end) | start < 0 | end <= start, .N] > 0L) {
  abort("Active-HOR file contains invalid BED intervals")
}
setorder(active, chr, start, end)

# Assign each bin to an active-HOR interval. This both validates the shell-stage
# extraction and prevents smoothing/calling across distinct active-HOR arrays.
dt[, bin_id := .I]
overlaps <- active[dt, on = .(chr, start <= start, end >= end), nomatch = 0L,
                   .(bin_id = i.bin_id, active_start = x.start, active_end = x.end)]
if (!nrow(overlaps)) abort("No methylation bins overlap the supplied active-HOR intervals")
if (overlaps[, anyDuplicated(bin_id)] > 0L) {
  abort("Overlapping active-HOR intervals assign at least one bin more than once")
}
dt <- overlaps[dt, on = "bin_id", nomatch = 0L]
setcolorder(dt, c("chr", "start", "end", "meth", "depth", "active_start", "active_end"))
setorder(dt, chr, active_start, start)

# Smoothing is used only for QC. Keep it within each active-HOR interval and
# break it at large gaps so the visual trace never bridges missing regions.
dt[, qc_segment := cumsum(
  is.na(shift(end)) | start - shift(end) > opt$bin_gap
), by = .(chr, active_start, active_end)]
dt[, meth_smooth := frollmean(
  meth,
  n = opt$smooth_k,
  align = "center",
  fill = NA_real_
), by = .(chr, active_start, active_end, qc_segment)]

call_one_region <- function(sub) {
  if (!nrow(sub)) return(NULL)

  region_median <- median(sub$meth, na.rm = TRUE)
  threshold <- region_median * (1 - opt$cutoff_frac)
  candidate_bins <- copy(sub[meth < threshold])
  if (!nrow(candidate_bins)) return(NULL)

  setorder(candidate_bins, start, end)
  candidate_bins[, gap := start - shift(end)]
  candidate_bins[, block := cumsum(is.na(gap) | gap > opt$bin_gap)]

  candidate_bins[, .(
    start = min(start),
    end = max(end),
    nbin = .N,
    mean_meth = mean(meth),
    median_active_hor_meth = region_median,
    cutoff = region_median * opt$cutoff_frac,
    threshold = threshold
  ), by = block][nbin >= opt$min_bins]
}

cdr <- dt[, call_one_region(.SD), by = .(chr, active_start, active_end)]

output_columns <- c(
  "chr", "start", "end", "sample", "nbin", "mean_meth",
  "median_active_hor_meth", "cutoff", "threshold"
)
if (!nrow(cdr)) {
  cdr <- data.table(
    chr = character(), start = integer(), end = integer(),
    nbin = integer(), mean_meth = numeric(),
    median_active_hor_meth = numeric(), cutoff = numeric(), threshold = numeric()
  )
}
cdr[, sample := opt$sample]

cdr_file <- file.path(opt$out_dir, paste0(opt$sample, ".identified_CDR.bed"))
summary_file <- file.path(opt$out_dir, paste0(opt$sample, ".identified_CDR.summary.tsv"))
qc_file <- file.path(opt$out_dir, paste0(opt$sample, ".CDR_QC.pdf"))

fwrite(cdr[, ..output_columns], cdr_file, sep = "\t", na = "NA",quote=F)

summary_dt <- data.table(
  sample = opt$sample,
  n_cdr = nrow(cdr),
  total_cdr_bp = if (nrow(cdr)) sum(cdr$end - cdr$start) else 0,
  mean_cdr_bp = if (nrow(cdr)) mean(cdr$end - cdr$start) else NA_real_
)
fwrite(summary_dt, summary_file, sep = "\t", na = "NA")

# Produce one QC page per chromosome.
grDevices::pdf(qc_file, width = 8, height = 4, onefile = TRUE)
on.exit(grDevices::dev.off(), add = TRUE)

for (ch in unique(dt$chr)) {
  sub <- copy(dt[chr == ch])
  sub[, in_cdr := FALSE]
  cdr_ch <- cdr[chr == ch]

  if (nrow(cdr_ch)) {
    for (i in seq_len(nrow(cdr_ch))) {
      sub[start >= cdr_ch$start[i] & end <= cdr_ch$end[i], in_cdr := TRUE]
    }
  }

  safe_mean <- function(x) if (length(x)) mean(x, na.rm = TRUE) else NA_real_
  cdr_mean <- safe_mean(sub[in_cdr == TRUE, meth])
  noncdr_mean <- safe_mean(sub[in_cdr == FALSE, meth])
  cdr_len_kb <- if (nrow(cdr_ch)) sum(cdr_ch$end - cdr_ch$start) / 1000 else 0
  fmt <- function(x) if (is.finite(x)) sprintf("%.1f", x) else "NA"
  delta <- noncdr_mean - cdr_mean
  label <- sprintf(
    "Non-CDR: %s; CDR: %s; delta: %s; CDR length: %.0f kb",
    fmt(noncdr_mean), fmt(cdr_mean), fmt(delta), cdr_len_kb
  )

  # Keep the familiar 0-100 percent scale when applicable, while avoiding
  # clipped data if another valid methylation scale is supplied.
  observed_range <- range(sub$meth, na.rm = TRUE)
  plot_min <- min(0, observed_range[1])
  plot_max <- max(100, observed_range[2])
  marker_height <- 0.03 * (plot_max - plot_min)

  p <- ggplot(sub, aes(x = start, y = meth, group = interaction(active_start, qc_segment))) +
    geom_line(linewidth = 0.25, alpha = 0.5) +
    geom_line(aes(y = meth_smooth), linewidth = 0.45, na.rm = TRUE) +
    geom_rect(
      data = cdr_ch,
      aes(xmin = start, xmax = end, ymin = plot_max, ymax = plot_max + marker_height),
      inherit.aes = FALSE,
      fill = "grey30"
    ) +
    annotate(
      "text", x = min(sub$start), y = plot_min + 0.05 * (plot_max - plot_min),
      label = label, hjust = 0, size = 3
    ) +
    scale_x_continuous(labels = function(x) paste0(format(x / 1e6, trim = TRUE), " Mb")) +
    scale_y_continuous(limits = c(plot_min, plot_max + marker_height)) +
    labs(
      title = paste0(opt$sample, " | ", ch, " active HOR"),
      x = NULL, y = "Methylation (%)"
    ) +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank())

  print(p)
}

grDevices::dev.off()
on.exit(NULL, add = FALSE)
message(sprintf("Wrote %d CDR candidate(s) for %s", nrow(cdr), opt$sample))

