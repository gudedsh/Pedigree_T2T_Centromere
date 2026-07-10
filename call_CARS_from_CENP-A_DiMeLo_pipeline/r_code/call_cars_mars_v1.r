#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(zoo)
  library(scales)
  library(patchwork) 
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx == length(args)) stop("Missing value for ", flag)
  args[idx + 1]
}

in_file  <- get_arg("--input")
cen_file <- get_arg("--cen")
out_file <- get_arg("--out", "CARS_MARS_results.tsv")
plotdir  <- get_arg("--plotdir", "CARS_MARS_QC_plots")

smooth_k <- as.integer(get_arg("--smooth_k", "20"))
bin_size <- as.integer(get_arg("--bin_size", "1000"))

cars_col <- get_arg("--cars_col", "m6A_per_AT")
mars_col <- get_arg("--mars_col", "mCG_per_CpG")

out_CDR_cutoff <- as.numeric(get_arg("--out_CDR_cutoff", "0.10"))
in_CDR_cutoff  <- as.numeric(get_arg("--in_CDR_cutoff", "0.05"))
cdr_expansion  <- as.integer(get_arg("--cdr_expansion", "0"))

if (is.null(in_file) || is.null(cen_file)) {
  stop("Usage: Rscript call_cars_mars.r --input bin_summary.tsv --cen activeHOR_CDR.tsv")
}

dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)

dt <- fread(in_file)
required_cols <- c("chr", "start", "end", "sample", cars_col, mars_col)

if (!all(required_cols %in% colnames(dt))) {
  stop("Input file must contain columns: ", paste(required_cols, collapse = ", "))
}

dt[, cars_signal := as.numeric(get(cars_col))]
dt[, mars_signal := as.numeric(get(mars_col))]
dt <- dt[is.finite(cars_signal) & is.finite(mars_signal)]


setorder(dt, chr, sample, start)

if (smooth_k > 1) {
  dt[, cars_smooth := zoo::rollmean(cars_signal, k = smooth_k, fill = NA, align = "center"),
     by = .(chr, sample)]
  dt[, mars_smooth := zoo::rollmean(mars_signal, k = smooth_k, fill = NA, align = "center"),
     by = .(chr, sample)]
} else {
  dt[, cars_smooth := cars_signal]
  dt[, mars_smooth := mars_signal]
}

dt <- dt[is.finite(cars_smooth) & is.finite(mars_smooth)]

cen <- fread(cen_file, header = FALSE)

if (ncol(cen) < 5) {
  stop("CEN file must have at least 5 columns: chr active_start active_end cdr_start cdr_end")
}

setnames(cen, 1:5, c("chr", "active_start", "active_end", "cdr_start_raw", "cdr_end_raw"))

cen <- cen[, .(
  chr = as.character(chr),
  active_start = as.integer(active_start),
  active_end = as.integer(active_end),
  cdr_start_raw = as.integer(cdr_start_raw),
  cdr_end_raw = as.integer(cdr_end_raw)
)]

cen[, cdr_start := pmax(active_start, cdr_start_raw - cdr_expansion)]
cen[, cdr_end := pmin(active_end, cdr_end_raw + cdr_expansion)]
cen[, cdr_len := cdr_end - cdr_start]

if (any(cen$cdr_end <= cen$cdr_start)) {
  stop("Invalid CDR after expansion.")
}

winsorize_vec <- function(x, trim_fraction) {
  x2 <- x
  if (is.null(trim_fraction) || trim_fraction <= 0) return(x2)
  finite_x <- x[is.finite(x)]
  if (length(finite_x) == 0) return(x2)
  
  low <- quantile(finite_x, trim_fraction, na.rm = TRUE)
  high <- quantile(finite_x, 1 - trim_fraction, na.rm = TRUE)
  
  x2[is.finite(x2) & x2 < low] <- low
  x2[is.finite(x2) & x2 > high] <- high
  x2
}

calc_one_score <- function(
    x,
    signal_col,
    cdr_start,
    cdr_end,
    active_start,
    active_end,
    out_cutoff,
    in_cutoff,
    direction = c("cdr_high", "cdr_low")
) {
  direction <- match.arg(direction)
  
  x <- x[start >= active_start & end <= active_end]
  
  empty_result <- function() {
    data.table(
      bg_signal = NA_real_,
      cdr_signal = NA_real_,
      bg_signal_raw = NA_real_,
      cdr_signal_raw = NA_real_,
      signal_difference = NA_real_,
      area_score = NA_real_,
      log2_area_score = NA_real_,
      mean_delta = NA_real_,
      n_active_bins = nrow(x),
      n_bg_bins_raw = 0,
      n_cdr_bins_raw = 0
    )
  }
  
  if (nrow(x) == 0) return(empty_result())
  
  bg <- x[end <= cdr_start | start >= cdr_end]
  cdr <- x[start >= cdr_start & end <= cdr_end]
  
  if (nrow(bg) == 0 || nrow(cdr) == 0) return(empty_result())
  
  bg_raw <- mean(bg[[signal_col]], na.rm = TRUE)
  cdr_raw <- mean(cdr[[signal_col]], na.rm = TRUE)
  
  bg_used <- copy(bg)
  cdr_used <- copy(cdr)
  
  bg_used[, signal_winsor := winsorize_vec(get(signal_col), out_cutoff)]
  cdr_used[, signal_winsor := winsorize_vec(get(signal_col), in_cutoff)]
  
  bg_mean <- mean(bg_used$signal_winsor, na.rm = TRUE)
  cdr_mean <- mean(cdr_used$signal_winsor, na.rm = TRUE)
  
  if (direction == "cdr_high") {
    cdr_used[, delta := signal_winsor - bg_mean]
    diff <- cdr_mean - bg_mean
  } else {
    cdr_used[, delta := bg_mean - signal_winsor]
    diff <- bg_mean - cdr_mean
  }
  
  score <- sum(cdr_used$delta, na.rm = TRUE)
  
  data.table(
    bg_signal = round(bg_mean, 6),
    cdr_signal = round(cdr_mean, 6),
    bg_signal_raw = round(bg_raw, 6),
    cdr_signal_raw = round(cdr_raw, 6),
    signal_difference = round(diff, 6),
    area_score = round(score, 6),
    log2_area_score = round(log2(pmax(score, 0) + 1e-9), 3),
    mean_delta = round(mean(cdr_used$delta, na.rm = TRUE), 6),
    n_active_bins = nrow(x),
    n_bg_bins_raw = nrow(bg),
    n_cdr_bins_raw = nrow(cdr)
  )
}

get_cell_group <- function(x) {
  fifelse(grepl("PBMC", x, ignore.case = TRUE), "PBMC",
  fifelse(grepl("iPSC", x, ignore.case = TRUE), "iPSC",
  fifelse(grepl("NPC", x, ignore.case = TRUE), "NPC",
  fifelse(grepl("Monocyte", x, ignore.case = TRUE), "Monocyte",
  fifelse(grepl("Macrophage", x, ignore.case = TRUE), "Macrophage", "Other")))))
}

sample_order_by_celltype <- function(samples) {
  tmp <- data.table(sample = unique(samples))
  tmp[, cell := get_cell_group(sample)]
  tmp[, rank := match(cell, c("PBMC", "iPSC", "NPC", "Monocyte", "Macrophage", "Other"))]
  tmp[order(rank, sample)]$sample
}



plot_one_chr_facet_cars_mars <- function(plot_dt, cen_row, res_chr, outfile) {
  if (nrow(plot_dt) == 0) return(NULL)

  plot_dt <- copy(plot_dt)
  plot_dt[, cell := get_cell_group(sample)]

  sample_order <- sample_order_by_celltype(plot_dt$sample)
  plot_dt[, sample := factor(sample, levels = sample_order)]

  res_chr <- copy(res_chr)
  res_chr[, sample := factor(sample, levels = sample_order)]

  label_dt <- res_chr[, .(
    sample,
    label_cars = paste0("log2(CARS) = ", round(log2_CARS, 3)),
    label_mars = paste0("log2(MARS) = ", round(log2_MARS, 3))
  )]

  bg_dt <- res_chr[, .(
    sample,
    bg_cars = bg_CARS_signal,
    bg_mars = bg_MARS_signal
  )]

  ## CARS area: CENP-A signal above background inside CDR
  cars_delta_list <- list()
  mars_delta_list <- list()

  for (sid in sample_order) {
    x <- plot_dt[sample == sid]

    cdr <- x[start >= cen_row$cdr_start & end <= cen_row$cdr_end]
    if (nrow(cdr) == 0) next

    bg_cars <- res_chr[sample == sid, bg_CARS_signal][1]
    bg_mars <- res_chr[sample == sid, bg_MARS_signal][1]

    cdr_cars <- copy(cdr)
    cdr_cars[, background := bg_cars]
    cdr_cars[, delta := cars_smooth - background]
    cdr_cars <- cdr_cars[is.finite(delta) & delta > 0]
    cars_delta_list[[as.character(sid)]] <- cdr_cars

    cdr_mars <- copy(cdr)
    cdr_mars[, background := bg_mars]
    cdr_mars[, delta := background - mars_smooth]
    cdr_mars <- cdr_mars[is.finite(delta) & delta > 0]
    mars_delta_list[[as.character(sid)]] <- cdr_mars
  }

  cars_delta_dt <- rbindlist(cars_delta_list, fill = TRUE)
  mars_delta_dt <- rbindlist(mars_delta_list, fill = TRUE)

  if (nrow(cars_delta_dt) > 0) {
    cars_delta_dt[, sample := factor(sample, levels = sample_order)]
  }
  if (nrow(mars_delta_dt) > 0) {
    mars_delta_dt[, sample := factor(sample, levels = sample_order)]
  }

  cars_rng <- quantile(plot_dt$cars_smooth, c(0.01, 0.99), na.rm = TRUE)
  mars_rng <- quantile(plot_dt$mars_smooth, c(0.01, 0.99), na.rm = TRUE)

  p_cars <- ggplot() +
    geom_rect(
      data = cen_row,
      aes(xmin = cdr_start, xmax = cdr_end, ymin = -Inf, ymax = Inf),
      fill = "grey85",
      alpha = 0.5,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = plot_dt,
      aes(x = start, y = cars_smooth),
      size = 0.1,
      color = "grey50",
      alpha = 0.5
    ) +
    geom_hline(
      data = bg_dt,
      aes(yintercept = bg_cars),
      linetype = "dashed",
      color = "#D55E00",
      linewidth = 0.35
    ) +
    geom_segment(
      data = cars_delta_dt,
      aes(
        x = start,
        xend = start,
        y = background,
        yend = cars_smooth
      ),
      color = "#D55E00",
      linewidth = 0.1,
      alpha = 0.5
    ) +
    geom_text(
      data = label_dt,
      aes(x = Inf, y = cars_rng[2], label = label_cars),
      hjust = 1.05,
      vjust = 1.1,
      size = 3,
      inherit.aes = FALSE
    ) +
    facet_wrap(~sample, ncol = 1) +
    scale_x_continuous(labels = label_number(scale = 1e-6)) +
 #   coord_cartesian(ylim = cars_rng) +
    labs(
      x = NULL,
      y = cars_col,
      title = paste0(cen_row$chr, ": CARS / CENP-A enrichment")
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      strip.text = element_text(size = 8, face = "bold"),
      plot.title = element_text(face = "bold")
    )

  p_mars <- ggplot() +
    geom_rect(
      data = cen_row,
      aes(xmin = cdr_start, xmax = cdr_end, ymin = -Inf, ymax = Inf),
      fill = "grey85",
      alpha = 0.5,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = plot_dt,
      aes(x = start, y = mars_smooth),
      size = 0.1,
      color = "grey50",
      alpha = 0.5
    ) +
    geom_hline(
      data = bg_dt,
      aes(yintercept = bg_mars),
      linetype = "dashed",
      color = "#0072B2",
      linewidth = 0.35
    ) +
    geom_segment(
      data = mars_delta_dt,
      aes(
        x = start,
        xend = start,
        y = mars_smooth,
        yend = background
      ),
      color = "#0072B2",
      linewidth = 0.1,
      alpha = 0.5
    ) +
    geom_text(
      data = label_dt,
      aes(x = Inf, y = 0.5, label = label_mars),
      hjust = 1.05,
      vjust = 1.1,
      size = 3,
      inherit.aes = FALSE
    ) +
    facet_wrap(~sample, ncol = 1) +
    scale_x_continuous(labels = label_number(scale = 1e-6)) +
  #  coord_cartesian(ylim = mars_rng) +
    labs(
      x = "Genomic position (Mb)",
      y = mars_col,
      title = paste0(cen_row$chr, ": MARS / mCG depletion")
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      strip.text = element_text(size = 8, face = "bold"),
      plot.title = element_text(face = "bold")
    )

  p <- p_cars / p_mars +
    plot_layout(heights = c(1, 1))

  print(p)

  ggsave(
    outfile,
    p,
    width = 8,
    height = max(6, 2.4 * uniqueN(plot_dt$sample))
  )
}


res_list <- list()

for (i in seq_len(nrow(cen))) {
  cr <- cen[i]
  
  chr_dt <- dt[
    chr == cr$chr &
      start >= cr$active_start &
      end <= cr$active_end
  ]
  
  if (nrow(chr_dt) == 0) next
  
  nsamp <- uniqueN(chr_dt$sample)
  
  chr_dt_common <- chr_dt[
    ,
    if (uniqueN(sample) == nsamp) .SD,
    by = .(chr, start, end)
  ]
  
  for (sid in unique(chr_dt_common$sample)) {
    x <- chr_dt_common[sample == sid]
    
    cars <- calc_one_score(
      x = x,
      signal_col = "cars_smooth",
      cdr_start = cr$cdr_start,
      cdr_end = cr$cdr_end,
      active_start = cr$active_start,
      active_end = cr$active_end,
      out_cutoff = out_CDR_cutoff,
      in_cutoff = in_CDR_cutoff,
      direction = "cdr_high"
    )
    
    mars <- calc_one_score(
      x = x,
      signal_col = "mars_smooth",
      cdr_start = cr$cdr_start,
      cdr_end = cr$cdr_end,
      active_start = cr$active_start,
      active_end = cr$active_end,
      out_cutoff = out_CDR_cutoff,
      in_cutoff = in_CDR_cutoff,
      direction = "cdr_low"
    )
    
    setnames(cars, names(cars), paste0(c(
      "bg", "cdr", "bg_raw", "cdr_raw", "diff", "score", "log2", "mean_delta",
      "n_active", "n_bg", "n_cdr"
    ), "_CARS"))
    
    setnames(mars, names(mars), paste0(c(
      "bg", "cdr", "bg_raw", "cdr_raw", "diff", "score", "log2", "mean_delta",
      "n_active", "n_bg", "n_cdr"
    ), "_MARS"))
    
    one <- cbind(
      data.table(
        sample = sid,
        cell = get_cell_group(sid),
        chr = cr$chr,
        active_start = cr$active_start,
        active_end = cr$active_end,
        cdr_start = cr$cdr_start,
        cdr_end = cr$cdr_end,
        cdr_start_raw = cr$cdr_start_raw,
        cdr_end_raw = cr$cdr_end_raw,
        cdr_len = cr$cdr_len,
        bin_size = bin_size,
        smooth_k = smooth_k,
        cars_col = cars_col,
        mars_col = mars_col,
        out_CDR_cutoff = out_CDR_cutoff,
        in_CDR_cutoff = in_CDR_cutoff,
        cdr_expansion = cdr_expansion,
        common_sample_n = nsamp
      ),
      cars,
      mars
    )
    
    res_list[[length(res_list) + 1]] <- one
  }
  
  res_tmp <- rbindlist(res_list, fill = TRUE)
  res_chr <- res_tmp[chr == cr$chr]
  
  if (nrow(res_chr) > 0) {
    setnames(
      res_chr,
      old = c("log2_CARS", "score_CARS", "bg_CARS", "cdr_CARS",
              "log2_MARS", "score_MARS", "bg_MARS", "cdr_MARS"),
      new = c("log2_CARS", "CARS", "bg_CARS_signal", "cdr_CARS_signal",
              "log2_MARS", "MARS", "bg_MARS_signal", "cdr_MARS_signal"),
      skip_absent = TRUE
    )
    
    plot_file <- file.path(plotdir, paste0(cr$chr, "_CARS_MARS_QC_facet.pdf"))
    plot_one_chr_facet_cars_mars(chr_dt_common, cr, res_chr, plot_file)
  }
}

res <- rbindlist(res_list, fill = TRUE)

setnames(
  res,
  old = c(
    "score_CARS", "log2_CARS", "diff_CARS", "bg_CARS", "cdr_CARS",
    "score_MARS", "log2_MARS", "diff_MARS", "bg_MARS", "cdr_MARS"
  ),
  new = c(
    "CARS", "log2_CARS", "CARS_difference", "bg_CARS_signal", "cdr_CARS_signal",
    "MARS", "log2_MARS", "MARS_difference", "bg_MARS_signal", "cdr_MARS_signal"
  ),
  skip_absent = TRUE
)

setcolorder(
  res,
  c(
    "sample", "cell", "chr",
    "active_start", "active_end",
    "cdr_start", "cdr_end", "cdr_len",
    "CARS", "log2_CARS", "CARS_difference",
    "bg_CARS_signal", "cdr_CARS_signal",
    "MARS", "log2_MARS", "MARS_difference",
    "bg_MARS_signal", "cdr_MARS_signal",
    "bin_size", "smooth_k", "cars_col", "mars_col",
    "out_CDR_cutoff", "in_CDR_cutoff", "cdr_expansion",
    "common_sample_n"
  ),
  skip_absent = TRUE
)

fwrite(res, out_file, sep = "\t")

res

library(data.table)
library(ggplot2)

##------------------------------------------
## Calculate correlation for each sample
##------------------------------------------
res[, cell := fifelse(

    grepl("iPSC", sample), "iPSC",

    fifelse(grepl("NPC", sample), "NPC", "Other")

)]

cor_dt <- res[, {
  ct <- cor.test(MARS, CARS, method = "pearson")

  .(
    R = unname(ct$estimate),
    P = ct$p.value
  )
}, by = cell]


cor_dt[, label :=

    sprintf("%s: R = %.2f, p = %.3g",

            cell, R, P)]

## label position (staggered)
cor_dt[, x := min(res$MARS) + 2]
cor_dt[, y := max(res$CARS) - c(0.15,0.55)[1:.N]]

cor_dt[, label := sprintf("R = %.2f, p = %.3g", R, P)]

##------------------------------------------
## Plot
##------------------------------------------

p_summary <- ggplot(res,aes(x = MARS,y = CARS,color = cell,group=cell)) +
  geom_point(size = 2.8,alpha = 0.9) +
  ## regression for each sample
  geom_smooth(method = "lm",se = FALSE,linewidth = 0.9) +
  ## correlation labels
  geom_text(data = cor_dt,aes(x = x,y = y,label = label,color = cell),hjust = 0,size = 4,show.legend = FALSE) +
  labs(
    x = "MARS: mCG depletion area score",
    y = "CARS: CENP-A enrichment area score",
    title = "Relationship between CARS and MARS",
    color = "cell") +
  theme_classic(base_size = 13) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(
  file.path(plotdir, "CARS_vs_MARS_summary.pdf"),
  p_summary,
  width = 7,
  height = 4.5
)

message("Done.")
message("Result written to: ", out_file)
message("QC plots written to: ", plotdir)







