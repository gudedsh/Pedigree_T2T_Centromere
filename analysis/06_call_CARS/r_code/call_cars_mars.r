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
out_file <- get_arg("--out", "call_CARS_results.tsv")
plotdir  <- get_arg("--plotdir", "plots_QC")

smooth_k <- as.integer(get_arg("--smooth_k", "10"))
bin_size <- as.integer(get_arg("--bin_size", "5000"))

cars_col <- get_arg("--cars_col", "m6A_per_AT")

out_CDR_cutoff <- as.numeric(get_arg("--out_CDR_cutoff", "0.10"))
in_CDR_cutoff  <- as.numeric(get_arg("--in_CDR_cutoff", "0.05"))
cdr_expansion  <- as.integer(get_arg("--cdr_expansion", "0"))

if (is.null(in_file) || is.null(cen_file)) {
  stop("Usage: Rscript call_cars_mars.r --input bin_summary.tsv --cen activeHOR_CDR.tsv")
}

dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)

out_summary <- paste0(plotdir,"/results_summary.tsv")

winsorize_vec <- function(x,p){
  if(is.null(p)||p<=0)return(x) 
  fx <- x[is.finite(x)]
  if(!length(fx))return(x)
  q <- quantile(fx,c(p,1-p),na.rm=TRUE)
  pmin(pmax(x,q[1]),q[2])
}

dt <- fread(in_file)
req <- c("chr","start","end","sample",cars_col)

if(!all(req %in% names(dt))) stop("Input file must contain: ",paste(req,collapse=", "))

dt[, cars_signal := as.numeric(get(cars_col))]

dt <- dt[is.finite(cars_signal)]

setorder(dt,chr,sample,start)

if(smooth_k>1){
  dt[, cars_smooth:=rollmean(cars_signal,k=smooth_k,fill=NA,align="center"),by=.(chr,sample)]
} else dt[, cars_smooth := cars_signal]

dt <- dt[is.finite(cars_smooth)]

cen <- fread(cen_file,header=FALSE)
if(ncol(cen)<5) stop("CEN file must have at least 5 columns: chr active_start active_end cdr_start cdr_end")
setnames(cen,1:5,c("chr","active_start","active_end","cdr_start_raw","cdr_end_raw"))
cen <- cen[,.(chr=as.character(chr),active_start=as.integer(active_start),active_end=as.integer(active_end),cdr_start_raw=as.integer(cdr_start_raw),cdr_end_raw=as.integer(cdr_end_raw))]
cen[, `:=`(cdr_start=pmax(active_start,cdr_start_raw-cdr_expansion),cdr_end=pmin(active_end,cdr_end_raw+cdr_expansion))]
cen[, cdr_len:=cdr_end-cdr_start]

if(any(cen$cdr_end<=cen$cdr_start)) stop("Invalid CDR after expansion.")

calc_one_score <- function(x,signal_col,cdr_start,cdr_end,active_start,active_end,out_cutoff,in_cutoff){
  x <- x[start>=active_start & end<=active_end]
  bg <- copy(x[end<=cdr_start | start>=cdr_end])
  cdr <- copy(x[start>=cdr_start & end<=cdr_end])
  
  empty <- data.table(
    CARS=NA_real_,
    log2_CARS=NA_real_,
    Fold_Enrichment=NA_real_,
    bg_signal=NA_real_,
    cdr_signal=NA_real_,
    bg_signal_raw=NA_real_,
    cdr_signal_raw=NA_real_,
    n_active_bins=nrow(x),
    n_bg_bins_raw=nrow(bg),
    n_cdr_bins_raw=nrow(cdr)
    )
  
  if(!nrow(x)||!nrow(bg)||!nrow(cdr)) return(empty)
  
  bg_raw <- mean(bg[[signal_col]],na.rm=TRUE)
  cdr_raw <- mean(cdr[[signal_col]],na.rm=TRUE)
  
  bg[, signal_w:=winsorize_vec(get(signal_col),out_cutoff)]
  cdr[, `:=`(signal_w=winsorize_vec(get(signal_col),in_cutoff),bin_width=end-start)]
  
  bg_mean <- mean(bg$signal_w,na.rm=TRUE)
  cdr_mean <- mean(cdr$signal_w,na.rm=TRUE)
  
  if(!is.finite(bg_mean)||!is.finite(cdr_mean)) return(empty)
  
  eps <- 1e-6
  cdr[, delta := (signal_w - bg_mean) / (bg_mean + eps) ]
  diff <- cdr_mean - bg_mean
  fold_enrichment <- cdr_mean/(bg_mean + eps)
  score <- sum(cdr$delta * cdr$bin_width,na.rm=TRUE)

  data.table(
    CARS=round(score,1),
    log2_CARS=round(sign(score)*log2(abs(score)+1),3),
    Fold_Enrichment = round(fold_enrichment, 2),
    bg_signal=round(bg_mean,4),
    cdr_signal=round(cdr_mean,4),
    bg_signal_raw=round(bg_raw,4),
    cdr_signal_raw=round(cdr_raw,4),
    n_active_bins=nrow(x),
    n_bg_bins_raw=nrow(bg),
    n_cdr_bins_raw=nrow(cdr))
}


plot_one_chr <- function(plot_dt,cr,res_chr){
  if(!nrow(plot_dt)) return(NULL)

  res_chr[,label := paste0("log2(CARS) = ",log2_CARS)]

  p1 <- ggplot()+
    geom_rect(data=cr,aes(xmin=cdr_start,xmax=cdr_end,ymin=-Inf,ymax=Inf),inherit.aes=FALSE,color="lightgrey",alpha=0.1)+
    geom_line(data=plot_dt,aes(x=start,y=cars_smooth_nor_bg,color=sample),linewidth=0.5,alpha=0.8)+
    geom_hline(aes(yintercept=0),linetype="dashed",linewidth=.35)+
    geom_text(data=res_chr,aes(x=Inf,y=Inf,label=label,color=sample),inherit.aes=FALSE,hjust=1.2,vjust=seq(1.5,3.5,length.out=nrow(res_chr)),size=2.8,show.legend=FALSE)+
    scale_x_continuous(labels=label_number(scale=1e-6))+
    labs(x=NULL,y="log2(CARS)",title=paste0(cr$asm," " ,cr$chr,": CARS / CENP-A enrichment"))+
    theme_bw(base_size=9)+
    theme(
      panel.grid=element_blank(),
      plot.title=element_text(face="bold"),
      legend.position="top")

    p1

}

res_dt <- data.table()
plot_list <- list()

for(i in seq_len(nrow(cen))){
  cr <- cen[i]
  chr_dt <- dt[chr==cr$chr & start>=cr$active_start & end<=cr$active_end]
  if(!nrow(chr_dt)) next
  nsamp <- uniqueN(chr_dt$sample)
  chr_dt <- chr_dt[, if(uniqueN(sample)==nsamp) .SD, by=.(chr,start,end)]
  eps <- 1e-6

  sid <- unique(chr_dt$sample)[1]
  res_chr <- data.frame()
  
  for(sid in unique(chr_dt$sample)){
    x <- chr_dt[sample==sid]
    cars <- calc_one_score(x,"cars_smooth",cr$cdr_start,cr$cdr_end,cr$active_start,cr$active_end,out_CDR_cutoff,in_CDR_cutoff)
    cars[,sample := sid]
    cars[,chr := cr$chr]
    res_chr <- rbind(res_chr,cars)
    chr_dt[sample==sid, cars_smooth_nor_bg := (cars_smooth - cars$bg_signal) / (cars$bg_signal + eps)]
  }

    res_dt <- rbind(res_dt, res_chr)

    p <- plot_one_chr(chr_dt,cr,res_chr)
    plot_list[[cr$chr]] <- p
}

## output results
setcolorder(res_dt,c("sample","chr","CARS","log2_CARS","Fold_Enrichment","bg_signal","cdr_signal","bg_signal_raw","cdr_signal_raw","n_active_bins","n_bg_bins_raw","n_cdr_bins_raw" ),skip_absent=TRUE) 

res_dt[,':='(bin_size=bin_size,
             smooth_k=smooth_k,
             out_CDR_cutoff=out_CDR_cutoff,
             in_CDR_cutoff=in_CDR_cutoff,
             cdr_expansion=cdr_expansion)]

fwrite(res_dt,out_file,sep="\t",col.names = T,row.names = F,quote = F)

summary_dt <- rbindlist(
  lapply(c("log2_CARS", "Fold_Enrichment"), function(x) {
    res_dt[, .(
      measure = x,
      n = .N,
      mean   = round(mean(get(x), na.rm = TRUE), 3),
      sd     = round(sd(get(x), na.rm = TRUE), 3),
      median = round(median(get(x), na.rm = TRUE), 3),
      min    = round(min(get(x), na.rm = TRUE), 3),
      max    = round(max(get(x), na.rm = TRUE), 3)
    ), by = sample]
  })
)
 
fwrite(summary_dt,out_summary,sep="\t",col.names = T,row.names = F,quote = F)

## plots
res_long <- melt(res_dt,id.vars = c("sample","chr"),measure.vars = c("log2_CARS","Fold_Enrichment"),variable.name = "sam_chr",value.name = "score")

library(data.table)
library(ggplot2)

res_long <- melt(
  res_dt,
  id.vars = c("sample", "chr"),
  measure.vars = c("log2_CARS", "Fold_Enrichment"),
  variable.name = "measure",
  value.name = "score"
)

res_long[, measure := factor(measure,levels = c("log2_CARS", "Fold_Enrichment"),  labels = c("log2(CARS)", "Fold enrichment"))]

p_summary <- ggplot(res_long, aes(x = sample, y = score, color = sample)) +
  geom_boxplot(width = 0.55, outlier.shape = NA,linewidth = 0.4) +
  geom_jitter(width = 0.18,size = 0.8,  alpha = 0.35) +
  facet_wrap(~measure, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = "Score") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

p_summary
p_by_chr <- ggplot(res_long, aes(x = chr, y = score, color = sample)) +
  geom_boxplot(aes(group = interaction(chr, sample)),width = 0.65,
    outlier.shape = NA,
    linewidth = 0.35,
    position = position_dodge(width = 0.75)
  ) +
  facet_wrap(~measure, scales = "free_y", ncol = 1) +
  labs(x = NULL, y = "Score") +
  theme_classic(base_size = 12) +
  theme(
      legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )

## output qc plot pdf
{
pdf(file.path(plotdir,"call_CARS_QC_plots.pdf"),width=8,height=6)
for(ch in names(plot_list)) print(plot_list[[ch]])
dev.off()
}

{
pdf(file.path(plotdir,"summary_plot.pdf"),width=8,height=6)
print(p_summary)
dev.off()
}

{
pdf(file.path(plotdir,"per_chr_plot.pdf"),width=12,height=6)
print(p_by_chr)
dev.off()
}

message("Done!")


