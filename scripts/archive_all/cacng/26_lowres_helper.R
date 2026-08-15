plot_locus_baseR <- function(pooled, CONDITIONS, region, winsize=300,
                              genes_df=NULL, dmr_start=NULL, dmr_end=NULL,
                              title="", out_path, width=11, height=6) {
  suppressPackageStartupMessages(library(DMRcaller))

  COND_COLOURS <- c(ASO_CTRL="#1F3A5F", Scramble_CTRL="#6B7280",
                    ASO_VPA="#C0392B",   Scramble_VPA="#D4A017")
  COND_LTY <- c(ASO_CTRL=1, Scramble_CTRL=2, ASO_VPA=1, Scramble_VPA=2)

  prof_list <- lapply(CONDITIONS, function(cond) {
    if (is.null(pooled[[cond]])) return(NULL)
    prof <- computeMethylationProfile(pooled[[cond]], region, winsize, "CG")
    df <- as.data.frame(prof)
    df$meth <- df$sumReadsM / df$sumReadsN
    df$pos  <- (df$start + df$end) / 2
    df[!is.na(df$meth) & df$sumReadsN >= 3, ]
  })
  names(prof_list) <- CONDITIONS
  prof_list <- Filter(Negate(is.null), prof_list)
  if (length(prof_list)==0) { message("no data"); return(NULL) }

  x_range <- range(unlist(lapply(prof_list, `[[`, "pos")), na.rm=TRUE)
  chr_label <- as.character(seqnames(region))

  cairo_pdf(out_path, width=width, height=height)

  if (!is.null(genes_df) && nrow(genes_df)>0) {
    layout(matrix(1:2, nrow=2), heights=c(0.12, 0.88))
    par(mar=c(0,5,2,1))
    plot(NULL, xlim=x_range, ylim=c(0,1),
         xaxt="n", yaxt="n", xlab="", ylab="", bty="n",
         main=title, cex.main=1.0, font.main=2)
    for (i in seq_len(nrow(genes_df))) {
      rect(genes_df$start[i], 0.1, genes_df$end[i], 0.9,
           col="#1F3A5F", border=NA)
      text((genes_df$start[i]+genes_df$end[i])/2, 0.5,
           genes_df$name[i], col="white", cex=0.7, font=2)
    }
    par(mar=c(4,5,0.5,1))
  } else {
    par(mar=c(4,5,3,1))
  }

  plot(NULL, xlim=x_range, ylim=c(0,1),
       xlab=paste0(chr_label," position (bp)"),
       ylab="CpG methylation",
       main=if(is.null(genes_df)||nrow(genes_df)==0) title else "")

  # DMR shading
  if (!is.null(dmr_start) && !is.null(dmr_end)) {
    rect(dmr_start, -0.05, dmr_end, 1.05,
         col=adjustcolor("#C0392B", 0.1), border=NA)
    abline(v=c(dmr_start, dmr_end), col="#C0392B", lty=3, lwd=0.8)
  }

  for (cond in names(prof_list)) {
    d <- prof_list[[cond]]
    lines(d$pos, d$meth, col=COND_COLOURS[cond],
          lwd=1.8, lty=COND_LTY[cond])
  }

  legend("topright", legend=names(prof_list),
         col=COND_COLOURS[names(prof_list)],
         lty=COND_LTY[names(prof_list)],
         lwd=1.8, bty="n", cex=0.85)

  dev.off()
  message("saved: ", basename(out_path))
}
