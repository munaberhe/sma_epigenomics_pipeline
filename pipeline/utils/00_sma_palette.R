# SMA colour palette
# Single source of truth for all active scripts.
# Conditions
#   ASO_CTRL      = #1F3A5F  (dark blue)
#   Scramble_CTRL = #6B7280  (grey)
#   ASO_VPA       = #C0392B  (red)       one red only
#   Scramble_VPA  = #F0A500  (amber)
# Direction (hypo/hyper) mirrors condition: blue = hypo, red = hyper
# All reds in active scripts must be #C0392B - not #C0392B, #A84B2F, or #D55E00

SMA_PALETTE <- list(
  conditions = c(
    ASO_CTRL      = "#1F3A5F",
    Scramble_CTRL = "#6B7280",
    ASO_VPA       = "#C0392B",
    Scramble_VPA  = "#F0A500"
  ),
  direction = c(hypo="#1F3A5F", hyper="#C0392B"),
  features = c(
    "Promoter (<=1kb)"="#1F3A5F","Promoter (1-2kb)"="#3A5A85",
    "Promoter (2-3kb)"="#5A7CA8","5' UTR"="#C0392B","Exon"="#F0A500",
    "Intron"="#6B7280","3' UTR"="#A0BFD9","Downstream"="#8B7355",
    "Distal Intergenic"="#DCD9D1","Other"="#CCCCCC"
  ),
  cpg_context = c(
    "Open_Sea"="#1F3A5F","Shore"="#5A7CA8",
    "Shelf"="#A0BFD9","Island"="#C0392B"
  ),
  bg="#F5F4EF", observed="#1F3A5F",
  null_fill="#DCD9D1", null_outline="#6B7280",
  highlight="#1F3A5F", other="#CCCCCC"
)
sma_heatmap_scale <- function(limits=c(-1.5,1.5)) {
  ggplot2::scale_fill_gradient2(
    low="#1F3A5F", mid="#F5F4EF", high="#C0392B",
    midpoint=0, limits=limits, oob=scales::squish,
    name=expression(log[2](obs/exp))
  )
}
