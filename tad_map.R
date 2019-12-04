library(rtracklayer)
library(EnrichedHeatmap)
library(purrr)
library(cowplot)

mapa <- import("hg19.wgEncodeCrgMapabilityAlign100mer.bigWig", format = "BigWig")

tads <- list(
    gm12 = import("hg19.TADs/GM12878_Lieberman-raw_TADs.bed", format = "BED", genome = "hg19"),
    h1es = import("hg19.TADs/H1-ESC_Dixon2015-raw_TADs.txt", format = "BED", genome = "hg19"),
    live = import("hg19.TADs/Liver_STL011_Leung2015-raw_TADs.txt", format = "BED", genome = "hg19")
)

mats <- purrr::map(
    tads,
    ~normalizeToMatrix(mapa, .x, extend = 1000000, w = 25000, value_column = "score", mean_mode = "w0")
)

titles <- c("GM12878_Lieberman", "H1-ESC_Dixon2015", "Liver_Leung2015")

png("p1.png", width = 200, height = 400)
EnrichedHeatmap(
    mats[[1]],
    row_order = order(enriched_score(mats[[1]]), decreasing = FALSE),
    col = rev(c("white", "orange", "black")),
    column_title = titles[[1]],
    name = "100mer mappability",
    show_heatmap_legend = FALSE
)
dev.off()
png("p2.png", width = 200, height = 400)
EnrichedHeatmap(
    mats[[2]],
    row_order = order(enriched_score(mats[[2]]), decreasing = FALSE),
    col = rev(c("white", "orange", "black")),
    column_title = titles[[2]],
    name = "100mer mappability",
    show_heatmap_legend = FALSE
)
dev.off()
png("p3.png", width = 300, height = 400)
EnrichedHeatmap(
    mats[[3]],
    row_order = order(enriched_score(mats[[3]]), decreasing = FALSE),
    col = rev(c("white", "orange", "black")),
    column_title = titles[[3]],
    name = "100mer mappability"
)
dev.off()

summary(width(tads))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 100000   375000   625000   842586  1050000 30550000





mat1 <- normalizeToMatrix(
    mapa,
    narrow(tads, start = 1, end = 2),
    extend = 1000000, w = 50000, value_column = "score", mean_mode = "w0"
)

EnrichedHeatmap(mat1, name = "Mappability at TAD border")


mat2 <- normalizeToMatrix(
    mapa,
    tads,
    extend = 1000000, w = 50000, value_column = "score", mean_mode = "w0"
)

EnrichedHeatmap(
    mat2,
    row_order = order(enriched_score(mat), decreasing = FALSE),
)
