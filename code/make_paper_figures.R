
# Load packages
library(here)
library(ggplot2)
library(patchwork)
library(rphylopic)
library(magick)

set.seed(123)

# Fig. 1: classification scheme + runtime benchmark ----
## 1A: classification scheme ----
fig_1a <- magick::image_read_pdf(
    here("products", "figs", "classification_scheme.pdf")
) 

## 1B: Runtime benchmark, classification ----
load(here("products", "result_files", "benchmark_classification.rda"))

### Get species silhouettes from phylopic
benchmark_classification$uuid <- sapply(
    benchmark_classification$species, function(x) {
        get_uuid(gsub(" ", "_", x))
    }
)
benchmark_classification$fig <- lapply(benchmark_classification$uuid, get_phylopic)

### Plot data
fig_1b <- ggplot(
    benchmark_classification, 
    aes(x = time_seconds, y = species)
) +
    geom_bar(stat = "identity", fill = "dodgerblue4", color = "black") +
    geom_text(aes(label = round(time_seconds, 2)), hjust = -0.2) +
    geom_phylopic(
        aes(img = fig, x = time_seconds + 1.4), 
        size = 0.5
    ) +
    theme_bw() +
    labs(
        x = "Time (in seconds)", y = "",
        title = "Runtime benchmark for *classify_gene_pairs()*"
    ) +
    theme(
        plot.title = ggtext::element_markdown(),
        panel.grid = element_blank(),
        axis.text.y = element_text(face = "italic", color = "black")
    ) +
    scale_x_continuous(
        limits = c(0, 10), breaks = seq(0, 10, by=2), expand = c(0, 0)
    )

fig_1b

## 1C: Runtime benchmark, substitution rate calculation ----
load(here("products", "result_files", "benchmark_kaks.rda"))

fig_1c <- ggplot(
    benchmark_kaks, aes(x = Time_minutes, y = Back.end)
) +
    geom_bar(stat = "identity", fill = "dodgerblue4", color = "black") +
    geom_text(
        aes(label = paste0(round(Pairs_per_second, 2), " pairs/second")), 
        hjust = -0.1
    ) +
    labs(
        x = "Time (in minutes)", y = "BiocParallel back-end",
        title = "Runtime benchmark for *pairs2kaks()*",
        subtitle = "*S. cerevisiae* paranome; *N* = 3588 gene pairs"
    ) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        plot.title = ggtext::element_markdown(),
        plot.subtitle = ggtext::element_markdown()
    ) +
    scale_x_continuous(
        limits = c(0, 60), breaks = seq(0, 60, by = 10), expand = c(0, 0)
    )

fig_1c

## Combining ABC ----
### Save Fig 1BC to a temporary PDF file and read it with magick
fig_1bc <- wrap_plots(fig_1b, fig_1c, nrow = 1, widths = c(1.5, 1)) +
    plot_annotation(tag_levels = list(c("B", "C"))) &
    theme(plot.tag = element_text(size = 20))
tempfig <- tempfile(fileext = ".pdf")
ggsave(
    plot = fig_1bc, 
    filename = tempfig,
    width = 12, height = 4
)
file2 <- image_read_pdf(tempfig)

### Combine PDFs with magick
fig1_final <- image_append(
    c(
        image_resize(fig_1a, paste0(image_info(file2)$width, "x")), 
        file2
    ), 
    stack = TRUE
) |>
    image_annotate("A", size = 100, location = "+40+40")
    

### Save figure
image_write(
    fig1_final, 
    path = here("products", "figs", "fig1.png"),
    density = 300
)


# Fig. 2: the duplicated gene repertoire across the Eukarya tree of life ----
load(here("products", "plots", "p_duplicates_all_ensembl.rda"))

ggsave(
    p_duplicates_all_ensembl,
    filename = here("products", "figs", "fig2.png"),
    width = 18, height = 14, dpi = 300
)


# Fig. 3: substitution rates by mode + rates in a phylogenetic context ----
load(here("products", "plots", "p_ks_legumes.rda"))
load(here("products", "plots", "p_rates_phylogeny.rda"))

fig3 <- wrap_plots(
    p_rates_phylogeny,
    p_ks_legumes, 
    nrow = 2,
    heights = c(1, 2)
) +
    plot_annotation(tag_levels = list(c("A", "", "", "", "B", "C")))


ggsave(
    fig3, 
    filename = here("products", "figs", "fig3.png"),
    width = 11, height = 10, dpi = 300
)


# Sup. Fig. S1: association between % BUSCOs and % SD-derived genes
load(here("products", "plots", "p_busco_association.rda"))

for(p in seq_along(p_busco_association)[-1]) {
    p_busco_association[[p]] <- p_busco_association[[p]] + labs(y = NULL)
}

ggsave(
    p_busco_association,
    file = here("products", "figs", "sup_fig1.png"),
    width = 12, height = 5, dpi = 300
)







