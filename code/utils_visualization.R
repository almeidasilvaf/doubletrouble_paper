
#' Wrapper function to plot tree with taxa highlighted
plot_tree_taxa <- function(
        tree, metadata, taxon, min_n = 2, 
        text_x = -2, text_size = 3, padding_text = 1, min_n_lab = 3,
        type = "rect", taxon_levels = NULL
) {
    
    # Filter metadata and clean taxon variable
    tokeep <- count(metadata, .data[[taxon]]) |> filter(n >= min_n) |> 
        pull(.data[[taxon]])
    
    tokeep2 <- count(metadata, .data[[taxon]]) |> filter(n >= min_n_lab) |> 
        pull(.data[[taxon]])
    
    meta <- metadata |>
        filter(species %in% tree$tip.label) |> 
        mutate(
            Taxon = ifelse(.data[[taxon]] %in% tokeep, .data[[taxon]], "Other"),
            Taxon2 = ifelse(.data[[taxon]] %in% tokeep2, .data[[taxon]], "Other")
        ) |>
        select(all_of(c("species", "Taxon", "Taxon2")))
    
    ## ## Plot species tree tips colored by taxon
    p_tree <- ggtree(tree, branch.length = "none") %<+% meta
    
    if(type == "tippoint") {
        p_tree <- p_tree +
            geom_tippoint(aes(color = .data$Taxon)) +
            ggsci::scale_color_d3("category20") +
            theme(legend.position = "left") +
            labs(color = str_to_title(taxon))
    } else {
        ## Get rectangle coordinates
        rect_data <- p_tree$data |>
            filter(isTip, Taxon != "Other")
        
        # Get text coordinates
        lev <- taxon_levels
        if(is.null(taxon_levels)) { lev <- unique(rect_data$Taxon2) }
        
        text_data <- rect_data |>
            dplyr::filter(Taxon2 %in% lev) |>
            dplyr::filter(Taxon2 != "Other") |>
            group_by(Taxon2) |>
            summarise(y = max(y) - padding_text) |>
            ungroup() |>
            mutate(x = text_x)
        
        p_tree <- p_tree + 
            geom_rect(
                data = rect_data,
                aes(
                    xmin = -Inf, 
                    xmax = Inf,
                    ymin = ifelse(isTip, y - 0.5, NA),
                    ymax = ifelse(isTip, y + 0.5, NA),
                    fill = Taxon
                ),
                alpha = 0.3,
                show.legend = FALSE
            ) +
            geom_label(
                data = text_data,
                aes(x = x, y = y, label = Taxon2),
                hjust = 0, size = text_size, label.size = NA
            ) +
            ggsci::scale_fill_d3("category20")
        
    }
    
    return(p_tree)
}
