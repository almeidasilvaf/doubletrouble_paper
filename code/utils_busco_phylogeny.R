
#' Get represented species for each BUSCO gene
#'
#' @param outdir Character indicating the path to the BUSCO output directory.
#'
get_busco_representation <- function(outdir) {
    
    species <- dir(outdir)
    species <- species[grepl("\\.fa$", species)]
    
    busco_rep <- lapply(species, function(x) {
        
        seq_dir <- list.dirs(file.path(outdir, x))
        seq_dir <- seq_dir[grep("single_copy_busco_sequences", seq_dir)]
        
        seq_files <- list.files(seq_dir, full.names = TRUE)
        seq_files <- gsub("\\.faa$", "", basename(seq_files))
        
        return(seq_files)
    })
    names(busco_rep) <- gsub("\\.fa", "", species)
    
    
    # Combine sequences from different species, but same BUSCO family
    buscos <- unlist(busco_rep)
    species_per_busco <- split(names(buscos), buscos)
    
    return(species_per_busco)
}


#' Create a list of AAStringSet objects with sequences for BUSCO families
#' 
#' @param outdir Path to directory containing BUSCO results.
#' @param conservation_freq Numeric indicating the minimum relative frequency
#' of conservation across species. Values must range from 0 to 1. Default: 1.
#' 
read_busco_sequences <- function(
        outdir, conservation_freq = 1, verbose = FALSE
) {
    
    species <- dir(outdir)
    species <- species[grepl("\\.fa$", species)]
    
    # Create a vector of BUSCOs to extract
    busco_rep <- get_busco_representation(outdir)
    perc <- round(length(species) * conservation_freq)
    selected_buscos <- names(busco_rep[lengths(busco_rep) >= perc])
    
    # Create a list (A) of lists B) containing seqs of BUSCOs B for species A
    if(verbose) { message("Getting BUSCO sequences...") }
    busco_seqs <- lapply(species, function(x) {
        
        seq_dir <- list.dirs(file.path(outdir, x))
        seq_dir <- seq_dir[grep("single_copy_busco_sequences", seq_dir)]
        
        seq_files <- list.files(seq_dir, full.names = TRUE)
        seq_files_names <- gsub("\\.faa$", "", basename(seq_files))
        seq_files <- seq_files[seq_files_names %in% selected_buscos]
        
        seqs <- lapply(seq_files, function(y) {
            s <- Biostrings::readAAStringSet(y)
            names(s) <- gsub("\\.fa$", "", x)
            
            return(s)
        })
        names(seqs) <- gsub("\\.faa$", "", basename(seq_files))
        
        return(seqs)
    })
    names(busco_seqs) <- gsub("\\.fa", "", species)
    
    # Combine sequences from different species, but same BUSCO family
    if(verbose) { message("Combining sequences from the same BUSCO family...") }
    buscos <- unique(unlist(lapply(busco_seqs, names)))
    busco_comb <- setNames(
        do.call(mapply, c(FUN = c, lapply(busco_seqs, `[`, buscos))), 
        buscos
    )
    busco_comb <- lapply(busco_comb, function(x) Reduce(c, x))
    
    # Remove BUSCOs that are not conserved
    perc <- length(species) * conservation_freq
    busco_comb <- busco_comb[lengths(busco_comb) >= perc]
    
    return(busco_comb)
}


#' Create presence/absence matrix for BUSCO genes in each species
#' 
#' @param busco_list A list of AAStringSet objects with sequences for each
#' BUSCO gene as returned by \code{read_busco_sequences()}.
#' 
#' 
get_busco_pav <- function(busco_list) {
    
    # Get vector of all possible species
    all_species <- sort(unique(unlist(lapply(busco_list, names))))
    
    # Create the PAV matrix
    pav_mat <- Reduce(cbind, lapply(seq_along(busco_list), function(x) {
        
        busco <- names(busco_list)[x]
        included_species <- data.frame(
            row.names = all_species,
            pav = ifelse(all_species %in% names(busco_list[[x]]), 1, 0)
        )
        names(included_species)[1] <- paste0("B_", busco)
        
        return(included_species)
    }))
    
    return(as.matrix(pav_mat))
}

#' Summarize a matrix of BUSCO presence/absence variation
#' 
#' @param busco_pav A presence/absence variations for BUSCO genes as returned
#' by \code{get_busco_pav}.
summarize_busco_pav <- function(busco_pav) {
    
    # Frequency of families by species
    fam_count <- rowSums(busco_pav)
    freq_fam_df <- data.frame(
        row.names = rownames(busco_pav),
        family_count = fam_count,
        family_percentage = round(fam_count / ncol(busco_pav) * 100)
    )
    
    # Frequency of species by families
    species_count <- colSums(busco_pav)
    freq_species_df <- data.frame(
        row.names = colnames(busco_pav),
        species_count = species_count,
        species_percentage = round(species_count / nrow(busco_pav) * 100)
    )
    
    # Combine results into a list
    summary_list <- list(
        families_by_species = freq_fam_df,
        species_by_family = freq_species_df
    )
    
    return(summary_list)
}


#' Simulate species included as a function of BUSCO conservation thresholds
#' 
#' @param busco_pav A presence/absence variations for BUSCO genes as returned
#' by \code{get_busco_pav}.
#' 
simulate_conservation_threshold <- function(
        busco_pav, thresholds = seq(50, 100, by = 5)
) {
    
    perc <- summarize_busco_pav(busco_pav)$species_by_family
    
    # For each threshold, get percentage of included species
    sim <- Reduce(rbind, lapply(thresholds, function(x) {
        
        tokeep <- rownames(perc[perc$species_percentage >= x, ])
        new_pav <- busco_pav[, tokeep]
        
        # Get % of species with at least 1 BUSCO gene
        sums <- rowSums(new_pav)
        pspecies <- length(sums[sums >= 1]) / nrow(busco_pav)
        
        sim_df <- data.frame(
            threshold = x,
            perc_species = round(pspecies * 100),
            n_BUSCOs = length(tokeep)
        )
        
        return(sim_df)
    }))
    
    # Plot simulations
    p1 <- ggplot(sim, aes(y = .data$perc_species, x = .data$threshold)) +
        geom_point(stat = "identity", color = "dodgerblue3") +
        geom_line(color = "dodgerblue4") +
        theme_bw() +
        labs(
            title = "Included species for each BUSCO conservation threshold",
            x = "BUSCO conservation threshold (%)",
            y = "Percentage of included species"
        )
    
    p2 <- ggplot(sim, aes(x = .data$threshold, y = .data$n_BUSCOs)) +
        geom_point(stat = "identity", color = "dodgerblue3") +
        geom_line(color = "dodgerblue4") +
        theme_bw() +
        labs(
            title = "Number of BUSCOs for each conservation threshold",
            x = "BUSCO conservation threshold (%)",
            y = "Number of BUSCO genes"
        )
        
    p_combined <- patchwork::wrap_plots(p1, p2)
    
    return(p_combined)
}



#' Run MAFFT to perform a multiple sequence alignment from the R session
#' 
#' @param seq_list A list of `AAStringSet` objects for each BUSCO family
#' @param threads Numeric indicating the number of threads to use.
#' Default: 2
#' 
align_sequences <- function(seq_list, threads = 2) {
    
    # Write sequences to tempdir()
    seq_paths <- unlist(lapply(seq_along(seq_list), function(x) {
        
        fname <- file.path(tempdir(), paste0(names(seq_list)[x], ".fa"))
        Biostrings::writeXStringSet(
            seq_list[[x]],
            filepath = fname
        )
        
        return(fname)
    }))
    
    # Align sequences
    aln_dir <- file.path(tempdir(), "aligned_seqs")
    dir.create(aln_dir)
    
    aln_seqs <- lapply(seq_paths, function(s) {
        
        outfile <- file.path(aln_dir, gsub("\\.fa", "_aln\\.fa", basename(s)))
        args <- c("--auto --thread", threads, s, " > ", outfile)
        system2("mafft", args = args)
        
        aama <- Biostrings::readAAStringSet(outfile)
        
        return(aama)
    })
    names(aln_seqs) <- gsub("\\.fa", "", basename(seq_paths))
    
    return(aln_seqs)
}


#' Trim MSA in an AAStringSet object
trim_alignment <- function(alignment, max_gap = 0.5) {
    # Calculate the total number of sequences and width of alignment
    total_seqs <- length(alignment)
    alignment_width <- unique(Biostrings::width(alignment))
    
    # Get gap percentages
    aln_mat <- as.matrix(alignment)
    gap_perc <- colSums(aln_mat == "-") / length(alignment)
    
    # Filter alignment
    aln_mat_filt <- aln_mat[, gap_perc <= max_gap]
    
    # Recreate AAStringSet object from matrix
    seqs <- apply(aln_mat_filt, 1, paste0, collapse = "")
    aa_seqs <- AAStringSet(seqs)
    
    return(aa_seqs)
}


#' Concatenate two or more multiple sequence alignments into one
#' 
#' @param aln_seqlist A list of `AAStringSet` objects with multiple
#' sequence alignments for each gene.
#' 
concatenate_alignments <- function(aln_seqlist) {
    
    species <- sort(unique(unlist(lapply(aln_seqlist, names))))
    
    # Create a 'complete' MSA by adding gaps for missing species
    complete_aln <- lapply(seq_along(aln_seqlist), function(x) {
        
        aln <- aln_seqlist[[x]]
        missing <- species[!species %in% names(aln)]
        
        # Create a 'dummy' AAStringSet object with only gaps
        awidth <- unique(Biostrings::width(aln))
        seq_string <- paste0(rep("-", awidth), collapse = "")
        dummy_aa <- AAStringSet(rep(seq_string, length(missing)))
        names(dummy_aa) <- missing
        
        # Combine real MSA with dummy alignment
        final_aa <- c(aln, dummy_aa)[species]
        
        return(final_aa)
    })
    
    # Concatenate alignments
    conc_aln <- Reduce(xscat, complete_aln)
    names(conc_aln) <- species
    
    return(conc_aln)
}


#' Infer a species tree from a multiple sequence alignment with IQ-TREE2
#' 
#' @param aln_seqlist A list of `AAStringSet` objects.
#' 
infer_species_tree <- function(aln_seqlist, outgroup, threads = 2) {
    
    trees <- lapply(seq_along(aln_seqlist), function(x) {
        
        message("Inferring tree for ", names(aln_seqlist)[x])
        
        # Write temporary file
        aln_file <- tempfile(fileext = ".fa")
        
        Biostrings::writeXStringSet(aln_seqlist[[x]], filepath = aln_file)
        
        # Run IQ-TREE2
        rgs <- c("-s", aln_file, "-o", outgroup, "-B 1000 -T", threads)
        system2("iqtree2", args = rgs)
        tree <- ape::read.tree(paste0(aln_file, ".treefile"))
        
        return(tree)
    })
    names(trees) <- names(aln_seqlist)
    
    return(trees)
}
