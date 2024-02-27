
#' Format a `classification` object (from taxize) into a data frame
format_classification <- function(class_list, verbose = FALSE) {
    
    class_df <- Reduce(rbind, lapply(seq_along(class_list), function(x) {
        
        if(verbose) {
            message("Creating taxonomy table for species ", x, "/", length(class_list))
        }

        cdf <- class_list[[x]]
        
        fn <- function(x) { return(ifelse(length(x) == 0, NA, x)) } 
        
        df <- data.frame(
            ncbi_species = fn(cdf$name[cdf$rank == "species"]),
            family = fn(cdf$name[cdf$rank == "family"]),
            order = fn(cdf$name[cdf$rank == "order"]),
            class = fn(cdf$name[cdf$rank == "class"]),
            superclass = fn(cdf$name[cdf$rank == "superclass"]),
            subphylum = fn(cdf$name[cdf$rank == "subphylum"]),
            phylum = fn(cdf$name[cdf$rank == "phylum"])
        )
        
        return(df)
    }))
    
    return(class_df)
}


#' Wrapper to get proteomes from Ensembl Genomes
get_proteomes <- function(species_metadata, ensembl = "plants") {
    
    species <- gsub("\\.", "_", species_metadata$species)
    
    # Get URL to proteome FTP directories
    burl <- "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-57/"
    url_seq <- paste0(burl, ensembl, "/fasta/", species, "/pep/")
    
    if(ensembl == "ensembl") {
        url_seq <- file.path(
            "https://ftp.ensembl.org/pub/release-110/fasta", species, "pep/"
        )
    }
    
    # Get filtered sequences (primary transcripts only)
    seqs <- lapply(url_seq, function(x) {
        
        # Get links to FASTA files
        seq_files <- XML::getHTMLLinks(RCurl::getURL(x, dirlistonly = TRUE))
        seq_files <- seq_files[grep("\\.fa\\.gz$", seq_files)]
        seq_files <- seq_files[!grepl("abinitio", seq_files)]
        seq_path <- paste0(x, seq_files)
        
        # Get FASTA file as an AAStringSet object
        s <- Biostrings::readAAStringSet(seq_path)
        names(s) <- gsub(" .*", "", gsub(".*gene:", "", names(s)))
        
        # Keep only longest transcript
        s <- s[order(Biostrings::width(s), decreasing = TRUE),]
        s <- s[!duplicated(names(s)), ]
        return(s)
    })
    names(seqs) <- gsub("_", ".", species)
    return(seqs)
}

#' Wrapper function to get CDS from Ensembl and Ensembl Genomes
get_cds_ensembl <- function(species, ensembl = "plants") {
    
    # Get URL to proteome FTP directories
    burl <- "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-57/"
    url_seq <- paste0(burl, ensembl, "/fasta/", species, "/cds/")
    
    if(ensembl == "ensembl") {
        url_seq <- file.path(
            "https://ftp.ensembl.org/pub/release-110/fasta", species, "cds/"
        )
    }
    
    # Get filtered sequences (primary transcripts only)
    seqs <- lapply(url_seq, function(x) {
        
        # Get links to FASTA files
        seq_files <- XML::getHTMLLinks(RCurl::getURL(x, dirlistonly = TRUE))
        seq_files <- seq_files[grep("\\.fa\\.gz$", seq_files)]
        seq_files <- seq_files[!grepl("abinitio", seq_files)]
        seq_path <- paste0(x, seq_files)
        
        # Get FASTA file as an AAStringSet object
        s <- Biostrings::readDNAStringSet(seq_path)
        names(s) <- gsub(" .*", "", gsub(".*gene:", "", names(s)))
        
        # Keep only longest transcript
        s <- s[order(Biostrings::width(s), decreasing = TRUE),]
        s <- s[!duplicated(names(s)), ]
        return(s)
    })
    
    names(seqs) <- species
    return(seqs)
}


#' Wrapper to get gene annotation from Ensembl Genomes
#' 
#' @param ensembl Character indicating the Ensembl/Ensembl Genomes instance
#' to download from. One of 'plants', 'fungi', 'metazoa', 'protists'
#' or 'ensembl'.
get_annotation <- function(species_metadata, ensembl = "plants") {
    
    # Wrapper to parse file sizes
    convb <- function(x) {
        ptn <- "(\\d*(.\\d+)*)(.*)"
        num  <- as.numeric(sub(ptn, "\\1", x))
        unit <- sub(ptn, "\\3", x)             
        unit[unit == ""] <- "1" 
        
        mult <- c("1" = 1, "K" = 1024, "M" = 1024^2, "G" = 1024^3)
        return(num * unname(mult[unit]))
    }
    
    species <- gsub("\\.", "_", species_metadata$species)
    
    # Get URL to proteomes and annotation FTP directories
    burl <- "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-57/"
    url_annot <- paste0(burl, ensembl, "/gff3/", species, "/")
    
    if(ensembl == "ensembl") {
        url_annot <- paste0(file.path(
            "https://ftp.ensembl.org/pub/release-110/gff3", species
        ), "/")
    }
    
    # Get gene ranges
    ranges <- lapply(url_annot, function(x) {
        
        # Get links to GFF files
        annot_files <- XML::readHTMLTable(
            RCurl::getURL(x, dirlistonly = TRUE)
        )[[1]]
        annot_files <- annot_files[endsWith(annot_files$Name, ".gff3.gz"), ] |>
            tidyr::drop_na()
        annot_files <- annot_files[!grepl("abinitio", annot_files$Name), ]
        annot_files$Size <- convb(annot_files$Size)
        annot_files <- annot_files[order(annot_files$Size, decreasing = TRUE), ]
        annot_files <- annot_files$Name[1]
        
        annot_path <- paste0(x, annot_files)
        
        # Get GFF files as GRanges objects
        outfile <- file.path(tempdir(), basename(annot_path))
        tfile <- tempfile(fileext = ".gff3")
        d <- download.file(annot_path, destfile = outfile)
        s <- system2("zcat", args = c(outfile, " | sed '/^##/d' > ", tfile))
        
        granges <- rtracklayer::import(tfile)
        
        # Remove temp files
        unlink(outfile)
        unlink(tfile)
        
        return(granges)
    })
    names(ranges) <- gsub("_", ".", species)
    return(ranges)
}


#' Wrapper function to only keep sequences of genes in annotation object
filter_sequences <- function(proteomes, annotation) {
    
    fseq <- lapply(seq_along(proteomes), function(x) {
        species <- names(proteomes)[x]
        annot_gene <- annotation[[species]][annotation[[species]]$type == "gene"]
        gene_ids <- as.character(annot_gene$gene_id)
        seq_ids <- as.character(names(proteomes[[species]]))
        
        if(!any(seq_ids %in% gene_ids)) {
            seq_ids <- gsub("\\.[0-9]*", "", seq_ids)
            names(proteomes[[species]]) <- seq_ids
        }
        
        seqs <- proteomes[[species]][seq_ids %in% gene_ids]
        return(seqs)
    })
    names(fseq) <- names(proteomes)
    return(fseq)
}


#' Wrapper to download genomes from Ensembl Genomes
download_filtered_proteomes <- function(
        species_metadata, ensembl = "plants", outdir
) {
    
    d <- unlist(lapply(seq_len(nrow(species_metadata)), function(x) {
        
        # Get filtered sequences
        proteome <- get_proteomes(species_metadata[x, ], ensembl = ensembl)
        annot <- get_annotation(species_metadata[x, ], ensembl = ensembl)
        proteome <- filter_sequences(proteome, annot)
        
        fname <- paste0(file.path(outdir, names(proteome)), ".fa")
        
        x <- Biostrings::writeXStringSet(
            proteome[[1]], filepath = fname
        )
        
        return(x)
    }))
    
    return(d)
}



#' Wrapper function to identify and classify duplicated genes from Ensembl 
#' and Ensembl Genomes one by one
ensembl2duplicates <- function(
        metadata, ensembl = NULL, outgroups = NULL,
        collinearity_dir = NULL, threads = NULL,
        tsv_dir = NULL, outdir_base = NULL
) {
    
    idx <- seq_len(nrow(metadata))
    
    # Get processed data for outgroups
    if(!is.null(outgroups)) {
        out <- unique(outgroups[, 2])
        out_meta <- metadata[metadata$species %in% out, ]
        
        out_annot <- get_annotation(out_meta, ensembl = ensembl)
        out_seq <- get_proteomes(out_meta, ensembl = ensembl)
        out_seq <- filter_sequences(out_seq, out_annot)
    }
    
    dup_list <- lapply(idx, function(x) {
        
        sp <- metadata$species[x]
        message("Working on species ", sp)
        
        # Obtaining annotation
        annotation <- get_annotation(metadata[x, ], ensembl = ensembl)
        
        # Obtaining proteomes
        proteomes <- get_proteomes(metadata[x, ], ensembl = ensembl)
        proteomes <- filter_sequences(proteomes, annotation)
        
        if(!is.null(outgroups)) {
            if(sp %in% outgroups[, 1]) {
                proteomes <- c(proteomes, out_seq)
                annotation <- c(annotation, out_annot)
                
                annotation <- annotation[!duplicated(annotation)]
                proteomes <- proteomes[!duplicated(proteomes)]
            }
        }
        
        # Process input
        processed_data <- process_input(
            proteomes, annotation, filter_annotation = TRUE
        )
        
        # Intraspecies DIAMOND searches
        bout <- ifelse(is.null(outdir_base), tempdir(), outdir_base)
        outdir <- file.path(
            bout, 
            paste0(sp, "_intra_", format(Sys.time(), "%d_%b_%Y_%Hh%M"))
        )
        diamond_intra <- run_diamond(
            seq = processed_data$seq[1],
            compare = "intraspecies",
            outdir = outdir,
            ... = "--sensitive",
            threads = threads
        )
        
        # Perform interspecies DIAMOND search in case outgroup exists
        scheme <- "standard"
        diamond_inter <- NULL
        intron_counts <- NULL
        if(!is.null(outgroups)) {
            if(sp %in% outgroups[, 1]) {
                
                out_df <- outgroups[outgroups[, 1] %in% sp, ]
                out_df[, 1] <- gsub("_", "\\.", out_df[, 1])
                out_df[, 2] <- gsub("_", "\\.", out_df[, 2])
                
                inter_dir <- file.path(
                    bout, 
                    paste0(sp, "_outgroup_", format(Sys.time(), "%d_%b_%Y_%Hh%M"))
                )
                diamond_inter <- run_diamond(
                    seq = processed_data$seq,
                    compare = out_df,
                    outdir = inter_dir,
                    ... = "--sensitive",
                    threads = threads
                )
                
                scheme <- "full"
                txdb_list <- lapply(annotation[1], GenomicFeatures::makeTxDbFromGRanges)
                intron_counts <- lapply(txdb_list, get_intron_counts)
                
                delete <- fs::dir_delete(inter_dir)
            }
        }
        
        # Classify duplicate pairs using the extended mode
        duplicate_pairs <- classify_gene_pairs(
            blast_list = diamond_intra,
            annotation = processed_data$annotation,
            blast_inter = diamond_inter,
            intron_counts = intron_counts,
            scheme = scheme,
            collinearity_dir = collinearity_dir
        )[[1]]
        
        if(!is.null(tsv_dir)) {
            if(!dir.exists(tsv_dir)) { dir.create(tsv_dir, recursive = TRUE) }
            fname <- file.path(tsv_dir, paste0(sp, ".tsv.gz"))
            w <- readr::write_tsv(duplicate_pairs, fname)
        }
        
        delete <- fs::dir_delete(outdir)
        return(duplicate_pairs)
    })
    names(dup_list) <- metadata$species
    
    return(dup_list)
}


#' Get a duplicate count matrix for each genome
#'
#' @param duplicate_list A list of data frames with the duplicated genes or
#' gene pairs and their modes of duplication as returned 
#' by \code{classify_gene_pairs()} or \code{classify_genes()}.
#' @param shape Character specifying the shape of the output data frame.
#' One of "long" (data frame in the long shape, in the tidyverse sense),
#' or "wide" (data frame in the wide shape, in the tidyverse sense).
#' Default: "long".
#' @param drop_zero_cols Logical indicating whether to drop columns with zero
#' counts in all rows. Default: TRUE.
#' @param drop_zero_rows Logical indicating whether to drop rows with zero
#' counts in all columns. Default: TRUE.
#' 
#' @return A count matrix containing the frequency of duplicated 
#' genes (or gene pairs) by mode for each speces, 
#' with species in rows and duplication modes in columns.
#' 
#' @examples
#' data(scerevisiae_kaks)
#' 
#' # Get unique duplicates
#' duplicate_list <- classify_genes(scerevisiae_kaks)
#'
duplicates2counts <- function(duplicate_list, shape = "long") {
    
    # Get factor levels for variable `type`
    tlevels <- lapply(duplicate_list, function(x) return(levels(x$type)))
    tlevels <- tlevels[[names(sort(lengths(tlevels), decreasing = TRUE)[1])]]
    
    counts <- Reduce(rbind, lapply(seq_along(duplicate_list), function(x) {
        
        species <- names(duplicate_list)[x]
        
        dup_table <- duplicate_list[[x]]
        dup_table$type <- factor(dup_table$type, levels = tlevels)
        
        if(shape == "long") {
            final_dups <- as.data.frame(table(dup_table$type))
            names(final_dups) <- c("type", "n")
            final_dups$species <- species
        } else if(shape == "wide") {
            final_dups <- t(as.matrix(table(dup_table$type)))
            final_dups <- cbind(species, as.data.frame(final_dups))
        } else {
            stop("Argument 'format' must be one of 'long' or 'wide'.")
        }
        
        return(final_dups)
    }))
    
    return(counts)
}





