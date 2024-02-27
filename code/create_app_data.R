
library(here)
library(tidyverse)
library(doubletrouble)

# metadata_all.rda ----
load(here("products", "result_files", "metadata_all.rda"))
names(metadata_all) <- c("Fungi", "Plants", "Metazoa", "Protists", "Vertebrates")
metadata_all$Plants <- metadata_all$Plants |>
    dplyr::filter(species != "triticum_aestivum_lancer")

save(
    metadata_all, 
    file = "~/Dropbox/package_benchmarks/doubletroubledb/data/metadata_all.rda",
    compress = "xz"
)


# trees.rda ----
load(here("products", "result_files", "trees", "fungi_busco_trees.rda"))
load(here("products", "result_files", "trees", "metazoa_busco_trees.rda"))
load(here("products", "result_files", "trees", "plants_busco_trees.rda"))
load(here("products", "result_files", "trees", "protists_busco_trees.rda"))
load(here("products", "result_files", "trees", "vertebrates_busco_trees.rda"))

## Rename tip labels of trees
tree_fungi <- fungi_busco_trees$conc
tree_fungi$tip.label <- gsub("\\.", "_", tree_fungi$tip.label)

tree_protists <- protists_busco_trees$conc
tree_protists$tip.label <- gsub("\\.", "_", tree_protists$tip.label)

tree_plants <- plants_busco_trees$conc
tree_plants$tip.label <- gsub("\\.", "_", tree_plants$tip.label)

tree_metazoa <- metazoa_busco_trees$conc
tree_metazoa$tip.label <- gsub("\\.", "_", tree_metazoa$tip.label)

tree_vertebrates <- vertebrates_busco_trees$conc
tree_vertebrates$tip.label <- gsub("\\.", "_", tree_vertebrates$tip.label)

trees <- list(
    Fungi = tree_fungi,
    Plants = tree_plants,
    Metazoa = tree_metazoa,
    Protists = tree_protists,
    Vertebrates = tree_vertebrates
)

save(
    trees, 
    file = "~/Dropbox/package_benchmarks/doubletroubledb/data/trees.rda",
    compress = "xz"
)

# dup_counts.rda ----
load(here("products", "result_files", "fungi_duplicates_unique.rda"))
load(here("products", "result_files", "protists_duplicates_unique.rda"))
load(here("products", "result_files", "plants_duplicates_unique.rda"))
load(here("products", "result_files", "vertebrates_duplicates_unique.rda"))
load(here("products", "result_files", "metazoa_duplicates_unique.rda"))

counts_fungi <- duplicates2counts(fungi_duplicates_unique)
counts_protists <- duplicates2counts(protists_duplicates_unique)
counts_plants <- duplicates2counts(plants_duplicates_unique)
counts_vertebrates <- duplicates2counts(vertebrates_duplicates_unique)
counts_metazoa <- duplicates2counts(metazoa_duplicates_unique)

dup_counts <- list(
    Fungi = counts_fungi,
    Plants = counts_plants,
    Metazoa = counts_metazoa,
    Protists = counts_protists,
    Vertebrates = counts_vertebrates
)

save(
    dup_counts,
    file = "~/Dropbox/package_benchmarks/doubletroubledb/data/dup_counts.rda",
    compress = "xz"
)


# Zenodo data: duplicated gene pairs and genes for each species ----
## Create temporary folder with data
zenodo_dir <- "~/Documents/doubletrouble_zenodo"
if(!dir.exists(zenodo_dir)) { dir.create(zenodo_dir, recursive = TRUE) }

# Define wrapper functions to write lists to files
write_dups <- function(dlist, outdir, type = c("genes", "pairs")) {
    
    w <- lapply(seq_along(dlist), function(x) {
        
        fname <- file.path(outdir, paste0(names(dlist)[x], "_pairs.tsv.gz"))
        if(type == "genes") { 
            fname <- file.path(outdir, paste0(names(dlist)[x], "_genes.tsv.gz")) 
        }
        
        readr::write_tsv(dlist[[x]], fname)
        
        return(NULL)
    })
    
    return(dir(outdir))
}

## Write duplicate pairs
load(here("products", "result_files", "fungi_duplicates.rda"))
load(here("products", "result_files", "protists_duplicates.rda"))
load(here("products", "result_files", "plants_duplicates.rda"))
load(here("products", "result_files", "vertebrates_duplicates.rda"))
load(here("products", "result_files", "metazoa_duplicates.rda"))

write_dups(fungi_duplicates, outdir = zenodo_dir, type = "pairs")
write_dups(protists_duplicates, outdir = zenodo_dir, type = "pairs")
write_dups(plants_duplicates, outdir = zenodo_dir, type = "pairs")
write_dups(vertebrates_duplicates, outdir = zenodo_dir, type = "pairs")
write_dups(metazoa_duplicates, outdir = zenodo_dir, type = "pairs")



## Write duplicated genes
load(here("products", "result_files", "fungi_duplicates_unique.rda"))
load(here("products", "result_files", "protists_duplicates_unique.rda"))
load(here("products", "result_files", "plants_duplicates_unique.rda"))
load(here("products", "result_files", "vertebrates_duplicates_unique.rda"))
load(here("products", "result_files", "metazoa_duplicates_unique.rda"))

write_dups(fungi_duplicates_unique, outdir = zenodo_dir, type = "genes")
write_dups(protists_duplicates_unique, outdir = zenodo_dir, type = "genes")
write_dups(plants_duplicates_unique, outdir = zenodo_dir, type = "genes")
write_dups(vertebrates_duplicates_unique, outdir = zenodo_dir, type = "genes")
write_dups(metazoa_duplicates_unique, outdir = zenodo_dir, type = "genes")

## Note: .collinearity files were manually moved to the directory

## Create .zip files for each species
files <- list.files(zenodo_dir, full.names = TRUE)

file_df <- data.frame(
    file = files,
    species = basename(files)
) |>
    mutate(
        species = str_replace_all(species, "_genes.*", ""),
        species = str_replace_all(species, "_pairs.*", ""),
        species = str_replace_all(species, "\\.collinearity.*", ""),
        species = str_replace_all(species, "\\.", "_")
    )

file_list <- split(file_df$file, file_df$species)

zips <- lapply(seq_along(file_list), function(x) {
    
    fname <- file.path(zenodo_dir, paste0(names(file_list)[x], ".zip"))
    fvec <- paste(file_list[[x]], collapse = " ")
    z <- system2("zip", args = c(fname, fvec))
    
    return(NULL)
})

## Move .zip files to subdirectories for each Instance
load(here("products", "result_files", "metadata_all.rda"))
names(metadata_all) <- c("Fungi", "Plants", "Metazoa", "Protists", "Vertebrates")

zip_files <- list.files(zenodo_dir, pattern = ".zip", full.names = TRUE)

dir.create(file.path(zenodo_dir, "Plants"))
zip_plants <- zip_files[gsub("\\.zip", "", basename(zip_files)) %in% metadata_all$Plants$species]

dir.create(file.path(zenodo_dir, "Fungi"))
zip_fungi <- zip_files[gsub("\\.zip", "", basename(zip_files)) %in% metadata_all$Fungi$species]

dir.create(file.path(zenodo_dir, "Metazoa"))
zip_metazoa <- zip_files[gsub("\\.zip", "", basename(zip_files)) %in% metadata_all$Metazoa$species]

dir.create(file.path(zenodo_dir, "Protists"))
zip_protists <- zip_files[gsub("\\.zip", "", basename(zip_files)) %in% metadata_all$Protists$species]

dir.create(file.path(zenodo_dir, "Vertebrates"))
zip_vertebrates <- zip_files[gsub("\\.zip", "", basename(zip_files)) %in% metadata_all$Vertebrates$species]
zip_vertebrates <- zip_vertebrates[file.exists(zip_vertebrates)]

fs::file_move(zip_plants, new_path = file.path(zenodo_dir, "Plants"))
fs::file_move(zip_fungi, new_path = file.path(zenodo_dir, "Fungi"))
fs::file_move(zip_metazoa, new_path = file.path(zenodo_dir, "Metazoa"))
fs::file_move(zip_protists, new_path = file.path(zenodo_dir, "Protists"))
fs::file_move(zip_vertebrates, new_path = file.path(zenodo_dir, "Vertebrates"))

