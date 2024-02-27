
library(here)
library(tidyverse)

# Supplementary Table S1: species metadata ----
load(here("products", "result_files", "metadata_all.rda"))
table_s1 <- Reduce(rbind, metadata_all)

readr::write_tsv(
    table_s1,
    file = here("products", "tables", "sup_table_S1.tsv")
)

# Supplementary Table S2: SD-rich species ----
load(here("products", "result_files", "sd_abundant_spp.rda"))
table_s2 <- sd_abundant_spp |>
    dplyr::select(
        Species = species,
        Instance = instance,
        SD_count = n,
        SD_percentage = percentage
    )

readr::write_tsv(
    table_s2,
    file = here("products", "tables", "sup_table_S2.tsv")
)


# Making Excel file
s1 <- read_tsv(here("products", "tables", "sup_table_S1.tsv"))
s2 <- read_tsv(here("products", "tables", "sup_table_S2.tsv"))

openxlsx::write.xlsx(
    x = list(S1 = s1, S2 = s2),
    file = here("products", "tables", "sup_tables.xlsx")
)

