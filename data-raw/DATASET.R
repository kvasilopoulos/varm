## code to prepare `DATASET` dataset goes here

yraw <- readr::read_tsv("data-raw/Yraw.dat", col_names = FALSE)

usethis::use_data(yraw, overwrite = TRUE)
