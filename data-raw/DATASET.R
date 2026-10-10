## code to prepare `yraw` dataset
yraw <- utils::read.table("data-raw/Yraw.dat", col.names = c("infl", "unemp", "rate"))
save(yraw, file = "data/yraw.rda", compress = "bzip2", version = 2)
