fs::dir_create("meaps-doc/output")
files <- fs::dir_ls("v2/output", glob = "*.png")
datas <- fs::dir_ls("v2/output", glob = "*.srda")
fs::file_copy(c(files, datas), "meaps-doc/output/", overwrite = TRUE)
