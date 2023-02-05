fs::dir_create("output")
files <- fs::dir_ls("../larochelle/v2/output", glob = "*.png")
datas <- fs::dir_ls("../larochelle/v2/output", glob = "*.srda")
dataqs <- fs::dir_ls("../larochelle/v2/output", glob = "*.qs")
datasqs <- fs::dir_ls("../larochelle/v2/output", glob = "*.sqs")

fs::file_copy(c(files, datas, dataqs, datasqs), "output/", overwrite = TRUE)
