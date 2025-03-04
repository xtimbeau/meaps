bibs <- fs::dir_ls(glob="*.bib")
fs::file_copy(bibs, "_presentations/", overwrite = TRUE)
fs::file_copy("r/rinit.r", "_presentations/rinit.r", overwrite = TRUE)
fs::dir_copy("secrets", "_presentations/secrets", overwrite = TRUE)
fs::dir_copy("output", "_presentations/output", overwrite = TRUE)
fs::dir_copy("_templates", "_presentations/_templates", overwrite = TRUE)
fs::dir_copy("images", "_presentations/images", overwrite = TRUE)
