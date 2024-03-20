bibs <- fs::dir_ls(glob="*.bib")
fs::file_copy(bibs, "_presentations/", overwrite = TRUE)
