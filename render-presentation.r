library(purrr)
library(stringr)

output <- Sys.getenv("QUARTO_PROJECT_OUTPUT_DIR")
pres <- c("MEAPS_p2203_fr")
walk(pres, ~{
  quarto::quarto_render("MEAPS_p2203_fr.qmd")
  file.copy(str_c(.x, ".html"), output, overwrite=TRUE)
  file.copy(str_c(.x, "_files"), output, overwrite=TRUE, recursive=TRUE)
  file.remove(str_c(.x, ".html"))
  unlink(str_c(.x, "_files"), recursive=TRUE)
    })
file.copy("images/", output, recursive = TRUE)
