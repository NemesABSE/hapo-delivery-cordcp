# Data location. The HAPO substudy dataset is controlled-access and is not part
# of this repository. Point HAPO_DATA at the directory that holds the three CSV
# files listed in the README, then source this file before running any script:
#
#   Sys.setenv(HAPO_DATA = "/path/to/hapo/csv")
#   source("config.R")

HAPO_DATA <- Sys.getenv("HAPO_DATA", unset = NA)

hapo_file <- function(name) {
  if (is.na(HAPO_DATA) || !nzchar(HAPO_DATA)) {
    stop("HAPO_DATA is not set. See the README: it must point at the directory\n",
         "holding the controlled-access HAPO substudy CSV files.", call. = FALSE)
  }
  path <- file.path(HAPO_DATA, name)
  if (!file.exists(path)) {
    stop("expected data file not found: ", path, call. = FALSE)
  }
  path
}
