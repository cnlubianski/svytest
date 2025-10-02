library(rpms)
library(dplyr)
library(tibble)

svytestCE <- rpms::CE |>
  dplyr::select(-c("INCNONWK", "IRAX", "LIQUIDX", "STOCKX", "STUDNTX")) |> # Mainly NA's
  dplyr::select(-c("VEHQL", "TOBACCCQ", "FOOTWRCQ")) |>  # Not really interesting
  dplyr::select(-PSU) |> # If PSU is available, it is always the same as CID
  dplyr::filter(TOTEXPCQ > 0, FINCBTAX > 0, SALARYX < 250000, SALARYX > 0) |>
  dplyr::select(-EARNER) |> # All are earners at this point so drop column
  tibble::as_tibble()

usethis::use_data(svytestCE, compress = "xz", overwrite = TRUE)

