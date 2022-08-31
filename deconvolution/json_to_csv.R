library(jsonlite)
library(tidyverse)

json_full <- readLines("clean_results.json")
parsed_full <- lapply(json_full, fromJSON)

jsonlist_to_row <- function(x) {
  if (isa(x$atc_codes, "list")) {
    atc <- x$atc_codes
  } else if (length(x$atc_codes == 2)) {
    atc <-  list(code1 = x$atc_codes[1, ],
                 code2 = x$atc_codes[2, ])
  }
  data.frame(
    X1 = x$table[1, 1],
    X2 = x$table[2, 1],
    N1 = x$table[1, 2],
    N2 = x$table[2, 2],
    unadjusted_Z = x$cox$unadjusted$coef,
    unadjusted_se = x$cox$unadjusted$se,
    unadjusted_p = x$cox$unadjusted$p,
    ipw_Z = x$cox$logistic_match$coef,
    ipw_se = x$cox$logistic_match$se,
    ipw_p = x$cox$logistic_match$p,
    trial = paste0(x$sub_infos$study, collapse = "-"),
    icd10 = paste0(x$icd10codes, collapse = "-"),
    atc1 = paste0(atc[[1]], collapse = "-"),
    atc2 = paste0(atc[[2]], collapse = "-")
  )
}

tbl_rows <- lapply(parsed_full, jsonlist_to_row)
tbl_full_df <- bind_rows(tbl_rows)
write_csv(tbl_full_df, "trialverify_full.csv")