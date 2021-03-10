#' @importFrom rlang .data

Clean_list <- function(x) {
  clean <- x %>%
    dplyr::select(.data$SubType, .data$SubName) %>%
    as.list() %>%
    purrr::map(.f = ~ stringr::str_replace_all(.x, pattern = "\\|", replacement = ";")) %>%
    purrr::map(.f = ~ stringr::str_split(.x, pattern = ";")) %>%
    purrr::pmap(~ paste(.x, .y, sep = ": ")) %>%
    purrr::map(.f = ~ paste(.x, collapse = " - ")) %>%
    unlist()
}

collapse_multi <- function(x) {
  if (stringr::str_detect(x, "\\+")) {
    site_split <- unlist(stringr::str_split(x, pattern = "\\+"))
    site_cmp <- Word_banks %>%
      purrr::keep(names(Word_banks) %in% site_split) %>%
      purrr::flatten() %>%
      purrr::map_chr(.data, .f = ~ paste("\\b", .x, "\\b", sep = "")) %>%
      paste(collapse = "|") %>%
      as.list()
  } else {
    as.list(purrr::pluck(Word_banks, x)) %>%
      purrr::map_chr(.data, .f = ~ paste("\\b", .x, "\\b", sep = "")) %>%
      paste(collapse = "|") %>%
      as.list()
  }
}

Word_banks <- list(
  human = c(
    "human", "humans", "men", "man", "woman", "women", "girl",
    "boy", "children", "child", "clinical sample", "clinical samples",
    "clinical strain", "clinical strains", "medical strain",
    "medical strains", "patient", "patients", "clinical isolate",
    "clinical isolates", "Homo sapiens", "Isolation_source: clinical"
  ),
  FRS = c(
    "vagina", "vaginal", "cervix", "uterus", "cervicovaginal",
    "female genital tract", "vaginitis", "follicular",
    "vulvovaginal candidiasis", "vaginosis", "labia", "amniotic fluid",
    "fornix", "endometrium", "endometrial", "fallopian", "ovary", "ovarian"
  ),
  vagina = c(
    "vagina", "vaginal", "cervix", "cervicovaginal",
    "female genital tract", "vaginitis", "vulvovaginal candidiasis",
    "vaginosis"
  ),
  gut = c(
    "gut", "GI tract", "gastrointestinal", "intestine", "intestinal",
    "urinary tract", "bowel", "abdomen", "abdominal",
    "stomach", "feces", "faeces", "stool", "fecal", "colon", "ileum"
  ),
  skin = c(
    "skin", "scalp", "head", "neck", "dermatitis", "eczema",
    "epidermis", "topical"
  ),
  oral = c(
    "oral cavity", "mouth", "saliva", "salivary glands", "nasal",
    "nasal cavity", "buccal cavity", "oral cavity", "dental",
    "oral bacteria", "respiratory secretion"
  ),
  clinical = c(
    "clinical sample", "clinical samples", "clinical strain",
    "clinical strains", "clinical isolate", "clinical isolates",
    "medical strain", "medical strains", "human sources",
    "Isolation_source: clinical"
  ),
  non_human = c(
    "Veterinary", "Animal clinic", "chicken",
    "Horticulture", "plant", "murine", "mice",
    "conservation science and restoration", "DOVE"
  )
)

