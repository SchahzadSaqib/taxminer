#' @importFrom rlang .data

check_align <- function(db_name, 
                        db_path, 
                        task, 
                        acsn_path, 
                        acsn_list, 
                        tab_out, 
                        tab_path) {
  
  if (is.null(db_path)) {
    stop("Please provide the full path to the database")
  } else if (!dir.exists(db_path)) {
    stop(paste("The directory: ", 
               db_path, 
               " does not exist", 
               sep = ""))
  }
  
  if (is.null(db_name)) {
    stop("Please provide the full name of the database")
  }
  
  if (!is.null(acsn_list)) { 
    if (is.null(acsn_path)) { 
      stop("Please provide a valid path to the accession ID list")
    } else if (!file.exists(here::here(acsn_path, 
                                       acsn_list))) {
      stop("Accession list does not exist in the specified path")
    }
  }
  
  task_list <- c("megablast", 
                 "blastn", 
                 "blastn-short", 
                 "dc-megablast")
  if (!task %in% task_list) {
    stop(paste("Please specify one of the following tasks: ", 
               paste(task_list, collapse = "; "),  
               sep = ""))
  }
  
  if (!dir.exists(tab_path)) {
    dir.create(tab_path, 
               recursive = T)
  }
}


check_ecosrc <- function(hit_tbl, 
                         filt_host, 
                         filt_site, 
                         filt_negt) {
  
  if (nrow(hit_tbl) == 0) {
    stop("Empty Data.frame - Aborting Relevance filtration")
  }
  
  filt_split <- unlist(
    stringr::str_split(c(filt_host, 
                         filt_site, 
                         filt_negt), 
                       pattern = "\\+")
  )
  
  filt_split <- filt_split[!is.na(filt_split)]
  
  wb_step <- 1
  while(wb_step <= length(filt_split)) {
    if (!filt_split[wb_step] %in% names(Word_banks)) {
      stop(paste("Invalid filtration entry ", 
                 filt_split[wb_step], 
                 sep = ""))
    } else {
      wb_step <- wb_step + 1
    }
  }
}

check_lineage <- function(taxids, 
                          bindtoAcc, 
                          asbd_tbl) {
  
  if (!any(names(taxids) %in% "TaxID")) {
    stop(paste("No column named 'TaxID' -", 
               "please provide a valid input column", 
               sep = ""))
  }
  
  if (bindtoAcc & !any(names(taxids) %in% "AccID")) {
    stop(paste("No columns named 'AccID' while bindtoAcc is T - ",
               "please provide a valid accession ID column"))
  }
  
  if (!is.na(asbd_tbl)) {
    if (!file.exists(asbd_tbl)) {
      stop(paste("Pre-assembled table does not exist in this location: ", 
                 asbd_tbl, 
                 sep = ""))
    }
  }
}


Clean_list <- function(x) {
  clean <- x %>%
    dplyr::select(.data$SubType, 
                  .data$SubName) %>%
    as.list() %>%
    purrr::map(
      .f = ~ stringr::str_replace_all(
        .x, pattern = "\\|", 
        replacement = ";")
      ) %>%
    purrr::map(
      .f = ~ stringr::str_split(
        .x, pattern = ";")
      ) %>%
    purrr::pmap(~ paste(.x, 
                        .y, 
                        sep = ": ")) %>%
    purrr::map(.f = ~ paste(.x, 
                            collapse = " - ")) %>%
    unlist()
}


collapse_multi <- function(x) {
  if (stringr::str_detect(x, "\\+")) {
    site_split <- unlist(
      stringr::str_split(x, 
                         pattern = "\\+"))
    site_cmp <- Word_banks %>%
      purrr::keep(names(Word_banks) %in% site_split) %>%
      purrr::flatten() %>%
      purrr::map_chr(.data, 
                     .f = ~ paste("\\b", 
                                  .x, 
                                  "\\b", 
                                  sep = "")) %>%
      paste(collapse = "|") %>%
      as.list()
  } else {
    as.list(purrr::pluck(Word_banks, x)) %>%
      purrr::map_chr(.data, 
                     .f = ~ paste("\\b", 
                                  .x, "\\b", 
                                  sep = "")) %>%
      paste(collapse = "|") %>%
      as.list()
  }
}


#' Word banks
#'
#' Full list of word banks being used for text mining based filtration
#'
#' @export
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
    "urinary tract", "urine", "bowel", "abdomen", "abdominal",
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
    "medical strain", "medical strains", "human sources", "blood",
    "Isolation_source: clinical"
  ),
  non_human = c(
    "Veterinary", "Animal clinic", "chicken", "man-made", "human-made",
    "Horticulture", "plant", "murine", "mice",
    "conservation science and restoration", "DOVE"
  )
)

