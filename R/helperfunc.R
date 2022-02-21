#' @importFrom rlang .data

check_align <- function(db_name, 
                        db_path, 
                        task, 
                        acsn_path, 
                        acsn_list, 
                        tab_out, 
                        tab_path, 
                        run_blst, 
                        alt_annot, 
                        alt_path) {
  
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
    } else if (!file.exists(
      here::here(acsn_path, 
                 acsn_list))) {
      stop("Accession list does not exist in the specified path")
    }
  }
  
  if (run_blst & 
      file.exists(
        here::here(tab_path,
                   paste("Alignment_",
                         tab_out,
                         ".csv",
                         sep = "")))) {
    overwrite <- utils::askYesNo(
      "Blast annotations exist in the specified directory. Overwrite?")
  } else {
    overwrite <- TRUE
  }
  
  if (run_blst & overwrite) {
    
    files_to_copy <- c(
      paste(db_path, "/taxdb.bti", sep = ""),
      paste(db_path, "/taxdb.btd", sep = "")
    )
    if (file.exists(files_to_copy[1]) & 
        file.exists(files_to_copy[2])) {
      file.copy(files_to_copy, 
                to = here::here(), 
                overwrite = TRUE)
    } else {
      taxdb_check <- utils::askYesNo(
        paste("No taxdb files found in ", 
              db_path, 
              ". Abort?", 
              sep = ""))
      if (taxdb_check) {
        stop("Aborting: please specify correct directory with taxdb files")
      }
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
               recursive = TRUE)
  }
  
  if (alt_annot) {
    list_dbs <- list.files(alt_path)
    
    to_find <- list("silva.*train", 
                    "silva.*species", 
                    "rdp.*train", 
                    "rdp.*species") %>%
      purrr::map(.f = function(x) {
        check <- str_subset(list_dbs, x) %>%
          purrr::is_empty()
        if (check) {
          stop(paste(x, "not found in the specified directory"))
        }
      })
    
    if (file.exists(
      here::here(tab_path, 
                 "Silva_annot.fst"))) {
      owrt_silva <- utils::askYesNo(
        "Silva annotations exist in the specified directory. Overwrite?"
      )
    } else {
      owrt_silva <- TRUE
    }
    
    if (file.exists(
      here::here(tab_path, 
                 "RDP_annot.fst"))) {
      owrt_RDP <- utils::askYesNo(
        "RDP annotations exist in the specified directory. Overwrite?"
      )
    } else {
      owrt_RDP <- TRUE
    }
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

add_score <- function(x, y, z) {
  
  x_lin <- flatten(x) %>%
    purrr::keep(names(.) %in% c("superkingdom", 
                                "phylum", 
                                "class",
                                "order",
                                "family",
                                "genus",
                                "species"))
  
  silva <- unlist(x_lin)[1:6] %>% 
    stringr::str_count(unlist(y)[1:6]) %>%
    tidyr::replace_na(0) %>%
    sum()
  
  silva_sp <- unlist(x_lin)[7] %>%
    stringr::str_count(unlist(y)[7]) %>%
    tidyr::replace_na(0)*2
  
  
  RDP <- unlist(x_lin)[1:6] %>%
    stringr::str_count(unlist(z)[1:6]) %>%
    tidyr::replace_na(0) %>% 
    sum()
  
  RDP_sp <- unlist(x_lin)[7] %>%
    stringr::str_count(unlist(z)[7]) %>%
    tidyr::replace_na(0)*2
  
  meta <- purrr::flatten(x) %>%
    purrr::keep(
      stringr::str_detect(
        names(.), 
        "PMIDs|host|Isolation_source")) %>%
    purrr::map(.f = ~ifelse(!is.na(.), 1, 0)) %>%
    unlist() %>%
    sum()*3
  
  out <- silva + silva_sp + RDP + RDP_sp + meta
}

meta_extr <- function(x) {
  df_extrt <- x %>%
    dtplyr::lazy_dt() %>%
    dplyr::mutate(
      dplyr::across(
        .cols = tidyselect::vars_select_helpers$where(is.factor), 
        .fns = ~as.character(.x))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(host = trimws(
      stringr::str_extract(
        .data$meta,
        pattern = stringr::regex("(?<=host:).*?(?=\\s-)")
      ))) %>%
    dplyr::mutate(Isolation_source = trimws(
      stringr::str_extract(
        .data$meta,
        pattern = stringr::regex(
          paste("(?<=isolation_source:).*?(?=\\s-)",
                "(?<=isolation_source:).*?(?=$)", 
                sep = "|"))))) %>%
    as.data.frame()
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

