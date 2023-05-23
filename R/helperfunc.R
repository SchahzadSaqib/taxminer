utils::globalVariables(c(
  "tab_path",
  "tab_out",
  "prgrs_bar"
))

#' @importFrom rlang .data


# Word banks ----
# Word banks
#
# Full list of word banks being used for text mining based filtration
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
    "bowel", "abdomen", "abdominal", "stomach", "feces", "faeces",
    "stool", "fecal", "colon", "ileum"
  ),
  ut = c("urine", "urinary tract", "urinary", "bladder", "sphincter"),
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
    "human origins", "Human milk", "Isolation_source: clinical"
  ),
  non_human = c(
    "Veterinary", "Animal clinic", "chicken", "man-made", "human-made",
    "Horticulture", "plant", "murine", "mice", "pig", "equine",
    "conservation science and restoration", "DOVE"
  )
)


# Input integrity checks ----
check_align <- function(db_name,
                        db_path,
                        task,
                        acsn_path,
                        acsn_list,
                        tab_out,
                        tab_path,
                        run_blst,
                        acsn_check,
                        alt_annot,
                        alt_path,
                        org) {
  if (run_blst) {
    if (is.null(db_path)) {
      stop("Please provide the full path to the database")
    } else if (!dir.exists(db_path)) {
      stop(paste("The directory: ",
        db_path,
        " does not exist",
        sep = ""
      ))
    }

    if (is.null(db_name)) {
      stop("Please provide the full name of the database")
    }
  }

  if (!is.null(acsn_list)) {
    if (is.null(acsn_path)) {
      stop("Please provide a valid path to the accession ID list")
    } else if (!file.exists(
      here::here(
        acsn_path,
        acsn_list
      )
    )) {
      stop("Accession list does not exist in the specified path")
    }
  }

  if (run_blst &&
    file.exists(
      here::here(
        tab_path,
        paste("Alignment_",
          tab_out,
          ".tsv",
          sep = ""
        )
      )
    )) {
    owrt <- utils::askYesNo(
      "Blast annotations exist in the specified directory. Overwrite?"
    )
    if (owrt) {
      unlink(
        here::here(
          tab_path,
          paste("Alignment_",
            tab_out,
            ".tsv",
            sep = ""
          )
        )
      )
    }
  } else {
    owrt <- TRUE
  }

  base::assign(
    "owrt",
    owrt,
    .GlobalEnv
  )

  if (run_blst & owrt) {
    files_to_copy <- c(
      paste(db_path, "/taxdb.bti", sep = ""),
      paste(db_path, "/taxdb.btd", sep = "")
    )

    if (file.exists(files_to_copy[1]) &
      file.exists(files_to_copy[2])) {
      file.copy(files_to_copy,
        to = here::here(),
        overwrite = TRUE
      )
    } else {
      taxdb_check <- utils::askYesNo(
        paste("No taxdb files found in ",
          db_path,
          ". Abort?",
          sep = ""
        )
      )

      if (taxdb_check) {
        stop("Aborting: please specify correct directory with taxdb files")
      }
    }

    if (!is.null(acsn_list)) {
      if (acsn_check) {
        cmd <- "blastdb_aliastool"
        seq_in_list <- paste("-seqid_file_in ",
          here::here(
            acsn_path,
            acsn_list
          ),
          sep = ""
        )

        trmnl_cmd <- paste(cmd,
          seq_in_list,
          sep = " "
        )

        acc_check <- rstudioapi::terminalExecute(
          trmnl_cmd,
          show = FALSE
        )

        while (is.null(rstudioapi::terminalExitCode(acc_check))) {
          Sys.sleep(0.1)
        }
        if (rstudioapi::terminalExitCode(acc_check) == 0) {
          print("Accession list check successful")
          rstudioapi::terminalKill(acc_check)
          acsn_list <- paste(acsn_list,
            ".bsl",
            sep = ""
          )
        } else {
          stop(print(
            paste(
              "Accession list check ran into an error - Code: ",
              rstudioapi::terminalExitCode(acc_check),
              " - Please check terminal for details"
            )
          ))
        }
      } else {
        acsn_list <- paste(acsn_list,
          ".bsl",
          sep = ""
        )

        if (!file.exists(
          here::here(acsn_path, acsn_list)
        )) {
          stop(print(
            paste("accession list not found in specified directory")
          ))
        }
      }
    }

    task_list <- c(
      "megablast",
      "blastn",
      "blastn-short",
      "dc-megablast"
    )
    if (!task %in% task_list) {
      stop(paste("Please specify one of the following tasks: ",
        paste(task_list, collapse = "; "),
        sep = ""
      ))
    }

    if (!dir.exists(tab_path)) {
      dir.create(tab_path,
        recursive = TRUE
      )
    }
  }

  if (alt_annot) {
    list_dbs <- list.files(alt_path)

    if (org == "bac") {
      to_find <- list(
        "silva.*train",
        "silva.*species",
        "rdp.*train",
        "rdp.*species"
      ) %>%
        purrr::map(.f = function(x) {
          check <- stringr::str_subset(list_dbs, x) %>%
            purrr::is_empty()
          if (check) {
            stop(paste(x, "not found in the specified directory"))
          }
        })

      if (file.exists(
        here::here(
          tab_path,
          paste0(
            tab_out,
            "_silva_annot.fst"
          )
        )
      )) {
        owrt_silva <- utils::askYesNo(
          "Silva annotations exist in the specified directory. Overwrite?"
        )
      } else {
        owrt_silva <- TRUE
      }

      assign(
        "owrt_silva",
        owrt_silva,
        .GlobalEnv
      )

      if (file.exists(
        here::here(
          tab_path,
          paste0(
            tab_out,
            "_rdp_annot.fst"
          )
        )
      )) {
        owrt_RDP <- utils::askYesNo(
          "RDP annotations exist in the specified directory. Overwrite?"
        )
      } else {
        owrt_RDP <- TRUE
      }

      assign(
        "owrt_RDP",
        owrt_RDP,
        .GlobalEnv
      )
    }

    if (org == "fungi") {
      to_find <- list(
        "sh_general.*fasta"
      ) %>%
        purrr::map(.f = function(x) {
          check <- stringr::str_subset(list_dbs, x) %>%
            purrr::is_empty()
          if (check) {
            stop(paste(x, "not found in the specified directory"))
          }
        })

      if (file.exists(
        here::here(
          tab_path,
          paste0(
            tab_out,
            "_UNITE_annot.fst"
          )
        )
      )) {
        owrt_unite <- utils::askYesNo(
          "UNITE annotations exist in the specified directory. Overwrite?"
        )
      } else {
        owrt_unite <- TRUE
      }

      assign(
        "owrt_unite",
        owrt_unite,
        .GlobalEnv
      )
    }
  }
}

check_ecosrc <- function(hit_tbl,
                         filt_host,
                         filt_site,
                         filt_negt,
                         alt_tbl_path,
                         alt_tbl_name,
                         org,
                         add_scrs,
                         asbd_tbl,
                         asgn_tbl) {
  if (nrow(hit_tbl) == 0) {
    stop("Empty Data.frame - Aborting Relevance filtration")
  }

  filt_split <- unlist(
    stringr::str_split(c(
      filt_host,
      filt_site,
      filt_negt
    ),
    pattern = "\\+"
    )
  )

  filt_split <- filt_split[!is.na(filt_split)]

  wb_step <- 1
  while (wb_step <= length(filt_split)) {
    if (!filt_split[wb_step] %in% names(Word_banks)) {
      stop(paste("Invalid filtration entry ",
        filt_split[wb_step],
        sep = ""
      ))
    } else {
      wb_step <- wb_step + 1
    }
  }


  if (add_scrs) {
    list_tbls <- list.files(alt_tbl_path)

    if (org == "bac") {
      to_find <- list(
        paste0(alt_tbl_name, "_silva_annot.fst"),
        paste0(alt_tbl_name, "_rdp_annot.fst")
      ) %>%
        purrr::map(.f = function(x) {
          check <- stringr::str_subset(list_tbls, x) %>%
            purrr::is_empty()
          if (check) {
            stop(paste(x, "not found in the specified directory"))
          }
        })
    }

    if (org == "fungi") {
      to_find <- list(
        paste0(
          alt_tbl_name,
          "_UNITE_annot.fst"
        )
      ) %>%
        purrr::map(.f = function(x) {
          check <- stringr::str_subset(list_tbls, x) %>%
            purrr::is_empty()
          if (check) {
            stop(paste(x, "not found in the specified directory"))
          }
        })
    }
  }

  if (!is.na(asbd_tbl)) {
    if (!file.exists(asbd_tbl)) {
      stop(paste("Pre-assembled table does not exist in this location: ",
        asbd_tbl,
        sep = ""
      ))
    }
  }

  if (is.na(asbd_tbl)) {
    if (file.exists(asgn_tbl)
    ) {
      stop(paste0(
        "The following table already exists: ",
        asgn_tbl,
        ". Please assign it as 'asbd_tbl'",
        "or define a different name/directory to create a new table"
      ))
    }
  }
}

check_lineage <- function(taxids,
                          asbd_tbl,
                          asgn_tbl) {
  if (!any(names(taxids) %in% "TaxID")) {
    stop(paste("No column named 'TaxID' -",
      "please provide a valid input column",
      sep = ""
    ))
  }

  if (!is.na(asbd_tbl)) {
    if (!file.exists(asbd_tbl)) {
      stop(paste("Pre-assembled table does not exist in this location: ",
        asbd_tbl,
        sep = ""
      ))
    }
  }

  if (is.na(asbd_tbl)) {
    if (file.exists(asgn_tbl)
    ) {
      stop(paste0(
        "The following table already exists: ",
        asgn_tbl,
        ". Please assign it as 'asbd_tbl'",
        "or define a different name/directory to create a new table"
      ))
    }
  }
}


# txm_align miscellaneous functions ----
mk_splts <- function(obj,
                     wth_each,
                     tab) {
  out_path <- here::here(
    temp,
    tab
  )

  base::system2(
    command = "seqkit",
    args = paste(
      "split",
      "-p",
      wth_each,
      "-O",
      out_path,
      obj
    ),
    wait = TRUE
  )

  splits <- list.files(
    out_path,
    full.names = TRUE
  ) %>%
    as.list()
}

mk_fasta <- function(x) {
  x <- x %>%
    dplyr::mutate(ID = paste(">",
      ID,
      sep = ""
    )) %>%
    base::as.list() %>%
    purrr::pmap(~ paste(.x,
      .y,
      sep = "\n"
    )) %>%
    base::unlist()
}

splt_sqcs <- function(seq,
                      pth,
                      bch,
                      out,
                      index = TRUE) {
  temp <- here::here(
    pth,
    "temp"
  )

  assign("temp", temp, envir = .GlobalEnv)

  if (!dir.exists(temp)) {
    dir.create(temp,
      recursive = TRUE
    )
  }

  if (is.character(seq) && file.exists(seq)) {
    seq_obj <- seq
  } else {
    fa <- seq %>%
      base::colnames() %>%
      base::data.frame() %>%
      purrr::set_names("Seq") %>%
      dplyr::mutate(ID = paste0("seq_1.", 1:nrow(.))) %>%
      dplyr::select(
        ID,
        Seq
      ) %>%
      mk_fasta() %>%
      readr::write_lines(
        file = here::here(
          temp,
          "to_align.fa"
        )
      )

    seq_obj <- list.files(
      path = temp,
      pattern = "*.fa",
      full.names = TRUE
    )
  }

  if (index) {
    index <- system2(
      command = "seqkit",
      args = paste0(
        "fx2tab ",
        seq_obj,
        " -i > ",
        here::here(
          pth,
          paste0(
            "Index_",
            out,
            ".txt"
          )
        )
      )
    )
  }

  seq_splt <- mk_splts(
    seq_obj,
    bch,
    "batch"
  )
}


sys_cmd <- function(task,
                    db_path,
                    db_name,
                    tab_path,
                    threads,
                    pctidt,
                    acsn_path,
                    acsn_list,
                    max_out,
                    run_prl) {
  prl <- paste(
    "ls",
    here::here(
      tab_path,
      "*.fa"
    ),
    "| parallel -a -"
  )

  blst_cmd <- paste(
    "blastn -task ",
    task,
    sep = ""
  )

  db <- paste(
    "-db ",
    here::here(
      db_path,
      db_name
    ),
    sep = ""
  )

  query <- paste("-query {}")

  output <- paste(
    "-out ",
    here::here(
      tab_path,
      "Alignment_{/.}.tsv"
    ),
    sep = ""
  )

  parameters <- paste(
    "-num_threads",
    threads,
    "-perc_identity",
    pctidt,
    sep = " "
  )

  acsn_limit <- paste(
    "-seqidlist ",
    here::here(
      acsn_path,
      acsn_list
    ),
    ".bsl",
    sep = ""
  )

  out_fmt <- paste(
    "-max_hsps 1 -max_target_seqs",
    max_out,
    "-outfmt",
    "\"'6 qacc saccver staxids sscinames bitscore evalue qcovs pident'\""
  )

  cmd <- paste(
    prl,
    blst_cmd,
    db,
    query,
    output,
    if (!is.null(acsn_list)) {
      acsn_limit
    },
    parameters,
    out_fmt,
    sep = " "
  )

  base::system(cmd, wait = TRUE)
}


check_files <- function(x) {
  if (file.size(x) == 0) {
    stop(print("Empty query object provided"))
  }
}


algn_blast <- function(seq,
                       task,
                       db_path,
                       db_name,
                       tab_path,
                       tab_out,
                       threads,
                       pctidt,
                       acsn_path,
                       acsn_list,
                       max_out,
                       chunks) {
  mk_chunks <- seq %>%
    mk_splts(.,
      wth_each = chunks,
      "chunk"
    )

  chnk_path <- here::here(
    tab_path,
    "temp",
    "chunk"
  )

  prgrs_bar$tick(len = -1)
  prgrs_bar$tick()

  fsts <- list.files(
    path = chnk_path,
    pattern = "*.fa",
    full.names = T
  ) %>%
    as.list() %>%
    purrr::iwalk(~ check_files(.x))

  run <- sys_cmd(
    task = task,
    db_path = db_path,
    db_name = db_name,
    tab_path = chnk_path,
    threads = threads,
    pctidt = pctidt,
    acsn_path = acsn_path,
    acsn_list = acsn_list,
    max_out = max_out
  )

  prgrs_bar$tick()

  concat <- base::system(
    paste(
      "cat ",
      chnk_path,
      "/*.tsv >> ",
      here::here(
        tab_path,
        paste("Alignment_",
          tab_out,
          sep = ""
        )
      ),
      ".tsv",
      sep = ""
    )
  )

  unlink(chnk_path,
    recursive = TRUE
  )
}

alt_taxa <- function(seq,
                     db,
                     db_tr,
                     db_sp = NA,
                     pth,
                     org) {
  prgrs_bar$tick(len = -1)
  prgrs_bar$tick()

  annot <- dada2::assignTaxonomy(
    seq,
    here::here(
      pth,
      stringr::str_subset(
        db,
        db_tr
      )
    ),
    minBoot = 80,
    multithread = T,
    tryRC = T
  )

  if (!org == "fungi") {
    annot <- annot %>%
      dada2::addSpecies(
        here::here(
          pth,
          stringr::str_subset(
            db,
            db_sp
          )
        ),
        allowMultiple = T,
        tryRC = T
      ) %>%
      base::data.frame()
  } else {
    annot <- annot %>% 
      data.frame()
  }

  prgrs_bar$tick()

  annot
}



# Only needed until Fannyhessea vaginae is updated in silva and RDP databases
FHV <- data.frame(
  "Bacteria",
  "Actinobacteria",
  "Actinomycetia",
  "Actinomycetales",
  "Actinomycetaceae",
  "Fannyhessea",
  "vaginae"
) %>%
  purrr::set_names(c(
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
  ))

rplc <- function(lin) {
  lin_det <- as.character(
    paste(
      lin[6],
      lin[7]
    )
  )
  if (stringr::str_detect(
    lin_det,
    "Atopobium vaginae|Fannyhessea"
  )) {
    out <- FHV
  } else {
    out <- lin
  }
}


# txm_lineage miscellaneous functions ----
get_lge <- function(x) {
  chunk_lin <- try({
    TaxIDs_post <- rentrez::entrez_post(
      db = "taxonomy",
      id = purrr::flatten(x)$TaxID
    )

    TaxIDs_fetch <- rentrez::entrez_fetch(
      db = "taxonomy",
      web_history = TaxIDs_post,
      rettype = "xml",
      retmode = "xml",
      parsed = F
    )

    lge_xml <- xml2::read_xml(TaxIDs_fetch) %>%
      xml2::xml_children() %>%
      xml2::as_list()

    lge_mod1 <- lge_xml %>%
      purrr::modify_depth(
        .depth = 1,
        .f = ~ purrr::keep(
          .x,
          names(.x) %in% c(
            "TaxId",
            "ScientificName",
            "LineageEx"
          )
        )
      )

    lge_mod2 <- lge_mod1 %>%
      purrr::modify_depth(
        .depth = 1,
        .f = ~ purrr::modify_at(
          .x,
          .at = "LineageEx",
          .f = function(x) {
            x <- x %>%
              purrr::set_names(
                stringr::str_replace_all(
                  make.unique(
                    as.character(
                      purrr::map(
                        .,
                        .f = ~ unlist(.x[["Rank"]])
                      )
                    )
                  ),
                  " ", ""
                )
              ) %>%
              purrr::keep(!names(.) %in% "species") %>%
              purrr::map_df(
                .,
                .f = ~ unlist(.x[["ScientificName"]])
              )
          }
        )
      )

    lge_mod3 <- lge_mod2 %>%
      purrr::modify_depth(
        .depth = 1,
        .f = ~ purrr::modify_at(
          .x,
          .at = "TaxId",
          .f = ~ purrr::set_names(.x, "TaxID")
        )
      ) %>%
      purrr::modify_depth(
        .depth = 1,
        .f = ~ purrr::modify_at(
          .x,
          .at = "ScientificName",
          .f = ~ purrr::set_names(.x, "species")
        )
      ) %>%
      base::list()
  })

  if (!class(chunk_lin) == "try-error") {
    prgrs_bar$tick()
    Sys.sleep(1)
    chunk_lin
  } else {
    stop()
  }
}

get_lge <- purrr::insistently(
  get_lge,
  rate = purrr::rate_delay(
    10,
    max_times = 4
  ),
  quiet = F
)

lge_cln <- function(x) {
  x <- x %>%
    ifelse(
      stringr::str_detect(
        .,
        paste(
          "environmental samples",
          "incertae sedis",
          "unidentified",
          "clone",
          "\\/",
          "group",
          "enrichment",
          "culture",
          "endosymbiont",
          "\\bsp\\..",
          "\\bsp\\.$",
          "\\bbacterium\\b",
          "unclassified",
          "uncultured",
          sep = "|"
        )
      ),
      NA,
      .
    ) %>%
    stringr::str_remove_all(
      .,
      "\\s\\w{2}\\.|\\[|\\]|\\'"
    ) %>%
    stringr::str_extract_all(
      .,
      "^[^ ]* [^ ]*|^[^ ]*"
    ) %>%
    unlist()
}

# txm_ecosrc miscellaneous functions ----
annot_score <- function(x, y, z, org) {
  x_lin <- purrr::flatten(x) %>%
    purrr::keep(names(.) %in% c(
      "superkingdom",
      "phylum",
      "class",
      "order",
      "family",
      "genus",
      "species"
    ))

  if (org == "bac") {
    silva <- unlist(x_lin)[1:6] %>%
      stringr::str_count(unlist(y)[1:6]) %>%
      tidyr::replace_na(0) %>%
      base::sum()

    silva_sp <- unlist(x_lin)[7] %>%
      stringr::str_count(
        paste("\\b",
          base::unlist(
            stringr::str_split(
              y[[1]][7],
              pattern = "/"
            )
          ),
          "\\b",
          sep = ""
        )
      ) %>%
      base::sum() %>%
      tidyr::replace_na(0) * 2


    RDP <- unlist(x_lin)[1:6] %>%
      stringr::str_count(unlist(z)[1:6]) %>%
      tidyr::replace_na(0) %>%
      sum()

    RDP_sp <- unlist(x_lin)[7] %>%
      stringr::str_count(
        paste("\\b",
          base::unlist(
            stringr::str_split(
              y[[1]][7],
              pattern = "/"
            )
          ),
          "\\b",
          sep = ""
        )
      ) %>%
      base::sum() %>%
      tidyr::replace_na(0) * 2

    meta_scr <- purrr::flatten(x) %>%
      purrr::keep(
        stringr::str_detect(
          names(.),
          "PMID|host|Isolation_source"
        )
      ) %>%
      purrr::map(.f = ~ ifelse(!is.na(.), 1, 0)) %>%
      unlist() %>%
      sum() * 3

    score <- silva + silva_sp + RDP + RDP_sp + meta_scr
    df <- data.frame(
      silva,
      silva_sp,
      RDP,
      RDP_sp,
      meta_scr,
      score
    )
  }

  if (org == "fungi") {
    unite <- unlist(x_lin)[1:6] %>%
      stringr::str_count(unlist(y)[1:6]) %>%
      tidyr::replace_na(0) %>%
      base::sum()

    unite_sp <- unlist(x_lin)[7] %>%
      stringr::str_count(
        paste("\\b",
          base::unlist(
            stringr::str_split(
              y[[1]][7],
              pattern = "/"
            )
          ),
          "\\b",
          sep = ""
        )
      ) %>%
      base::sum() %>%
      tidyr::replace_na(0) * 2

    meta_scr <- purrr::flatten(x) %>%
      purrr::keep(
        stringr::str_detect(
          names(.),
          "PMID|host|Isolation_source"
        )
      ) %>%
      purrr::map(.f = ~ ifelse(!is.na(.), 1, 0)) %>%
      unlist() %>%
      sum() * 3

    score <- unite + unite_sp + meta_scr

    df <- data.frame(
      unite,
      unite_sp,
      meta_scr,
      score
    )
  }

  df
}

new_bar <- function(len,
                    msg = "downloading") {
  progress::progress_bar$new(
    format = paste(msg, " (:spin) [:bar]",
      ":current/:total",
      "(:percent) eta:",
      ":eta elapsed:",
      ":elapsedfull",
      sep = " "
    ),
    total = len,
    clear = FALSE,
    width = 90,
    show_after = 0
  )
}

save_temp <- function(x) {
  chunk <- x
  if (file.exists(here::here(
    "temp_files",
    "pmids_temp.rds"
  ))) {
    to_sv <- readRDS(here::here(
      "temp_files",
      "pmids_temp.rds"
    ))

    add <- append(to_sv, list(chunk))
    saveRDS(add, here::here(
      "temp_files",
      "pmids_temp.rds"
    ))
  } else {
    nwsv <- list()
    nwsv <- append(nwsv, list(chunk))
    saveRDS(nwsv, here::here(
      "temp_files",
      "pmids_temp.rds"
    ))
  }
}

rm_prv <- function(x) {
  if (file.exists(here::here(
    "temp_files",
    "pmids_temp.rds"
  ))) {
    to_rm <- readRDS(here::here(
      "temp_files",
      "pmids_temp.rds"
    ))
    assign("to_rm", to_rm, envir = .GlobalEnv)
    pmids <- x[-c(1:length(to_rm))]
  } else {
    assign("to_rm", list(), envir = .GlobalEnv)
    pmids <- x
  }
}

get_esumm <- function(accID,
                      sys.break = 1) {
  chunk_ID <- try({
    ids_post <- rentrez::entrez_post(
      db = "nuccore",
      id = paste(unlist(accID),
        collapse = ","
      )
    )

    ids_summ <- rentrez::entrez_summary(
      db = "nucleotide",
      web_history = ids_post,
      version = "2.0",
      retmode = "xml"
    )

    to_sub <- c(
      "AccessionVersion",
      "Strain",
      "Title",
      "TaxId",
      "SubType",
      "SubName"
    )
    
    ids_extr <- t(rentrez::extract_from_esummary(
      ids_summ,
      elements = c(to_sub)
    )) %>%
      base::as.data.frame() %>%
      purrr::set_names(c(to_sub)) %>% 
      dplyr::mutate(
        dplyr::across(
          tidyselect::vars_select_helpers$where(is.list),
          .fns = ~ as.character(.x)
        )
      )
  })

  if (!class(chunk_ID) == "try-error") {
    prgrs_bar$tick()
    Sys.sleep(sys.break)
  } else {
    stop()
  }
  return(chunk_ID)
}


get_esumm <- purrr::insistently(
  get_esumm,
  rate = purrr::rate_delay(10,
    max_times = 4
  ),
  quiet = F
)

get_pbids <- function(accID,
                      sys.break = 1) {
  chunk_link <- try({
    ids_elink <- rentrez::entrez_link(
      db = "pubmed",
      dbfrom = "nuccore",
      id = purrr::flatten(accID)$AccID,
      cmd = "neighbor",
      by_id = TRUE,
      rettype = "native",
      idtype = "acc"
    )

    if (length(purrr::flatten(accID)$AccID) == 1) {
      ids_elink <- list(ids_elink)
    }

    ids_elink <- ids_elink %>%
      purrr::set_names(
        purrr::flatten(accID)$AccID
      ) %>%
      purrr::map(
        .f = function(x) {
          x$links[["nuccore_pubmed"]]
        }
      ) %>%
      purrr::compact()

    if (length(ids_elink) > 0) {
      ids_elink <- ids_elink %>%
        purrr::map(.f = ~ purrr::pluck(.x, 1)) %>%
        dplyr::bind_rows() %>%
        tidyr::pivot_longer(
          cols = tidyselect::everything(),
          names_to = "AccID",
          values_to = "PMID"
        ) %>%
        base::list()
    } else {
      ids_elink <- list()
    }
  })

  if (!class(chunk_link) == "try-error") {
    prgrs_bar$tick()
    purrr::iwalk(chunk_link, ~ save_temp(.x))
    Sys.sleep(sys.break)
    chunk_link
  } else {
    stop()
  }
}


get_pbids <- purrr::insistently(
  get_pbids,
  rate = purrr::rate_delay(10,
    max_times = 4
  ),
  quiet = F
)

get_pbdt <- function(pmid,
                     sys.break = 1) {
  chunk_pub <- try({
    pmids_post <- rentrez::entrez_post(
      db = "pubmed",
      id = purrr::flatten(pmid)$PMID
    )

    pmids_fetch <- rentrez::entrez_fetch(
      db = "pubmed",
      web_history = pmids_post,
      rettype = "xml",
      retmode = "xml",
      parsed = F
    )

    pmids_cn1 <- xml2::read_xml(pmids_fetch) %>%
      xml2::xml_children() %>%
      xml2::as_list() %>%
      purrr::modify_depth(
        .depth = 1,
        .f = ~ purrr::keep(
          .x,
          names(.x) %in% "MedlineCitation"
        )
      ) %>%
      purrr::modify_depth(
        .depth = 2,
        .f = ~ modify_at(
          .x,
          .at = "PMID",
          .f = ~ purrr::set_names(.x, "PMID")
        )
      )

    pmids_cn2 <- pmids_cn1 %>%
      purrr::modify_depth(
        .depth = 2,
        .f = ~ purrr::keep(
          .x,
          names(.) %in% c(
            "PMID",
            "Article",
            "MeshHeadingList"
          )
        )
      )


    pmids_cn3 <- pmids_cn2 %>%
      purrr::modify_depth(
        .depth = 2,
        .f = ~ purrr::modify_at(
          .x,
          .at = "Article",
          .f = ~ purrr::keep(
            .x,
            names(.x) %in% c(
              "Journal",
              "ArticleTitle",
              "Abstract"
            )
          )
        )
      )

    pmids_cn4 <- pmids_cn3 %>%
      purrr::modify_depth(
        .depth = 2,
        .f = ~ purrr::modify_at(
          .x,
          .at = "Article",
          .f = ~ purrr::modify_at(
            .x,
            .at = "Journal",
            .f = ~ purrr::pluck(
              .x,
              "Title"
            )
          )
        )
      )

    pmids_cn5 <- pmids_cn4 %>%
      purrr::modify_depth(
        .depth = 2,
        .f = ~ purrr::modify_at(
          .x,
          .at = "MeshHeadingList",
          .f = function(x) {
            x <- x %>%
              purrr::flatten() %>%
              purrr::reduce(paste)
          }
        )
      )

    pmids_cn6 <- pmids_cn5 %>%
      purrr::modify_depth(
        .depth = 2,
        .f = ~ purrr::modify_at(
          .x,
          .at = "Article",
          .f = ~ purrr::modify_depth(
            .x,
            .depth = 1,
            .f = function(x) {
              x <- x %>%
                purrr::flatten() %>%
                purrr::reduce(paste)
            }
          )
        )
      ) %>%
      purrr::map_dfr(.f = dplyr::bind_cols)

    if (!any(names(pmids_cn6) == "MeshHeadingList")) {
      pmids_cn6 <- pmids_cn6 %>%
        tibble::add_column("Mesh" = NA)
    }

    pmids_cn7 <- pmids_cn6 %>%
      purrr::set_names(
        c(
          "PMID",
          "Journal",
          "Title",
          "Abstract",
          "Mesh"
        )
      ) %>%
      dplyr::mutate(Article = paste(
        "Journal: ",
        Journal,
        " - Title: ",
        Title,
        " - Abstract: ",
        Abstract,
        sep = ""
      )) %>%
      dplyr::select(-c(
        Journal,
        Title,
        Abstract
      )) %>%
      base::list()
  })

  if (!class(chunk_pub) == "try-error") {
    prgrs_bar$tick()
    Sys.sleep(sys.break)
    chunk_pub
  } else {
    stop()
  }
}

get_pbdt <- purrr::insistently(
  get_pbdt,
  rate = purrr::rate_delay(10,
    max_times = 4
  ),
  quiet = F
)

host_kp <- function(x, y) {
  host_filt <- x %>%
    data.frame() %>%
    dplyr::mutate(
      dplyr::across(
        .cols = tidyselect::vars_select_helpers$where(is.factor),
        .fns = ~ as.character(.x)
      )
    ) %>%
    dplyr::filter(
      stringr::str_detect(
        pooled_data,
        stringr::regex(
          as.character(
            unlist(y)
          ),
          ignore_case = T
        )
      )
    ) %>%
    dplyr::filter(
      is.na(host) | stringr::str_detect(
        host,
        pattern = stringr::regex(
          as.character(
            unlist(y)
          ),
          ignore_case = T
        )
      )
    )
}

negt <- function(x, y) {
  negt_filt <- x %>%
    dplyr::mutate(
      dplyr::across(
        .cols = tidyselect::vars_select_helpers$where(is.factor),
        .fns = ~ as.character(.x)
      )
    ) %>%
    dplyr::filter(
      stringr::str_detect(
        pooled_data,
        stringr::regex(
          as.character(
            unlist(y)
          ),
          ignore_case = T
        ),
        negate = T
      ) |
        stringr::str_detect(
          pooled_data,
          "ISHAM"
        ) |
        stringr::str_detect(
          host,
          as.character(
            unlist(y)
          ),
          negate = T
        )
    )
}

site_kp <- function(x, y) {
  site_filt <- x %>%
    dplyr::mutate(
      dplyr::across(
        .cols = tidyselect::vars_select_helpers$where(is.factor),
        .fns = ~ as.character(.x)
      )
    ) %>%
    dplyr::filter(
      stringr::str_detect(
        pooled_data,
        stringr::regex(
          as.character(
            unlist(y)
          ),
          ignore_case = T
        )
      )
    ) %>%
    dplyr::filter(
      is.na(Isolation_source) |
        stringr::str_detect(
          Isolation_source,
          pattern = stringr::regex(
            as.character(
              unlist(y)
            ),
            ignore_case = T
          )
        )
    )
}

concat <- function(x) {
  clean <- x %>%
    dplyr::select(
      SubType,
      SubName
    ) %>%
    as.list() %>%
    purrr::map(
      .f = ~ stringr::str_replace_all(
        .x,
        pattern = "\\|",
        replacement = ";"
      )
    ) %>%
    purrr::map(
      .f = ~ stringr::str_split(
        .x,
        pattern = ";"
      )
    ) %>%
    purrr::pmap(~ paste(.x,
      .y,
      sep = ": "
    )) %>%
    purrr::map(.f = ~ paste(.x,
      collapse = " - "
    )) %>%
    unlist()
}

clpse_wbks <- function(x) {
  if (stringr::str_detect(x, "\\+")) {
    site_split <- unlist(
      stringr::str_split(
        x,
        pattern = "\\+"
      )
    )
    site_cmp <- Word_banks %>%
      purrr::keep(names(Word_banks) %in% site_split) %>%
      purrr::flatten() %>%
      purrr::map_chr(
        .data,
        .f = ~ paste("\\b",
          .x,
          "\\b",
          sep = ""
        )
      ) %>%
      paste(collapse = "|") %>%
      as.list()
  } else {
    as.list(purrr::pluck(Word_banks, x)) %>%
      purrr::map_chr(
        .data,
        .f = ~ paste("\\b",
          .x, "\\b",
          sep = ""
        )
      ) %>%
      paste(collapse = "|") %>%
      as.list()
  }
}

clpse_lists <- function(x, y) {
  z <- x %>%
    dplyr::bind_rows(
      y <- y %>%
        dplyr::filter(!ID %in% x$ID)
    )
}

hm_wbk <- purrr::map(
  "human",
  .f = clpse_wbks
) %>%
  purrr::set_names("human")

meta_extr <- function(x) {
  
  fct2chr <- function(x) {
    if (is.factor(x)) {
      x <- as.character(x)
    } else {
      x
    }
  }
  
  df_extrt <- x %>%
    dtplyr::lazy_dt() %>%
    dplyr::mutate(
      dplyr::across(
        .cols = everything(),
        ## "where" function breaks inside dtplyr. 
        .fns = ~ fct(.x)
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(host = trimws(
      stringr::str_extract(
        .data$meta,
        pattern = stringr::regex("(?<=host:).*?(?=\\s-)|(?<= - host: ).*")
      )
    )) %>%
    dplyr::mutate(Isolation_source = trimws(
      stringr::str_extract(
        .data$meta,
        pattern = stringr::regex(
          paste("(?<=isolation_source:).*?(?=\\s-)",
            "(?<=isolation_source:).*?(?=$)",
            sep = "|"
          )
        )
      )
    )) %>%
    dplyr::mutate(
      host = dplyr::case_when(
        (is.na(host) & stringr::str_detect(
          .data$Isolation_source,
          pattern = unlist(hm_wbk)
        )) ~ "Homo sapiens",
        TRUE ~ host
      )
    ) %>%
    as.data.frame()
}

save_dscr <- function(x, y, z) {
  t <- y %>%
    dplyr::filter(!.$AccID %in% x$AccID) %>%
    fst::write_fst(
      path = z,
      compress = 100
    )

  print(paste0(
    nrow(t),
    " discarded hits written to file: ",
    z
  ))
}
