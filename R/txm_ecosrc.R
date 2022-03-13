utils::globalVariables(c(
  "AccessionVersion",
  "TaxId",
  "AccID",
  "Article",
  "pooled_data",
  "host",
  "Organism",
  "Isolation_source",
  "ID",
  ".",
  "Journal",
  "Mesh",
  "data",
  "silva",
  "RDP",
  "score",
  "to_rm",
  "ArticleTitle",
  "Abstract"
))

#' Text mining and filtration
#'
#' Extract information from NCBI nucleotide and PubMed databases, attaching 
#' ecosystem specificity to each accession ID. Different combinations of word 
#' banks are used to scan through this data and apply the filtration 
#' criteria. The accession IDs are split into batches of 200, 
#' and \link[rentrez]{rentrez} is used to communicate with NCBI.
#' An annotation score is calculated using this information and the consensus
#' between the annotations methods/databases. 
#'
#' @param hit_tbl (Required) Default NULL. Output table obtained 
#' from \link[taxminer]{txm_align}. 
#' @param alt_tbl_path (Required) Default NULL. Full path to silva and RDP
#' annotations.
#' @param filt_host (Optional) Default NA. Filter annotations by host
#' @param filt_site (Optional) Default NA. Filter annotations by body site or 
#' environment.
#' @param filt_negt (Optional) Default NA. Disregard annotations that contain 
#' these terms.
#' @param do_filt (Logical) Default TRUE. Perform filtration using the word 
#' banks. If FALSE the output will contain all accession IDs and the extracted 
#' information associated to them.
#' @param asbd_tbl (Optional) Default NA. Specify the name of the pre-compiled 
#' database present within the directory.
#' @param asgn_tbl (Optional) Default Dataset + system date. Name of a new 
#' compiled database to be assigned
#' @param sys.break (integer) Default 1. Amount of time, in seconds, that the
#'                  system is paused between iterations. This is handy with larger
#'                  queries to reduce to load on the NCBI servers.
#'
#' @export
txm_ecosrc <- function(hit_tbl,
                       filt_host = NA,
                       filt_site = NA,
                       filt_negt = NA,
                       alt_tbl_path = NULL,
                       do_filt = T,
                       asbd_tbl = NA,
                       asgn_tbl = paste("Dataset_",
                         Sys.Date(),
                         ".fst",
                         sep = ""
                       ),
                       sys.break = 1) {
  check_ecosrc(
    hit_tbl,
    filt_host,
    filt_site,
    filt_negt,
    alt_tbl_path
  )

  ##### Input preparation -----
  ids_src <- hit_tbl %>%
    dplyr::distinct(.data$AccID) %>%
    dplyr::filter(stringr::str_detect(.data$AccID, "\\."))

  if (!is.na(asbd_tbl)) {
    print("Reading in dataset and searching for existing Accession ids")

    asbd_tbl_sub <- fst::read_fst(asbd_tbl) %>%
      dtplyr::lazy_dt() %>%
      dplyr::inner_join(ids_src, by = "AccID") %>%
      base::as.data.frame()

    ids_src <- ids_src %>%
      dplyr::filter(!.data$AccID %in% asbd_tbl_sub$AccID) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(.data$AccID)

    print(paste(nrow(ids_src),
      "of",
      sum(
        nrow(asbd_tbl_sub),
        nrow(ids_src)
      ),
      "will be searched for",
      sep = " "
    ))
  }

  if (nrow(ids_src) > 0) {
    if (is.na(asbd_tbl)) {
      print(paste("No dataset provided -", nrow(ids_src),
        "accession ids will be searched for",
        sep = " "
      ))
    }

    splits <- round(nrow(ids_src) / 200)
    if (splits < 1) splits <- 1

    ids_src <- ids_src %>%
      dplyr::mutate(chunks = rep(0:splits,
        each = 200,
        length.out = nrow(.)
      )) %>%
      dplyr::group_by(.data$chunks) %>%
      tidyr::nest() %>%
      dplyr::pull()

    if (!dir.exists(
      here::here("temp_files")
    )) {
      dir.create(
        here::here("temp_files"),
        recursive = T
      )
    }

    ##### Accession ID retrieval -----
    if (!file.exists(
      here::here(
        "temp_files",
        "AccID_temp.fst"
      )
    )) {
      prgrs_bar <- new_bar(length(ids_src))
      base::assign(
        "prgrs_bar",
        prgrs_bar,
        .GlobalEnv
      )

      print("Retrieving Accession ID data")

      df_summ <- ids_src %>%
        purrr::map_dfr(~ get_esumm(.x, 
                                   sys.break)) %>%
        dplyr::mutate(meta = concat(.)) %>%
        dplyr::select(-c(
          .data$SubType,
          .data$SubName
        )) %>%
        dplyr::mutate(chunks = rep(0:nrow(.),
          each = 100000,
          length.out = nrow(.)
        )) %>%
        dplyr::group_by(.data$chunks) %>%
        tidyr::nest() %>%
        dplyr::pull() %>%
        purrr::map_dfr(.f = ~ meta_extr(.x))

      fst::write_fst(
        df_summ,
        here::here(
          "temp_files",
          "AccID_temp.fst"
        ),
        100
      )
    } else {
      df_summ <- fst::read_fst(
        here::here(
          "temp_files",
          "AccID_temp.fst"
        )
      )
    }


    ##### PMIDs retrieval -----
    pmids <- ids_src %>%
      rm_prv()

    prgrs_bar <- new_bar(length(pmids))
    base::assign(
      "prgrs_bar",
      prgrs_bar,
      .GlobalEnv
    )

    print("Retrieving PMIDs")

    pmids <- pmids %>%
      purrr::lmap(~ get_pbids(.x, 
                              sys.break)) %>%
      purrr::prepend(to_rm) %>%
      dplyr::bind_rows() %>%
      dplyr::group_by(.data$AccID) %>%
      dplyr::distinct(.data$PMID,
        .keep_all = T
      ) %>%
      dplyr::ungroup() %>%
      base::data.frame()

    rm(to_rm, envir = .GlobalEnv)

    ##### PubMed data retrieval -----
    if (nrow(pmids) > 0) {
      pmids_src <- pmids %>%
        dplyr::select(.data$PMID) %>%
        tidyr::drop_na() %>%
        dplyr::distinct()

      splits <- round(nrow(pmids_src) / 200)
      if (splits < 1) splits <- 1

      pmids_src <- pmids_src %>%
        dplyr::mutate(chunks = rep(0:splits,
          each = 200,
          length.out = nrow(.)
        )) %>%
        dplyr::group_by(.data$chunks) %>%
        tidyr::nest() %>%
        dplyr::pull()

      prgrs_bar <- new_bar(length(pmids_src))
      base::assign(
        "prgrs_bar",
        prgrs_bar,
        .GlobalEnv
      )

      print("Retrieving PubMed data")

      pbdt <- pmids_src %>%
        purrr::lmap(~ get_pbdt(.x, 
                               sys.break)) %>%
        purrr::map_dfr(.f = dplyr::bind_rows) %>%
        dplyr::right_join(pmids,
          by = "PMID"
        )

      df_summ <- df_summ %>%
        dplyr::rename("AccID" = AccessionVersion) %>%
        dplyr::left_join(pbdt,
          by = "AccID"
        )
    } else {
      df_summ <- df_summ %>%
        dplyr::rename("AccID" = AccessionVersion) %>%
        dplyr::mutate(Mesh = NA) %>%
        dplyr::mutate(Article = NA)
    }

    if (!is.null(df_summ)) {
      if (!is.na(asbd_tbl)) {
        if (!is.null(ids_src)) {
          print("Writing new data to Pre-compiled dataset")
          asbd_full <- fst::read_fst(asbd_tbl) %>%
            dplyr::bind_rows(df_summ) %>%
            dplyr::distinct(.data$AccID,
              .keep_all = T
            ) %>%
            fst::write_fst(
              asbd_tbl,
              100
            )
        }
      } else {
        print("Creating dataset")
        fst::write_fst(
          df_summ,
          asgn_tbl,
          100
        )
      }
    }

    if (!is.na(asbd_tbl)) {
      df_summ <- df_summ %>%
        dplyr::bind_rows(asbd_tbl_sub)
    }
  } else {
    df_summ <- asbd_tbl_sub
  }

  unlink("temp_files",
    recursive = T
  )

  df_summ <- df_summ %>%
    {
      if ("TaxID" %in% names(hit_tbl)) {
        dplyr::select(., -.data$TaxId)
      } else {
        dplyr::rename(., "TaxID" = .data$TaxId)
      }
    } %>%
    dplyr::inner_join(hit_tbl) %>%
    dplyr::arrange(.data$ID) %>%
    dplyr::select(
      .data$ID,
      .data$AccID,
      tidyselect::everything()
    )


  ##### filter preparation ----
  if (do_filt) {
    if (!is.na(filt_host[1])) {
      host_wbk <- purrr::map(
        filt_host,
        .f = clpse_wbks
      ) %>%
        purrr::set_names(filt_host)
    } else {
      host_wbk <- NA
    }

    if (!is.na(filt_site[1])) {
      site_wbk <- purrr::map(
        filt_site,
        .f = clpse_wbks
      ) %>%
        purrr::set_names(filt_site)
    } else {
      site_wbk <- NA
    }

    if (!is.na(filt_negt[1])) {
      negt_wbk <- purrr::map(
        filt_negt,
        .f = clpse_wbks
      ) %>%
        purrr::set_names(filt_negt)
    } else {
      negt_wbk <- NA
    }

    df_filt <- df_summ %>%
      tidyr::unite("pooled_data",
        c(
          .data$Mesh,
          .data$Article,
          .data$meta
        ),
        na.rm = T
      )

    ##### filtration -----
    if (!is.na(negt_wbk[1])) {
      print(paste(
        "Removing",
        paste(names(negt_wbk),
          collapse = " & "
        ),
        "sites"
      ))

      df_filt <- df_filt %>%
        base::list() %>%
        purrr::map2(
          .,
          negt_wbk,
          negt
        ) %>%
        purrr::reduce(.f = clpse_lists) %>%
        dplyr::group_by(ID) %>%
        dplyr::distinct(AccID,
          .keep_all = T
        ) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(ID) %>%
        purrr::walk(save_dscr(
          .,
          df_filt,
          here::here("negt_discard.xlsx")
        ))
    }

    if (!is.na(host_wbk[1])) {
      print(paste(
        "Searching through",
        paste(names(host_wbk),
          collapse = " & "
        ),
        "site banks"
      ))

      df_filt <- df_filt %>%
        base::list() %>%
        purrr::map2(
          .,
          host_wbk,
          host_kp
        ) %>%
        purrr::reduce(.f = clpse_lists) %>%
        dplyr::group_by(ID) %>%
        dplyr::distinct(AccID,
          .keep_all = T
        ) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(ID) %>%
        purrr::walk(save_dscr(
          .,
          df_filt,
          here::here("host_discard.xlsx")
        ))
    }

    if (!is.na(site_wbk[1])) {
      print(paste(
        "Searching through",
        paste(names(site_wbk),
          collapse = " & "
        ),
        "site banks"
      ))

      df_filt <- df_filt %>%
        base::list() %>%
        purrr::map2(
          .,
          site_wbk,
          .f = site_kp
        ) %>%
        purrr::reduce(.f = clpse_lists) %>%
        dplyr::group_by(ID) %>%
        dplyr::distinct(AccID,
          .keep_all = T
        ) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(ID) %>%
        purrr::walk(save_dscr(
          .,
          df_filt,
          here::here("site_discard.xlsx")
        ))
    }
    df_out <- df_summ %>%
      dplyr::filter(AccID %in% df_filt$AccID)
  } else {
    df_out <- df_summ
  }


  ##### add scores for annotations -----
  silva <- fst::read_fst(
    here::here(
      alt_tbl_path,
      "Silva_annot.fst"
    )
  ) %>%
    dplyr::group_by(ASVs) %>%
    tidyr::nest(silva = !ASVs)

  rdp <- fst::read_fst(
    here::here(
      alt_tbl_path,
      "RDP_annot.fst"
    )
  ) %>%
    dplyr::group_by(ASVs) %>%
    tidyr::nest(RDP = !ASVs)

  prgrs_bar <- new_bar(
    nrow(df_out),
    "adding score"
  )
  assign(
    "prgrs_bar",
    prgrs_bar,
    .GlobalEnv
  )

  blst <- df_out %>%
    dplyr::group_by(
      .data$ID,
      .data$ASVs,
      .data$AccID
    ) %>%
    dplyr::distinct(.data$AccID,
      .keep_all = TRUE
    ) %>%
    tidyr::nest() %>%
    dplyr::left_join(silva,
      by = "ASVs"
    ) %>%
    dplyr::left_join(rdp,
      by = "ASVs"
    ) %>%
    dplyr::mutate(
      score = purrr::pmap(
        .l = list(
          data,
          silva,
          RDP
        ),
        .f = ~ annot_score(
          data,
          silva,
          RDP
        )
      )
    ) %>%
    dplyr::select(-c(silva, RDP)) %>%
    tidyr::unnest(cols = -c(
      ID,
      ASVs,
      AccID
    )) %>%
    dplyr::arrange(
      ID,
      dplyr::desc(score)
    )
  
  rm(prgrs_bar, envir = .GlobalEnv)
  
  blst
}
