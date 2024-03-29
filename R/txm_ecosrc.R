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
#' banks are used to scan through the data and apply the filtration
#' criteria. The accession IDs are split into batches of 200,
#' and \link[rentrez]{rentrez} is used to communicate with NCBI.
#' For each annotation hit a "score" is calculated based on the consensus
#' between the superkingdom-species level lineage obtained from
#' BLAST + Silva + RDP and the completeness of the available
#' ecosystem/publication information.
#' Uses a NCBI api key is highly recommended.
#' \href{  https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html#rate-limiting-and-api-keys}{rentrez vignette}
#'
#' @param hit_tbl (data.frame) **Required** Default NULL. Output table obtained
#' from \link[taxminer]{txm_align}.
#' @param alt_tbl_path (String) **Required** Default NULL. Full path to silva
#' and RDP annotation outputs.
#' @param org (String) Default "bac". The organism being investigated.
#' \itemize{
#'   \item "bac" - Bacteria. Silva and RDP databases are used for alternative
#'   annotations to BLAST.
#'   \item "fungi" - Fungi. UNITE database is used for alternative annotations
#'   to BLAST.
#' }
#' @param filt_host (String) Default NA. Filter annotations by host
#' @param filt_site (String) Default NA. Filter annotations by body site or
#' environment.
#' @param filt_negt (String) Default NA. Disregard annotations that contain
#' these terms.
#' @param do_filt (Logical) Default FALSE. Perform filtration using the word
#' banks. If FALSE the output will contain all accession IDs and the extracted
#' information associated to them.
#' @param add_scrs (Logical) Default FALSE. Add scores to each hit based on
#' alternative alignments to silva, RDP, or UNITE databases.
#' @param asbd_tbl (String) Default NA. Specify the name of the pre-compiled
#' database present within the directory.
#' @param asgn_tbl (String) Default Dataset + system date. Name of a new
#' compiled database to be assigned
#' @param sys.break (>0) Default 1. Amount of time, in seconds, that the
#'                  system is paused between iterations. This is handy with larger
#'                  queries to reduce the load on the NCBI servers.
#'
#' @export
txm_ecosrc <- function(hit_tbl,
                       filt_host = NA,
                       filt_site = NA,
                       filt_negt = NA,
                       alt_tbl_path = NULL,
                       alt_tbl_name = NULL,
                       org = "bac",
                       do_filt = FALSE,
                       add_scrs = FALSE,
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
    alt_tbl_path,
    alt_tbl_name,
    org,
    add_scrs,
    asbd_tbl,
    asgn_tbl
  )

  # Input preparation ----
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

    # Accession ID retrieval ----
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
        purrr::map_dfr(~ get_esumm(
          .x,
          sys.break
        )) %>%
        dplyr::mutate(meta = concat(.)) %>%
        dplyr::select(-c(
          SubType,
          SubName
        )) %>%
        dplyr::mutate(chunks = rep(0:nrow(.),
          each = 100000,
          length.out = nrow(.)
        )) %>%
        dplyr::group_by(.data$chunks) %>%
        tidyr::nest() %>%
        dplyr::pull() %>%
        purrr::map_dfr(.f = ~ meta_extr(.x))

      rm(prgrs_bar, envir = .GlobalEnv)

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


    # PMIDs retrieval ----
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
      purrr::lmap(~ get_pbids(
        .x,
        sys.break
      )) %>%
      base::append(to_rm, after = 0) %>%
      dplyr::bind_rows()

    if (nrow(pmids) > 0) {
      pmids <- pmids %>%
        dplyr::group_by(.data$AccID) %>%
        dplyr::distinct(.data$PMID,
          .keep_all = T
        ) %>%
        dplyr::ungroup() %>%
        base::data.frame()
    }

    rm(to_rm, prgrs_bar, envir = .GlobalEnv)

    # PubMed data retrieval ----
    if (nrow(pmids) > 0) {
      pmids_src <- pmids %>%
        dplyr::select(PMID) %>%
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
        purrr::lmap(~ get_pbdt(
          .x,
          sys.break
        )) %>%
        purrr::map_dfr(.f = dplyr::bind_rows) %>%
        dplyr::right_join(pmids,
          by = "PMID"
        )

      df_summ <- df_summ %>%
        dplyr::rename("AccID" = AccessionVersion) %>%
        dplyr::left_join(pbdt,
          by = "AccID"
        )

      rm(prgrs_bar, envir = .GlobalEnv)
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
        dplyr::select(., -TaxId)
      } else {
        dplyr::rename(., "TaxID" = .data$TaxId)
      }
    } %>%
    dplyr::inner_join(hit_tbl) %>%
    dplyr::arrange(.data$ID) %>%
    dplyr::select(
      ID,
      AccID,
      tidyselect::everything()
    )


  # filter preparation ----
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
          Mesh,
          Article,
          meta
        ),
        na.rm = T
      )

    # filtration ----
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
          here::here("negt_discard.fst")
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
          here::here("host_discard.fst")
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
          here::here("site_discard.fst")
        ))
    }

    print("Creating final output")
    df_out <- df_summ %>%
      dtplyr::lazy_dt() %>%
      dplyr::group_by(ID, AccID) %>%
      dplyr::inner_join(df_filt) %>%
      dplyr::select(-pooled_data) %>%
      dplyr::ungroup() %>%
      data.frame()
  } else {
    print("Creating final output")
    df_out <- df_summ
  }

  if (add_scrs) {
    # add scores for annotations ----
    if (org == "bac") {
      silva <- fst::read_fst(
        here::here(
          alt_tbl_path,
          paste0(
            alt_tbl_name,
            "_silva_annot.fst"
          )
        )
      ) %>%
        dplyr::group_by(Seq) %>%
        tidyr::nest(silva = !Seq)

      rdp <- fst::read_fst(
        here::here(
          alt_tbl_path,
          paste0(
            alt_tbl_name,
            "_rdp_annot.fst"
          )
        )
      ) %>%
        dplyr::group_by(Seq) %>%
        tidyr::nest(RDP = !Seq)

      print("Adding scores...")
      blst <- df_out %>%
        dplyr::group_by(
          .data$ID,
          .data$Seq,
          .data$AccID
        ) %>%
        dplyr::distinct(.data$AccID,
          .keep_all = TRUE
        ) %>%
        tidyr::nest() %>%
        dplyr::left_join(silva,
          by = "Seq"
        ) %>%
        dplyr::left_join(rdp,
          by = "Seq"
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
              RDP,
              org = org
            )
          )
        ) %>%
        dplyr::select(-c(silva, RDP)) %>%
        tidyr::unnest(cols = -c(
          ID,
          Seq,
          AccID
        )) %>%
        dplyr::arrange(
          ID,
          dplyr::desc(score)
        )
    }

    if (org == "fungi") {
      unite <- fst::read_fst(
        here::here(
          alt_tbl_path,
          paste0(
            alt_tbl_name,
            "_UNITE_annot.fst"
          )
        )
      ) %>%
        dplyr::group_by(Seq) %>%
        tidyr::nest(unite = !Seq)

      print("Adding scores...")
      blst <- df_out %>%
        dplyr::group_by(
          .data$ID,
          .data$Seq,
          .data$AccID
        ) %>%
        dplyr::distinct(.data$AccID,
          .keep_all = TRUE
        ) %>%
        tidyr::nest() %>%
        dplyr::left_join(unite,
          by = "Seq"
        ) %>%
        dplyr::mutate(
          score = purrr::pmap(
            .l = list(
              data,
              unite
            ),
            .f = ~ annot_score(
              data,
              unite,
              org = org
            )
          )
        ) %>%
        dplyr::select(-unite) %>%
        tidyr::unnest(cols = -c(
          ID,
          Seq,
          AccID
        )) %>%
        dplyr::arrange(
          ID,
          dplyr::desc(score)
        )
    }
  } else {
    blst <- df_out
  }

  blst
}
