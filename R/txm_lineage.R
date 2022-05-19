utils::globalVariables(c(
  "TaxID",
  "norank.1",
  "species"
))

#' Extract lineage
#'
#' Obtain full lineage from taxonomic IDs. Splits TaxIDs into batches of 200 and uses
#' \link[rentrez:entrez_post]{entrez_post} to post them onto the NCBI
#' web history server, followed by using
#' \link[rentrez:entrez_fetch]{entrez_fetch} to obtain xml files.
#'
#' @param taxids (data.frame) **Required** Output table obtained
#' from \link[taxminer]{txm_ecosrc}. Alternatively, a data.frame with a single
#' column named "TaxID" and optionally a second column named "AccID".
#' @param asbd_tbl (String) Default NA. Specify the name of the
#' pre-assembled database present within the directory.
#' @param asgn_tbl (String) Default Dataset_lge + system date. Name of a new
#' database to be assembled
#'
#' @export

txm_lineage <- function(taxids,
                        asbd_tbl = NA,
                        asgn_tbl = paste("Dataset_lge_",
                          Sys.Date(),
                          ".fst",
                          sep = ""
                        )) {
  check_lineage(
    taxids,
    asbd_tbl,
    asgn_tbl
  )


  # Data preparation -----
  taxids_src <- taxids %>%
    dplyr::distinct(TaxID)

  if (!is.na(asbd_tbl)) {
    print("Reading in dataset and searching for existing lineage")
    asbd_tbl_sub <- fst::read_fst(asbd_tbl) %>%
      dplyr::inner_join(taxids_src, by = "TaxID")

    taxids_src <- taxids_src %>%
      dplyr::filter(!.data$TaxID %in% asbd_tbl_sub$TaxID) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(.data$TaxID)
    print(paste(nrow(taxids_src), "of",
      sum(nrow(asbd_tbl_sub), nrow(taxids_src)),
      "will be searched for",
      sep = " "
    ))
  }

  if (nrow(taxids_src) > 0) {
    if (is.na(asbd_tbl)) {
      print(paste("No dataset provided -", nrow(taxids_src),
        "taxids will be searched for",
        sep = " "
      ))
    }


    ## split data into chunks -----
    splits <- round(nrow(taxids_src) / 200)
    if (splits < 1) splits <- 1

    taxids_src <- taxids_src %>%
      dplyr::mutate(chunks = rep(0:splits,
        each = 200,
        length.out = nrow(.)
      )) %>%
      dplyr::group_by(.data$chunks) %>%
      tidyr::nest() %>%
      dplyr::pull()

    prgrs_bar <- new_bar(length(taxids_src))
    assign(
      "prgrs_bar",
      prgrs_bar,
      .GlobalEnv
    )

    lge <- c(
      "superkingdom",
      "phylum",
      "class",
      "order",
      "family",
      "genus",
      "species"
    ) %>%
      purrr::map_dfr(~ tibble(!!.x := logical()))

    print("Retrieving Lineage")


    ## extract lineage -----
    get <- taxids_src %>%
      purrr::lmap(~ get_lge(.x))

    clean <- get %>%
      purrr::flatten() %>%
      purrr::map_dfr(.f = dplyr::bind_cols) %>%
      dplyr::bind_rows(lge[setdiff(names(lge), names(.))]) %>%
      dplyr::select(TaxID, names(lge), everything()) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        norank.1 = ifelse(
          "norank.1" %in% names(.),
          norank.1,
          NA
        ),
        dplyr::across(
          .cols = tidyselect::everything(),
          .f = ~ ifelse(stringr::str_detect(
            .x,
            "environmental samples"
          ),
          NA,
          .x
          )
        ),
        norank.1 = ifelse(is.na(norank.1),
          .data$species,
          norank.1
        ),
        norank.1 = ifelse(stringr::str_starts(norank.1, "bacterium|unclassified"),
          NA,
          norank.1
        )
      ) %>%
      dplyr::select(-tidyselect::contains("species")) %>%
      dplyr::rename("species" = norank.1) %>%
      dplyr::select(tidyselect::all_of(
        c(
          "TaxID",
          names(lge)
        )
      )) %>%
      dplyr::mutate(
        species = stringr::str_extract(
          species,
          "^[^ ]* [^ ]*"
        )
      )

    rearrange <- clean %>%
      tidyr::pivot_longer(-TaxID,
        names_to = "level",
        values_to = "names"
      ) %>%
      dplyr::group_by(TaxID) %>%
      dplyr::mutate(names = ifelse(is.na(names),
        paste(
          "unclassified",
          dplyr::nth(rev(stats::na.omit(names)), 1)
        ),
        names
      )) %>%
      tidyr::pivot_wider(
        id_cols = TaxID,
        names_from = .data$level,
        values_from = .data$names
      )

    lineage <- rearrange %>%
      dplyr::ungroup() %>%
      dplyr::mutate(TaxID = as.numeric(TaxID)) %>%
      as.data.frame() %>%
      dplyr::mutate(
        dplyr::across(
          .cols = -TaxID,
          ~ as.character(.x)
        )
      ) %>%
      dplyr::mutate(
        species = stringr::str_remove_all(
          species,
          "'"
        )
      )


    ## save -----
    if (!is.null(lineage)) {
      if (!is.na(asbd_tbl)) {
        if (!is.null(taxids_src)) {
          rm(prgrs_bar, envir = .GlobalEnv)
          print("Writing new data to Pre-compiled dataset")
          asbd_full <- fst::read_fst(asbd_tbl) %>%
            dplyr::bind_rows(lineage) %>%
            dplyr::arrange(.data$TaxID) %>%
            dplyr::distinct(.data$TaxID,
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
          lineage,
          asgn_tbl,
          100
        )
      }
    }

    if (!is.na(asbd_tbl)) {
      lineage <- lineage %>%
        dplyr::bind_rows(asbd_tbl_sub) %>%
        dplyr::distinct(TaxID,
          .keep_all = T
        )
    }
  } else {
    lineage <- asbd_tbl_sub
  }


  ## bind to accessions -----
  if (any(names(taxids) == "Species")) {
    df_out <- taxids %>%
      dplyr::select(-.data$Species) %>%
      dplyr::left_join(lineage, by = "TaxID")
  } else {
    df_out <- lineage
  }
}
