utils::globalVariables(c(
  "TaxID"
))

#' Extract lineage
#'
#' Obtain full lineage from taxonomic IDs. Splits TaxIDs into batches of 200 and uses
#' \link[rentrez:entrez_post]{entrez_post} to post them onto the NCBI 
#' web history server, followed by using
#' \link[rentrez:entrez_fetch]{entrez_fetch} to obtain xml files.
#'
#' @param taxids (Required) Output table obtained 
#' from \link[taxminer]{txm_ecosrc}. Alternatively, a data.frame with a single 
#' column named "TaxID" and optionally a second column named "AccID".
#' @param bindtoAcc (Logical) Default TRUE. Bind lineage to accession IDs.
#' @param asbd_tbl (Optional) Default NA. Specify the name of the 
#' pre-assembled database present within the directory.
#' @param asgn_tbl (Optional) Default Dataset + system date. Name of a new 
#' database to be assembled
#' @param savedata (Logical) Default TRUE. Should a compiled database be 
#' saved to directory?
#'
#' @export

txm_lineage <- function(
  taxids, 
  bindtoAcc = T, 
  asbd_tbl = NA,
  asgn_tbl = paste("Dataset_", Sys.Date(), ".rds", sep = ""), 
  savedata = T) {
  
  check_lineage(taxids, 
               bindtoAcc, 
               asbd_tbl)

  ##### Data preparation -----
  taxids_src <- taxids %>%
    dplyr::distinct(TaxID)

  if (!is.na(asbd_tbl)) {
    print("Reading in dataset and searching for existing lineage")
    asbd_tbl_sub <- readRDS(asbd_tbl) %>%
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

  splits <- round(nrow(taxids_src) / 200)
  if (splits < 1) splits <- 1

  taxids_src <- taxids_src %>%
    dplyr::mutate(chunks = rep(0:splits, 
                               each = 200, 
                               length.out = nrow(.))) %>%
    dplyr::group_by(.data$chunks) %>%
    tidyr::nest() %>%
    dplyr::pull()

  ##### extract lineage -----
  lge <- 1
  lineage <- data.frame()
  pb_assign <- progress::progress_bar$new(
    format = paste("  downloading [:bar]", 
                   ":current/:total", 
                   "(:percent) eta:", 
                   ":eta elapsed:", 
                   ":elapsed", 
                   sep = " "),
    total = splits, 
    clear = FALSE, 
    width = 60)
  print("Retrieving Lineage")
  while (lge <= round(splits)) {
    chunk_lin <- try({
      TaxIDs_post <- rentrez::entrez_post(
        db = "taxonomy",
        id = paste(as.character(unlist(taxids_src[[lge]])),
                   collapse = ","
        )
      )
      
      TaxIDs_fetch <- rentrez::entrez_fetch(
        db = "taxonomy", web_history = TaxIDs_post,
        rettype = "xml", retmode = "xml", parsed = T
      )
      
      xml_step <- 1
      while (xml_step <= nrow(XML::xmlToDataFrame(TaxIDs_fetch))) {
        xml_lineage <- XML::xmlRoot(TaxIDs_fetch) %>%
          .[[xml_step]] %>%
          .[["LineageEx"]] %>%
          XML::xmlChildren() %>%
          purrr::set_names(
            make.unique(
              as.character(
                purrr::map(
            .x = .,
            .f = ~ XML::xmlValue(.x[[3]]))))
          ) %>%
          purrr::map_dfr(
            .x = .,
            .f = ~ XML::xmlValue(.x[[2]])
          ) %>%
          as.data.frame() %>%
          dplyr::mutate(
            TaxID = XML::xmlValue(
              XML::xmlRoot(TaxIDs_fetch)[[xml_step]][[1]])
            ) %>%
          dplyr::mutate(
            species = XML::xmlValue(
              XML::xmlRoot(TaxIDs_fetch)[[xml_step]][["ScientificName"]][[1]])
            ) %>%
          dplyr::select(.data$TaxID, 
                        tidyselect::everything())

        lineage <- lineage %>%
          dplyr::bind_rows(xml_lineage)
        xml_step <- xml_step + 1
      }
    })
    if(!class(chunk_lin) == "try-error") {
      pb_assign$tick()
      lge <- lge + 1
      if (lge > round(splits)) {
        message("...Done!")
      }
      Sys.sleep(1)
    } else {
      lge <- lge
      Sys.sleep(1)
    }
  }
  lineage <- lineage %>%
    dplyr::mutate(TaxID = as.numeric(TaxID)) %>%
    dplyr::select(tidyselect::any_of(c("TaxID", 
                                       "superkingdom", 
                                       "phylum",
                                       "class",
                                       "order",
                                       "family",
                                       "genus",
                                       "species")))

  if (savedata) {
    if (!is.null(lineage)) {
      if (!is.na(asbd_tbl)) {
        if (!is.null(taxids_src)) {
          print("Writing new data to Pre-compiled dataset")
          asbd_full <- readRDS(asbd_tbl) %>%
            dplyr::bind_rows(lineage) %>%
            dplyr::arrange(.data$TaxID) %>%
            saveRDS(file = asbd_tbl)
        }
      } else {
        print("Creating dataset")
        saveRDS(lineage, file = asgn_tbl)
      }
    }
  }
  if (!is.na(asbd_tbl)) {
    lineage <- lineage %>%
      dplyr::bind_rows(asbd_tbl_sub)
  }
  } else {
    lineage <- asbd_tbl_sub
  }


  ##### bind to accessions -----
  if (bindtoAcc) {
    taxids <- taxids %>%
      dplyr::rename_with(~paste("species"), 
                         tidyselect::starts_with("Organism")) %>%
      dplyr::select(-tidyselect::starts_with("species")) %>%
      dplyr::mutate(TaxID = as.numeric(TaxID)) %>%
      dplyr::inner_join(lineage, by = "TaxID")
  } else {
    lineage <- lineage
  }
}
