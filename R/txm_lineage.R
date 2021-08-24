#' Extract lineage
#'
#' Obtain full lineage from taxonomic IDs. Splits TaxIDs into batches of 200 and uses
#' \link[rentrez:entrez_post]{entrez_post} to post them onto the NCBI web history server, followed by using
#' \link[rentrez:entrez_fetch]{entrez_fetch} to obtain xml files.
#'
#' @param taxids (Required) Output table obtained from \link[taxminer]{txm_ecosrc}. Alternatively, a data.frame with
#'                a single column named "TaxID" and optionally a second column named "AccID".
#' @param bindtoAcc (Logical) Default TRUE. Bind lineage to accession IDs.
#' @param Precomp_tbl (Optional) Default NA. Specify the name of the pre-compiled database
#'                     present within the directory.
#' @param Precomp_tbl_assign (Optional) Default Dataset + system date. Name of a new compiled database
#'                            to be assigned
#' @param savedata (Logical) Default TRUE. Should a compiled database be saved to directory?
#'
#' @export

txm_lineage <- function(taxids, bindtoAcc = T, Precomp_tbl = NA,
                        Precomp_tbl_assign = paste("Dataset_", Sys.Date(), ".rds", sep = ""), savedata = T) {

  if (!any(names(taxids) %in% "TaxID")) {
    stop("No column named 'TaxID' - please provide a valid input column")
  }

  if (bindtoAcc&!any(names(taxids) %in% "AccID")) {
    stop("No columns named 'AccID' while bindtoAcc is T - please provide a valid accession ID column")
  }

  ##### Lineage #####
  TaxID_to_src <- taxids %>%
    dplyr::distinct(TaxID)

  ##### Read in pre-compiled database #####
  if (!is.na(Precomp_tbl)) {
    print("Reading in dataset and searching for existing lineage")
    Precomp_tbl_sub <- readRDS(Precomp_tbl) %>%
      dplyr::inner_join(TaxID_to_src, by = "TaxID")

    TaxID_to_src <- TaxID_to_src %>%
      dplyr::filter(!.data$TaxID %in% Precomp_tbl_sub$TaxID) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(.data$TaxID)
    print(paste(nrow(TaxID_to_src), "of",
                sum(nrow(Precomp_tbl_sub), nrow(TaxID_to_src)),
                "will be searched for",
                sep = " "
    ))
  }


  ##### Check whether there are any IDs left to search for #####
  if (nrow(TaxID_to_src) > 0) {
    if (is.na(Precomp_tbl)) {
      print(paste("No dataset provided -", nrow(TaxID_to_src),
                  "taxids will be searched for",
                  sep = " "
      ))
    }


  ##### Define split ####
  num_groups <- nrow(TaxID_to_src) / 200
  if (num_groups < 1) num_groups <- 1

  TaxID_to_src <- TaxID_to_src %>%
    dplyr::group_by((dplyr::row_number() - 1) %/% (dplyr::n() / num_groups)) %>%
    tidyr::nest() %>%
    dplyr::pull()

  i <- 1
  lineage <- data.frame()
  pb_assign <- progress::progress_bar$new(
    format = "  downloading [:bar] :current/:total (:percent) eta: :eta elapsed: :elapsed",
    total = round(num_groups), clear = FALSE, width= 60)
  print("Retrieving Lineage")
  while (i <= round(num_groups)) {
    chunk_lin <- try({
      # Post TaxIDs
      TaxIDs_post <- rentrez::entrez_post(
        db = "taxonomy",
        id = paste(as.character(unlist(TaxID_to_src[[i]])),
                   collapse = ","
        )
      )
      # Fetch full taxonomy from web history
      TaxIDs_fetch <- rentrez::entrez_fetch(
        db = "taxonomy", web_history = TaxIDs_post,
        rettype = "xml", retmode = "xml", parsed = T
      )
      j <- 1
      while (j <= nrow(XML::xmlToDataFrame(TaxIDs_fetch))) {
        xml_lineage <- XML::xmlRoot(TaxIDs_fetch) %>%
          .[[j]] %>%
          .[["LineageEx"]] %>%
          XML::xmlChildren() %>%
          purrr::set_names(make.unique(as.character(purrr::map(
            .x = .,
            .f = ~ XML::xmlValue(.x[[3]]))))
          ) %>%
          purrr::map_dfr(
            .x = .,
            .f = ~ XML::xmlValue(.x[[2]])
          ) %>%
          as.data.frame() %>%
          dplyr::mutate(
            TaxID = XML::xmlValue(XML::xmlRoot(TaxIDs_fetch)[[j]][[1]])
            ) %>%
          dplyr::mutate(
            species = XML::xmlValue(XML::xmlRoot(TaxIDs_fetch)[[j]][["ScientificName"]][[1]])
            ) %>%
          dplyr::select(TaxID, tidyselect::everything())

        lineage <- lineage %>%
          dplyr::bind_rows(xml_lineage)
        j <- j + 1
      }
    })
    if(!class(chunk_lin) == "try-error") {
      pb_assign$tick()
      i <- i + 1
      if (i > round(num_groups)) {
        message("...Done!")
      }
      Sys.sleep(1)
    } else {
      i <- i
      Sys.sleep(1)
    }
  }
  lineage <- lineage %>%
    dplyr::select(.data$TaxID, .data$superkingdom, .data$phylum, .data$class,
                  .data$order, .data$family, .data$genus, .data$species)

  if (savedata) {
    if (!is.null(lineage)) {
      if (!is.na(Precomp_tbl)) {
        if (!is.null(TaxID_to_src)) {
          print("Writing new data to Pre-compiled dataset")
          Precomp_full <- readRDS(Precomp_tbl) %>%
            dplyr::bind_rows(lineage) %>%
            dplyr::arrange(TaxID) %>%
            saveRDS(file = Precomp_tbl)
        }
      } else {
        print("Creating dataset")
        saveRDS(lineage, file = Precomp_tbl_assign)
      }
    }
  }
  if (!is.na(Precomp_tbl)) {
    lineage <- lineage %>%
      dplyr::bind_rows(Precomp_tbl_sub)
  }
  } else {
    lineage <- Precomp_tbl_sub
  }


  if (bindtoAcc) {
    taxids <- taxids %>%
      dplyr::rename_with(~paste("species"), tidyselect::starts_with("Organism")) %>%
      dplyr::select(-tidyselect::starts_with("species")) %>%
      dplyr::mutate(TaxID = as.numeric(TaxID)) %>%
      dplyr::inner_join(lineage, by = "TaxID")
  } else {
    lineage <- lineage
  }
}
