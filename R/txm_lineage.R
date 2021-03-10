#' Extract lineage
#'
#' Obtain full lineage from taxonomic IDs. Splits TaxIDs into batches of 200 and uses
#' \link[rentrez:entrez_post]{entrez_post} to post them onto the NCBI web history server, followed by using
#' \link[rentrez:entrez_fetch]{entrez_fetch} to obtain xml files.
#'
#' @param taxids (Required) Output table obtained from \link[taxminer]{txm_ecosrc}. Alternatively, a data.frame with
#'                a single column named "TaxID" and optionally a second column named "AccID".
#' @param bindtoAcc (Logical) Default TRUE. Bind lineage to accession IDs.
#'
#' @export

txm_lineage <- function(taxids, bindtoAcc = T) {

  if (!any(names(taxids) %in% "TaxId")) {
    stop("No column named 'TaxId' - please provide a valid input column")
  }

  if (bindtoAcc&!any(names(taxids) %in% "AccID")) {
    stop("No columns named 'AccID' while bindtoAcc is T - please provide a valid accession ID column")
  }


  ##### Lineage #####
  TaxID_to_src <- taxids %>%
    dplyr::distinct(TaxId)

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
    pb_assign$tick(1)
    TaxIDs_post <- rentrez::entrez_post(
      db = "taxonomy",
      id = paste(as.character(unlist(TaxID_to_src[[i]])),
                 collapse = ","
      )
    )

    TaxIDs_fetch <- rentrez::entrez_fetch(
      db = "taxonomy", web_history = TaxIDs_post,
      rettype = "xml", retmode = "xml", parsed = T
    )

    if (!class(try({XML::xmlToDataFrame(TaxIDs_fetch)}, silent = T)) == "try-error") {
      if (nrow(
        XML::xmlToDataFrame(TaxIDs_fetch)) == length(as.character(unlist(TaxID_to_src[[i]])))) {

        j <- 1
        while (j <= nrow(XML::xmlToDataFrame(TaxIDs_fetch))) {
          xml <- XML::xmlRoot(TaxIDs_fetch)
          xmlnodes <- XML::xmlSApply(xml[[j]], XML::xmlName)
          xmlsub <- which(xmlnodes == "LineageEx")
          xmlsubnode <- xml[[j]][[xmlsub]]
          nodeChildren <- XML::xmlChildren(xmlsubnode)

          names(nodeChildren) <- make.unique(as.character(purrr::map(
            .x = nodeChildren,
            .f = ~ XML::xmlValue(.x[[3]])
          )))


          lineage_step <- purrr::map_dfr(.x = nodeChildren, .f = ~ XML::xmlValue(.x[[2]])) %>%
            as.data.frame() %>%
            dplyr::mutate(taxid = XML::xmlValue(xml[[j]][[1]]))

          lineage <- lineage %>%
            dplyr::bind_rows(lineage_step)
          j <- j + 1
        }
        if (i >= round(num_groups)) {
          message("...Done!")
        }
        i <- i + 1
        Sys.sleep(1)
      }
    } else {
      Sys.sleep(1)
    }
  }

  lineage <- lineage %>%
    dplyr::rename("TaxId" = .data$taxid) %>%
    dplyr::select(.data$TaxId, .data$superkingdom, .data$phylum, .data$class,
                  .data$order, .data$family, .data$genus)

  if (bindtoAcc) {
    taxids <- taxids %>%
      dplyr::rename("species" = .data$Organism) %>%
      dplyr::left_join(lineage, by = "TaxId")
  } else {
    lineage
  }
}
