#' Fetch accession ID list
#'
#' Fetcha a list of accession IDs that will be used for restricting subsequent BLAST alignments. This uses \link[rentrez]{rentrez}
#' to communicate with NCBI databases. Specifically, \link[rentrez:entrez_search]{entrez_search} to obtain all relevant UIDs for
#' the user defined "Text query", with the "use_history" option set to TRUE, thus posting all UIDs onto the
#' NCBI web history server. Next, \link[rentrez:entrez_fetch]{entrez_fetch} is used to extract accession IDs in batches of 10,000, appending each cycle to an
#' output file.
#'
#' @name txm_accIDs
#' @param text_query (Required) Default NA. Text query specifying what should be looked for.
#' @param db2src (Required) Default "nucleotide". Name of NCBI database to search.
#' @param out_name (Required) Default NA. Name of the output file.
#' @examples
#' \dontrun{
#' taxminer::txm_accIDs(text_query = "Fungi [organism] AND vagina)",
#'                      db2src = "nucleotide",
#'                      out_name = "Fungi_noEnv.seq")
#'}
#' @return A list of accession IDs, one per line, corresponding to the user defined text query and NCBI database
#' @importFrom magrittr %>%
#' @export

txm_accIDs <- function(
  text_query = NA,
  db2src = "nucleotide",
  out_name = NA
) {


  ##### Search for IDs and store in web history #####

  AccIDs_esearch <- rentrez::entrez_search(db = db2src,
                                           term = text_query,
                                           use_history = T)


  ##### Divided into chunks of 10,000 to fetch in batches #####
  total_ids <- AccIDs_esearch$count
  message(paste("Found ", AccIDs_esearch$count, " IDs", sep = ""))
  start_chunk <- 0
  pb_assign <- progress::progress_bar$new(
    format = "  downloading [:bar] :current/:total (:percent) eta: :eta elasped: :elapsed",
    total = total_ids/10000, clear = FALSE, width= 60)
  while (start_chunk <= total_ids) {
    seqname <- out_name
    return_chunk <- try({
      returns <- rentrez::entrez_fetch(db = "nuccore", web_history = AccIDs_esearch$web_history,
                                     rettype = "acc",
                                     retstart = start_chunk, retmax = 10000)
      readr::write_lines(returns, seqname, append = T, sep = "")
    })
    if (!class(return_chunk) == "try-error") {
      start_chunk <- start_chunk + 10000
      pb_assign$tick()
    } else {
      Sys.sleep(1)
      start_chunk <- start_chunk
    }
    if (start_chunk >= total_ids) {
      message(paste("Done!"))
    }
  }
}
