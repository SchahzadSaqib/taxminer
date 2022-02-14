utils::globalVariables(c(
  "ID", "."
))

#' Customized local BLAST
#'
#' Build and run a customized BLAST alignment using local cmd line BLAST 
#' and locally downloaded databases from within RStudio. User defined seq_ins 
#' and specifications are used to build a BLAST cmd, which is transferred 
#' to RStudio terminal. The output is formatted as BLASTn tabular output 
#' format 6.
#'
#' @name txm_align
#' @param seq_in (Required) All seq_in types accepted by 
#' \link[dada2]{getSequences}.
#' @param db_path (Required) Full path for local BLAST database being used for 
#' alignment.
#' @param db_name (Required) Name of BLAST database.
#' @param task (Optional) Default "megablast". Other possible tasks can be 
#' found here: 
#' \href{https://www.ncbi.nlm.nih.gov/books/NBK569839/#usrman_BLAST_feat.Tasks}{Tasks}
#' @param tab_out (Optional) Default "Output" + System date. Specify name of 
#' the output files.
#' @param tab_path (Optional) Default "." (current working directory). Specify 
#' full output path for results.
#' @param threads (Optional) Default 1. Number of threads assigned to BLAST 
#' alignment.
#' @param acsn_list (Optional) Default NULL. Full path to a list of accession 
#' IDs that will be used to restrict the BLAST database (-seqidlist), which 
#' is highly recommended for large databases such as "nt/nr". This can be 
#' obtained using the \link[taxminer]{txm_accIDs} function provided within this
#' package. Further information regarding search limitations is available on
#' \href{https://www.ncbi.nlm.nih.gov/books/NBK279673/}{BLAST user manual}
#' @param acsn_path (Optional) Default "." (current working directory). Specify 
#' full path to the accession ID list.
#' @param acsn_check (Logical) Default FALSE. If an accession ID list is 
#' provided, BLASTDB v5 requires it to be pre-processed (blastdb_aliastool) 
#' prior to being used for restricting the database. Set this to TRUE if an 
#' unprocessed accession ID list is specified.
#' @param show (Logical) Default FALSE. Switch from console to terminal?
#' @param run_blst (Logical) Default TRUE. Set to FALSE if an existing 
#' comma-delimited alignment file is present in the directory.
#' @param qcvg (Optional) Default 98 (%). Threshold for query coverage. 
#' @param pctidt (Optional) Default 98 (%). Threshold for percentage identity.
#' @param max_out (Optional) Default 500. Maximum number of alignments 
#' ("-max_target_seqs"). Higher values are recommended when using large 
#' databases such as "nt/nr".
#' @export
#' @importFrom rlang .data
#' @importFrom data.table fread

txm_align <- function(
  seq_in,
  db_path = NULL,
  db_name = NULL,
  tab_out = paste("Output", Sys.Date(), sep = ""),
  tab_path = ".",
  acsn_list = NULL,
  acsn_path = ".",
  task = "megablast",
  threads = 1,
  acsn_check = F,
  show = F,
  run_blst = T,
  qcvg = 98,
  pctidt = 98,
  max_out = 500
) {
  
  ##### Data prep -----
  check_align(db_name, 
               db_path, 
               task, 
               acsn_path, 
               acsn_list, 
               tab_out, 
               tab_path)
  
  ASVs <- dada2::getSequences(seq_in) %>%
    base::as.data.frame() %>%
    purrr::set_names("ASVs") %>%
    dplyr::mutate(ID = 1:nrow(.)) %>% 
    dplyr::select(.data$ID, .data$ASVs)
  
  print(paste(nrow(ASVs), 
              " will be aligned", 
              sep = ""))
  
  FASTA_file <- ASVs %>%
    dplyr::mutate(ID = paste(">", ID, sep = "")) %>%
    base::as.list() %>%
    purrr::pmap(~ paste(.x, .y, sep = "\n")) %>%
    base::unlist()
  
  readr::write_lines(
    FASTA_file, 
    file = paste(
      here::here(tab_path, 
                 tab_out), 
      ".fa", 
      sep = ""))
  
  if (run_blst & 
      file.exists(
        here::here(tab_path,
                   paste("Alignment_",
                         tab_out,
                         ".csv",
                         sep = "")))) {
    overwrite <- utils::askYesNo(
      paste(
        basename(
          here::here(tab_path,
                     paste("Alignment_",
                           tab_out,
                           ".csv",
                           sep = ""))), 
        "Overwrite existing file?", 
        sep = " "))
  } else {
    overwrite <- T
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
                overwrite = T)
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
    
    ##### Accession ID prep -----
    if (!is.null(acsn_list)) {
      if (acsn_check == T) {
        cmd <- "blastdb_aliastool"
        seq_in_list <- paste("-seqid_file_in ", 
                             here::here(acsn_path,
                                        acsn_list), 
                             sep = "")
        
        trmnl_cmd <- paste(cmd,
                           seq_in_list,
                           sep = " ")
        
        acc_check <- rstudioapi::terminalExecute(trmnl_cmd, 
                                                 show = show)
        
        while (is.null(rstudioapi::terminalExitCode(acc_check))) {
          Sys.sleep(0.1)
        }
        if (rstudioapi::terminalExitCode(acc_check) == 0) {
          print("Accession list check successful")
          rstudioapi::terminalKill(acc_check)
          acsn_list <- paste(acsn_list, 
                             ".bsl", 
                             sep = "")
        } else {
          stop(print(
            paste("Accession list check ran into an error - Code: ",
                  rstudioapi::terminalExitCode(acc_check),
                  " - Please check terminal for details")))
        }
      } else {
        acsn_list <- paste(acsn_list, 
                           ".bsl", 
                           sep = "")
        if (!file.exists(
          here::here(acsn_path, acsn_list))) {
          stop(print(
            paste("accession list not found in specified directory")))
        }
      }
    }
    
    if (file.size(
      paste(
        here::here(tab_path, 
                   tab_out), 
        ".fa", 
        sep = "")) == 0) {
      stop(print("Empty query object provided"))
    }
    
    ##### Create BLAST command and run -----
    blst_cmd <- paste("blastn -task ", 
                      task, 
                      sep = "")
    db <- paste("-db ", 
                here::here(db_path, 
                           db_name), 
                sep = "")
    query <- paste("-query ", 
                   here::here(tab_path,
                              paste(tab_out, 
                                    ".fa", 
                                    sep = "")), 
                   sep = "")
    output <- paste("-out ", 
                    here::here(tab_path,
                               paste("Alignment_",
                                     tab_out,
                                     ".csv",
                                     sep = "")), 
                    sep = "")
    parameters <- paste("-num_threads", 
                        threads, 
                        "-perc_identity", 
                        pctidt, 
                        sep = " ")
    accession_limit <- paste("-seqidlist ", 
                             here::here(acsn_path, 
                                        acsn_list), 
                             sep = "")
    output_format <- paste(
      "-max_target_seqs",
      max_out,
      "-outfmt '6 qacc sseqid staxids sscinames bitscore qcovs evalue pident'")
    
    trmnl_cmd <- paste(blst_cmd,
                       db,
                       query,
                       output,
                       if (!is.null(acsn_list)) {
                         accession_limit
                       },
                       parameters,
                       output_format,
                       sep = " ")
    
    Blast <- rstudioapi::terminalExecute(trmnl_cmd, 
                                         show = show)
    print(paste("Running Blast - ", 
                trmnl_cmd, 
                sep = ""))
    while (is.null(rstudioapi::terminalExitCode(Blast))) {
      Sys.sleep(0.1)
    }
    
    if (rstudioapi::terminalExitCode(Blast) == 0) {
      print("Blast successful")
    } else {
      stop(print(
        paste("Blast ran into an error - Code: ", 
              rstudioapi::terminalExitCode(Blast),
              " - Please check terminal for details")))
    }
    rstudioapi::terminalKill(Blast)
  }
  
  ##### Process alignment results -----
  align_path <- here::here(tab_path, 
                           paste("Alignment_", 
                                 tab_out,
                                 ".csv", 
                                 sep = ""))
  
  if (nrow(
    data.table::fread(align_path)) > 0) {
    outtab <- data.table::fread(
      paste(here::here(tab_path, 
                       paste("Alignment_", 
                             tab_out,
                             ".csv", 
                             sep = "")))) %>%
      magrittr::set_colnames(c("ID", 
                               "SeqID", 
                               "TaxID", 
                               "Species",
                               "bitscore", 
                               "qcovs", 
                               "Evalue", 
                               "Pct")) %>% 
      dtplyr::lazy_dt() %>%
      dplyr::filter(.data$qcovs >= qcvg) %>%
      dplyr::mutate(
        SeqID = stringr::str_replace_all(
          .data$SeqID,
          pattern = "\\|", 
          replacement = ";"), 
        AccID = stringr::str_extract(
          .data$SeqID, 
          "(?<=gb;).*(?=;)"),
        GiID = stringr::str_extract(
          .data$SeqID,
          "(?<=gi;).*(?=;gb)"),
        TaxID = as.numeric(
          stringr::str_extract(
            .data$TaxID,
            "^[^;]*")),
        Species = stringr::str_extract(
          .data$Species, 
          "^[^;]*")) %>%
      base::as.data.frame() %>%
      dplyr::select(.data$ID, 
                    .data$TaxID, 
                    .data$AccID, 
                    .data$GiID, 
                    tidyselect::everything(), 
                    -.data$SeqID)
  } else {
    print("No alignments")
  }
}
