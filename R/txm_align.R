utils::globalVariables(c(
  "ID",
  ".",
  "owrt",
  "owrt_silva",
  "owrt_RDP",
  "ASVs",
  "bitscore",
  "Pct",
  "qcovs",
  "Evalue"
))

#' Customized local BLAST
#'
#' Build and run a customized BLAST alignment using local cmd line BLAST
#' and locally downloaded databases from within RStudio. User defined seq_ins
#' and specifications are used to build a BLAST cmd, which is transferred
#' to RStudio terminal. The output is formatted as BLASTn tabular output
#' format 6.
#' In parallel, SILVA and RDP databases will be used to assign alternative
#' annotations that are further used for annotations scores.
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
#' @param alt_annot (Logical) Default TRUE. Alternative annotations with SILVA
#' and RDP databases using \link[dada2]{assignTaxonomy} for up to genus level
#' annotations followed by \link[dada2]{addSpecies} for assigning species.
#' @param alt_path (Required) Default NULL. Path to SILVA and RDP databases.
#' Must identify a folder that contains at least 4 files:
#' \itemize{
#'   \item silva training set (regex: "silva.*train")
#'   \item silva species assignment (regex: "silva.*species")
#'   \item rdp train set (regex: "rdp.*train")
#'   \item rdp species assignment (regex: "rdp.*species")
#' }
#' Pre-formatted dada2 databases are recommended:
#' \href{https://benjjneb.github.io/dada2/training.html}{dada2 fasta files}
#' @param asbd_tbl_lge (Optional) Default NA. Full path and name of a
#' pre-compiled table to be used for lineage assignment
#' @param asgn_tbl_lge (Optional) Default Dataset_lge + system date. Full path
#' and name of a new data set to be created.
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

txm_align <- function(seq_in,
                      db_path = NULL,
                      db_name = NULL,
                      tab_out = paste("Output", Sys.Date(), sep = ""),
                      tab_path = ".",
                      acsn_list = NULL,
                      acsn_path = ".",
                      alt_path = NULL,
                      alt_annot = T,
                      asbd_tbl_lge = NA,
                      asgn_tbl_lge = paste("Dataset_lge",
                        Sys.Date(),
                        ".fst",
                        sep = ""
                      ),
                      task = "megablast",
                      threads = 1,
                      acsn_check = F,
                      show = F,
                      run_blst = T,
                      qcvg = 98,
                      pctidt = 98,
                      max_out = 500) {

  ##### Data prep -----
  check_align(
    db_name,
    db_path,
    task,
    acsn_path,
    acsn_list,
    tab_out,
    tab_path,
    run_blst,
    acsn_check,
    alt_annot,
    alt_path
  )

  mk_fasta(
    seq_in,
    tab_path,
    tab_out
  )

  if (run_blst & owrt) {

    ##### Create BLAST command and run -----
    blst_cmd <- paste("blastn -task ",
      task,
      sep = ""
    )

    db <- paste("-db ",
      here::here(
        db_path,
        db_name
      ),
      sep = ""
    )

    query <- paste("-query ",
      here::here(
        tab_path,
        paste(tab_out,
          ".fa",
          sep = ""
        )
      ),
      sep = ""
    )

    output <- paste("-out ",
      here::here(
        tab_path,
        paste("Alignment_",
          tab_out,
          ".csv",
          sep = ""
        )
      ),
      sep = ""
    )

    parameters <- paste("-num_threads",
      threads,
      "-perc_identity",
      pctidt,
      sep = " "
    )

    accession_limit <- paste("-seqidlist ",
      here::here(
        acsn_path,
        acsn_list
      ),
      sep = ""
    )

    output_format <- paste(
      "-max_target_seqs",
      max_out,
      "-outfmt '6 qacc sgi saccver staxids sscinames bitscore qcovs evalue pident'"
    )

    trmnl_cmd <- paste(blst_cmd,
      db,
      query,
      output,
      if (!is.null(acsn_list)) {
        accession_limit
      },
      parameters,
      output_format,
      sep = " "
    )

    Blast <- rstudioapi::terminalExecute(trmnl_cmd,
      show = show
    )

    print(paste("Running Blast - ",
      trmnl_cmd,
      sep = ""
    ))

    while (is.null(rstudioapi::terminalExitCode(Blast))) {
      Sys.sleep(0.1)
    }

    if (rstudioapi::terminalExitCode(Blast) == 0) {
      print("Blast successful")
    } else {
      stop(print(
        paste(
          "Blast ran into an error - Code: ",
          rstudioapi::terminalExitCode(Blast),
          " - Please check terminal for details"
        )
      ))
    }
    rstudioapi::terminalKill(Blast)
  }

  ##### Alternative annotations (Silva + RDP) -----
  if (alt_annot) {
    list_dbs <- list.files(alt_path)
    set.seed(10)

    if (owrt_silva) {
      print("assigning taxonomies from silva")
      annot_silva_sp <- dada2::assignTaxonomy(
        seq_in,
        here::here(
          alt_path,
          stringr::str_subset(
            list_dbs,
            "silva.*train"
          )
        ),
        minBoot = 80,
        multithread = T,
        tryRC = T
      ) %>%
        dada2::addSpecies(
          here::here(
            alt_path,
            stringr::str_subset(
              list_dbs,
              "silva.*species"
            )
          ),
          allowMultiple = T,
          tryRC = T
        ) %>%
        data.frame() %>%
        tibble::rownames_to_column(
          var = "ASVs"
        ) %>%
        dplyr::group_by(ASVs) %>%
        tidyr::nest() %>%
        dplyr::mutate(data = purrr::map(
          data,
          rplc
        )) %>%
        tidyr::unnest(cols = data) %>%
        dplyr::mutate(
          # modify upper lineage for select taxa to align with 
          # current NCBI taxonomic levels to facilitate annotation scores. 
          # Potentially update to TaxIDs 
          Phylum = dplyr::case_when(
            Phylum == "Actinobacteriota" ~ "Actinobacteria",
            Phylum == "Fusobacteriota" ~ "Fusobacteria",
            Phylum == "Bacteroidota" ~ "Bacteroidetes", 
            TRUE ~ Phylum
          ), 
          Class = dplyr::case_when(
            Class == "Actinobacteria" ~ "Actinomycetia",
            TRUE ~ Class
          )
        )

      fst::write_fst(annot_silva_sp,
        path = here::here(
          tab_path,
          "Silva_annot.fst"
        )
      )
    }

    if (owrt_RDP) {
      print("assigning taxonomies from rdp")
      annot_rdp_sp <- dada2::assignTaxonomy(
        seq_in,
        here::here(
          alt_path,
          stringr::str_subset(
            list_dbs,
            "rdp.*train"
          )
        ),
        minBoot = 80,
        multithread = T,
        tryRC = T
      ) %>%
        dada2::addSpecies(
          here::here(
            alt_path,
            stringr::str_subset(
              list_dbs,
              "rdp.*species"
            )
          ),
          allowMultiple = T,
          tryRC = T
        ) %>%
        data.frame() %>%
        tibble::rownames_to_column(
          var = "ASVs") %>%
        dplyr::group_by(ASVs) %>%
        tidyr::nest() %>%
        dplyr::mutate(data = purrr::map(
          data,
          rplc
        )) %>%
        tidyr::unnest(cols = data) %>%
        # modify upper lineage for select taxa to align with 
        # current NCBI taxonomic levels to facilitate annotation scores. 
        # Potentially update to TaxIDs 
        dplyr::mutate( 
          Class = dplyr::case_when(
            Class == "Actinobacteria" ~ "Actinomycetia", 
            TRUE ~ Class
          )
        )

      fst::write_fst(annot_rdp_sp,
        path = here::here(
          tab_path,
          "RDP_annot.fst"
        )
      )
    }
  }

  rm(
    list = stringr::str_subset(
      ls(pos = .GlobalEnv),
      "owrt"
    ),
    envir = .GlobalEnv
  )

  ##### Process alignment results -----
  align_path <- here::here(
    tab_path,
    paste("Alignment_",
      tab_out,
      ".csv",
      sep = ""
    )
  )

  if (nrow(
    data.table::fread(align_path)
  ) > 0) {
    outtab <- data.table::fread(
      paste(here::here(
        tab_path,
        paste("Alignment_",
          tab_out,
          ".csv",
          sep = ""
        )
      ))
    ) %>%
      magrittr::set_colnames(c(
        "ID",
        "GiID",
        "AccID",
        "TaxID",
        "Species",
        "bitscore",
        "qcovs",
        "Evalue",
        "Pct"
      )) %>%
      dtplyr::lazy_dt() %>%
      dplyr::filter(.data$qcovs >= qcvg) %>%
      dplyr::mutate(
        TaxID = as.numeric(
          stringr::str_extract(
            .data$TaxID,
            "^[^;]*"
          )
        ),
        Species = stringr::str_extract(
          .data$Species,
          "^[^;]*"
        )
      ) %>%
      dplyr::left_join(ASVs) %>%
      dplyr::group_by(ID, AccID) %>%
      dplyr::filter(bitscore == max(bitscore)) %>%
      dplyr::filter(Pct == max(Pct)) %>%
      dplyr::filter(qcovs == max(qcovs)) %>%
      dplyr::filter(Evalue == min(Evalue)) %>%
      dplyr::distinct(AccID,
        .keep_all = T
      ) %>%
      dplyr::group_by(ID) %>%
      dplyr::add_tally() %>%
      dplyr::filter(n > 1) %>%
      dplyr::ungroup() %>%
      base::as.data.frame() %>%
      dplyr::select(
        .data$ID,
        .data$ASVs,
        .data$TaxID,
        .data$AccID,
        .data$GiID,
        tidyselect::everything()
      ) %>%
      txm_lineage(
        taxids = .,
        asbd_tbl = asbd_tbl_lge,
        asgn_tbl = asgn_tbl_lge
      )
  } else {
    print("No alignments")
  }
}
