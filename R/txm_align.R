utils::globalVariables(c(
  "ID",
  ".",
  "n",
  "owrt",
  "owrt_silva",
  "owrt_RDP",
  "owrt_unite",
  "ASVs",
  "bitscore",
  "Pct",
  "family",
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
#' @param seq_in **Required** All input types accepted by
#' \link[dada2]{getSequences}.
#' @param db_path (String) **Required** Full path for local BLAST database being used for
#' alignment.
#' @param db_name (String) **Required** Name of BLAST database.
#' @param task (String) Default "megablast". Other possible tasks can be
#' found here:
#' \href{https://www.ncbi.nlm.nih.gov/books/NBK569839/#usrman_BLAST_feat.Tasks}{Tasks}
#' @param tab_out (String) Default "Output" + System date. Specify name of
#' the output files.
#' @param tab_path (String) Default "." (current working directory). Specify
#' full output path for results.
#' @param threads (>0) Default 1. Number of threads assigned to BLAST
#' alignment. Set these in accordance to your hardware.
#' @param acsn_list (String) Default NULL. Full path to a list of accession
#' IDs that will be used to restrict the BLAST database (-seqidlist), which
#' is highly recommended for large databases such as "nt/nr". This can be
#' obtained using the \link[taxminer]{txm_accIDs} function provided within this
#' package. Further information regarding search limitations is available on
#' \href{https://www.ncbi.nlm.nih.gov/books/NBK279673/}{BLAST user manual}
#' @param acsn_path (String) Default "." (current working directory). Specify
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
#'   \item UNITE general release fasta (regex: "sh_general.*fasta")
#' }
#' Pre-formatted dada2 databases are recommended:
#' \href{https://benjjneb.github.io/dada2/training.html}{dada2 fasta files}
#' @param org (String) Default "bac". The organism being investigated.
#' \itemize{
#'   \item "bac" - Bacteria. Silva and RDP databases are used for alternative
#'   annotations to BLAST.
#'   \item "fungi" - Fungi. UNITE database is used for alternative annotations
#'   to BLAST.
#' }
#' @param asbd_tbl_lge (String) Default NA. Full path and name of a
#' pre-compiled table to be used for lineage assignment
#' @param asgn_tbl_lge (String) Default Dataset_lge + system date. Full path
#' and name of a new data set to be created.
#' @param show (Logical) Default FALSE. Switch from console to terminal?
#' @param run_blst (Logical) Default TRUE. Run BLAST alignment?
#' @param batches (>0) Default 1. Should the sequences be aligned in batches?
#' \itemize{
#'   \item 1 = no batches
#'   \item >1 = total sequences/batches
#' }
#' @param chunks (>0) Default 1. Should each batch of sequences be split into
#' further chunks to be aligned in parallel?
#' \itemize{
#'   \item 1 = no further chunks
#'   \item >1 = each batch is further divided to the closest number of chunks
#'   possible and ran in parallel using GNU parallel
#' }
#' @param qcvg (0-100) Default 98 (%). Threshold for query coverage.
#' @param pctidt (0-100) Default 98 (%). Threshold for percentage identity.
#' @param max_out (>0) Default 500. Maximum number of alignments
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
                      org = "bac",
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
                      batches = 1,
                      chunks = 1,
                      qcvg = 98,
                      pctidt = 98,
                      max_out = 500) {

  ##### Input checks -----
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
    alt_path,
    org
  )

  ##### BLAST prep and run -----
  ASVs <- dada2::getSequences(seq_in) %>%
    base::as.data.frame() %>%
    purrr::set_names("ASVs") %>%
    dplyr::mutate(ID = 1:nrow(.)) %>%
    dplyr::select(
      .data$ID,
      .data$ASVs
    )

  if (run_blst & owrt) {
    ASV_splt <- batch(
      ASVs,
      batches
    )

    print(
      paste(
        nrow(ASVs),
        " sequences will be aligned - batches: ",
        length(ASV_splt),
        ", chunks: ",
        chunks,
        sep = ""
      )
    )

    prgrs_bar <- new_bar(
      length(ASV_splt),
      "aligning"
    )
    assign("prgrs_bar",
      prgrs_bar,
      envir = .GlobalEnv
    )

    blst_run <- ASV_splt %>%
      purrr::iwalk(~ blast_run(
        .x,
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
        chunks
      ))
  }


  ##### Alternative annotations (Silva + RDP + UNITE) -----
  if (alt_annot) {
    list_dbs <- list.files(alt_path)
    set.seed(10)

    if (org == "bac") {
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
            var = "ASVs"
          ) %>%
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

    if (org == "fungi") {
      if (owrt_unite) {
        print("assigning taxonomies from UNITE")
        annot_unite <- dada2::assignTaxonomy(
          seq_in,
          here::here(
            alt_path,
            stringr::str_subset(
              list_dbs,
              "sh_general.*fasta"
            )
          ),
          minBoot = 80,
          multithread = T,
          tryRC = T
        ) %>%
          data.frame() %>%
          dplyr::rename_with(tolower) %>%
          dplyr::mutate(kingdom = "Eukaryota") %>%
          tibble::rownames_to_column(
            var = "ASVs"
          ) %>%
          dplyr::mutate(
            dplyr::across(
              tidyselect::everything(),
              .fns = ~ stringr::str_remove(.x, ".*__")
            ),
            dplyr::across(
              tidyselect::everything(),
              .fns = ~ ifelse(
                stringr::str_detect(.x, "_fam_.*"),
                paste(
                  "unknown",
                  stringr::str_remove(.x, "_fam_.*")
                ),
                .x
              )
            ),
            family = ifelse(
              stringr::str_detect(
                family,
                "unknown Saccharomycetales"
              ),
              "Debaryomycetaceae",
              family
            )
          )

        fst::write_fst(annot_unite,
          path = here::here(
            tab_path,
            "UNITE_annot.fst"
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
  }


  ##### Process alignment results -----
  align_path <- here::here(
    tab_path,
    paste("Alignment_",
      tab_out,
      ".tsv",
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
          ".tsv",
          sep = ""
        )
      ))
    ) %>%
      magrittr::set_colnames(c(
        "ID",
        "AccID",
        "TaxID",
        "Species",
        "bitscore",
        "Evalue",
        "qcovs",
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
