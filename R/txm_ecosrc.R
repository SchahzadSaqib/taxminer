utils::globalVariables(c(
  "AccessionVersion", "TaxId", "pooled_data", "host", "Isolation_source", "ID", ".",
  "Journal", "ArticleTitle", "Abstract"
))

#' Text mining and filtration
#'
#' Extract information from NCBI nucleotide and PubMed databases, attaching ecosystem specificity to each
#' accession ID. Different combinations of word banks are used to scan through this data and apply the filtration
#' criteria. The accession IDs are split into batches of 200, and \link[rentrez]{rentrez} is used to communicate with the NCBI databases.
#'
#'@param input_table (Required) Default NULL. Output table obtained from \link[taxminer]{txm_align}. Alternatively,
#'                    A data.frame with at least 3 columns
#' \itemize{
#'  \item ID: sequence/hit ID numbers.
#'  \item AccID: Accession IDs.
#'  \item TaxID': taxonomic IDs.
#'}
#' @param filter_host (Optional) Default NA. Filter annotations by host
#' @param filter_site (Optional) Default NA. Filter annotations by body site or environment.
#' @param filter_negate (Optional) Default NA. Disregard annotations that contain these terms.
#' @param do_filter (Logical) Default TRUE. Perform filtration using the word banks. If FALSE the output will
#'                   contain all accession IDs and the extracted information associated to them.
#' @param Precomp_tbl (Optional) Default NA. Specify the name of the pre-compiled database
#'                     present within the directory.
#' @param Precomp_tbl_assign (Optional) Default Dataset + system date. Name of a new compiled database
#'                            to be assigned
#' @param savedata (Logical) Default TRUE. Should a compiled database be saved to directory?
#'
#' @export
txm_ecosrc <- function(
  input_table,
  filter_host = NA,
  filter_site = NA,
  filter_negate = NA,
  do_filter = T,
  Precomp_tbl = NA,
  Precomp_tbl_assign = paste("Dataset_", Sys.Date(), ".rds", sep = ""),
  savedata = T) {

  ##### Check whether the input contains data #####
  if (is.na(input_table[[1]][[1]])) {
    stop("Empty Data.frame -  Aborting Relevance filtration")
  }

  filt_split <- unlist(stringr::str_split(c(filter_host, filter_site, filter_negate), pattern = "\\+"))
  filt_split <- filt_split[!is.na(filt_split)]

  i <- 1
  while(i <= length(filt_split)) {
    if (!filt_split[i] %in% names(Word_banks)) {
      stop(paste("Invalid filtration entry ", filt_split[i], sep = ""))
    } else {
      i <- i + 1
    }
  }


  ##### Extract distinct IDs from the input table #####
  AccIDs_to_src <- input_table %>%
    dplyr::distinct(.data$AccID) %>%
    dplyr::filter(str_detect(AccID, "\\."))

  ##### Read in pre-compiled database #####
  if (!is.na(Precomp_tbl)) {
    print("Reading in dataset and searching for existing Accession ids")
    Precomp_tbl_sub <- readRDS(Precomp_tbl) %>%
      dplyr::inner_join(AccIDs_to_src, by = "AccID")

    AccIDs_to_src <- AccIDs_to_src %>%
      dplyr::filter(!.data$AccID %in% Precomp_tbl_sub$AccID) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(.data$AccID)
    print(paste(nrow(AccIDs_to_src), "of",
                sum(nrow(Precomp_tbl_sub), nrow(AccIDs_to_src)),
                "will be searched for",
                sep = " "
    ))
  }


  ##### Check whether there are any IDs left to search for #####
  if (nrow(AccIDs_to_src) > 0) {
    if (is.na(Precomp_tbl)) {
      print(paste("No dataset provided -", nrow(AccIDs_to_src),
                  "accession ids will be searched for",
                  sep = " "
      ))
    }

    ##### Define split ####
    num_groups <- nrow(AccIDs_to_src) / 200
    if (num_groups < 1) num_groups <- 1
    AccIDs_to_src <- AccIDs_to_src %>%
      dplyr::group_by((dplyr::row_number() - 1) %/% (dplyr::n() / num_groups)) %>%
      tidyr::nest() %>%
      dplyr::pull()

    # Create temp files directory
    dir.create(here::here("temp_files"), recursive = T)

    ##### Accession ID retrieval #####
    if (!file.exists(here::here("temp_files", "AccID_temp.rds"))) {
    i <- 1
    DocSum <- data.frame()
    pb_assign <- progress::progress_bar$new(
      format = "  downloading [:bar] :current/:total (:percent) eta: :eta elapsed: :elapsed",
      total = round(num_groups), clear = FALSE, width= 60)
    print("Retrieving Accession ID data")
    while (i <= round(num_groups)) {
      chunk_ID <- try({
        # Post IDs
        AccIDs_post <- rentrez::entrez_post(
          db = "nuccore",
          id = paste(as.character(unlist(AccIDs_to_src[[i]])),
                     collapse = ","))
        # Download summary
        AccIDs_summary <- rentrez::entrez_summary(
          db = "nucleotide",
          web_history = AccIDs_post,
          version = "2.0",
          retmode = "xml")
        # Extract elements
        AccIDs_extract <- t(rentrez::extract_from_esummary(
          AccIDs_summary,
          elements = c("AccessionVersion", "Organism",
                       "Strain", "Title", "TaxId",
                       "SubType", "SubName"))) %>%
          as.data.frame() %>% dplyr::mutate(dplyr::across(tidyselect::vars_select_helpers$where(is.list),
                                                          .fns = ~as.character(.x)))
      })
      if(!class(chunk_ID) == "try-error") {
        # Bind documents
        DocSum <- DocSum %>%
          dplyr::bind_rows(AccIDs_extract)
        i <- i + 1
        pb_assign$tick()
        if (i > round(num_groups)) {
          message("Done!")
        }
        Sys.sleep(1)
      } else {
        print("trying again")
        i <- i
        Sys.sleep(1)
      }
    }


    # write temp file
    readr::write_rds(DocSum, file = here::here("temp_files", "AccID_temp.rds"))
    } else {
      DocSum <- readr::read_rds(here::here("temp_files", "AccID_temp.rds"))
    }

    # Clean final list
    DocSum <- DocSum %>%
      dplyr::mutate(meta = Clean_list(DocSum)) %>%
      dplyr::select(-c(.data$SubType, .data$SubName))



    ###### Elink extraction ######
    if (!file.exists(here::here("temp_files", "PMIDS_temp.rds"))) {
    i <- 1
    PMIDs <- data.frame()
    pb_assign <- progress::progress_bar$new(
      format = "  downloading [:bar] :current/:total (:percent) eta: :eta elapsed: :elapsed",
      total = round(num_groups), clear = FALSE, width= 60)
    print("Retrieving PubMed links")
    while (i <= round(num_groups)) {
      chunk_link <- try({
        # Find links
        AccIDs_elink <- rentrez::entrez_link(
          db = "pubmed", dbfrom = "nuccore",
          id = c(paste(as.character(unlist(AccIDs_to_src[[i]])))),
          cmd = "neighbor", by_id = T,
          rettype = "native", idtype = "acc"
      )
      })
      if(is.list(chunk_link)) {
      # Extract Links
      if (length(AccIDs_elink) > 0) {
        names(AccIDs_elink) <- c(paste(as.character(unlist(AccIDs_to_src[[i]]))))
        AccIDs_elink <- purrr::map(AccIDs_elink,
                                   .f = function(x) x$links[["nuccore_pubmed"]]) %>%
          purrr::compact() %>%
          purrr::map(.f = ~ paste(.x, collapse = ",")) %>%
          purrr::map_df(.f = ~.x) %>%
          {
            if (nrow(.) > 0)
              tidyr::pivot_longer(., cols = tidyselect::everything(),
                                  names_to = "AccID",
                                  values_to = "PMIDs") %>%
              dplyr::group_by(.data$AccID) %>%
              dplyr::distinct(.data$PMIDs, .keep_all = T)
          }
      }
      if (!is.null(AccIDs_elink)) {
        PMIDs <- PMIDs %>%
          dplyr::bind_rows(AccIDs_elink)
      }
        i <- i + 1
        pb_assign$tick()
        if (i > round(num_groups)) {
          message("Done!")
        }
        Sys.sleep(1)
      } else {
        print("trying again")
        i <- i
        Sys.sleep(1)
      }
    }

    # write temp PMIDs
    readr::write_rds(PMIDs, file = here::here("temp_files", "PMIDs_temp.rds"))
    } else {
      PMIDs <- readr::read_rds(here::here("temp_files", "PMIDs_temp.rds"))
    }

    ##### PubMed data #####
    if (nrow(PMIDs) > 0) {
      PMIDs_to_src <- PMIDs %>%
        dplyr::distinct(PMIDs)


      ##### Define split ####
      num_groups <- nrow(PMIDs_to_src) / 200
      if (num_groups < 1) num_groups <- 1

      PMIDs_to_src <- PMIDs_to_src %>%
        dplyr::group_by((dplyr::row_number() - 1) %/% (dplyr::n() / num_groups)) %>%
        tidyr::nest() %>%
        dplyr::pull()

      i <- 1
      PMIDs_data <- data.frame()
      pb_assign <- progress::progress_bar$new(
        format = "  downloading [:bar] :current/:total (:percent) eta: :eta elapsed: :elapsed",
        total = round(num_groups), clear = FALSE, width= 60)
      print("Retrieving PubMed data")
      while (i <= round(num_groups)) {
        chunk_pub <- try({
          PMIDs_post <- rentrez::entrez_post(
            db = "pubmed",
            id = paste(as.character(unlist(PMIDs_to_src[[i]])),
                      collapse = ","
          )
        )
        PMIDs_fetch <- rentrez::entrez_fetch(
          db = "pubmed", web_history = PMIDs_post,
          rettype = "native", retmode = "xml", parsed = T
        )
        j <- 1
        while (j <= nrow(XML::xmlToDataFrame(PMIDs_fetch))) {
          Med_sub <- XML::xmlRoot(PMIDs_fetch) %>%
            .[[j]] %>%
            .[["MedlineCitation"]]

          PMIDs_MESH <- Med_sub %>%
            XML::xmlChildren() %>%
            purrr::set_names(as.character(make.unique(names(.)))) %>%
            purrr::keep(.x = ., names(.) %in% c("PMID", "MeshHeadingList")) %>%
            purrr::map_dfr(.x = ., .f = XML::xmlValue) %>%
            as.data.frame() %>%
            dplyr::rename("PMIDs" = .data$PMID)

          Title_Abstract <- Med_sub %>%
            .[["Article"]] %>%
            XML::xmlChildren() %>%
            purrr::set_names(as.character(make.unique(names(.)))) %>%
            purrr::keep(.x = ., names(.) %in% c("Journal", "ArticleTitle", "Abstract")) %>%
            purrr::map_dfr(.x = ., .f = XML::xmlValue) %>%
            as.data.frame() %>%
            dplyr::mutate(Abstract = ifelse("Abstract" %in% names(.), Abstract, NA)) %>%
            dplyr::mutate(Journal = ifelse("Journal" %in% names(.), Journal, NA)) %>%
            dplyr::mutate(ArticleTitle = ifelse("ArticleTitle" %in% names(.),
                                                ArticleTitle, NA)) %>%
            dplyr::mutate(Article = paste0(
              "Journal: ", Journal, " - ArticleTitle: ", ArticleTitle, " - Abstract: ",
              Abstract, sep = "")
            ) %>%
            dplyr::bind_cols(PMIDs_MESH) %>%
            dplyr::select(dplyr::any_of(c("PMIDs", "Article", "MeshHeadingList")))

          PMIDs_data <- PMIDs_data %>%
            dplyr::bind_rows(Title_Abstract) %>%
            dplyr::mutate(MeshHeadingList = ifelse(
              "MeshHeadingList" %in% names(.), MeshHeadingList, NA)) %>%
            dplyr::mutate(Article = ifelse(
              "Article" %in% names(.), Article, NA))
          j <- j + 1
        }
        })
        if(!class(chunk_pub) == "try-error") {
          i <- i + 1
          pb_assign$tick()
          if (i > round(num_groups)) {
            message("Done!")
          }
          Sys.sleep(1)
        } else {
          print("trying again")
          i <- i
          Sys.sleep(1)
        }
      }

      PMIDs_data <- PMIDs_data %>%
        dplyr::left_join(PMIDs, by = "PMIDs") %>%
        dplyr::select(.data$AccID, .data$PMIDs, .data$Article,
                      .data$MeshHeadingList)


      DocSum <- DocSum %>%
        dplyr::rename("AccID" = AccessionVersion) %>%
        dplyr::left_join(PMIDs_data, by = "AccID")
    } else {
      DocSum <- DocSum %>%
        dplyr::rename("AccID" = AccessionVersion) %>%
        dplyr::mutate(MeshHeadingList = NA) %>%
        dplyr::mutate(Article = NA)
    }

    if (savedata) {
      if (!is.null(DocSum)) {
        if (!is.na(Precomp_tbl)) {
          if (!is.null(AccIDs_to_src)) {
            print("Writing new data to Pre-compiled dataset")
            Precomp_full <- readRDS(Precomp_tbl) %>%
              dplyr::bind_rows(DocSum) %>%
              dplyr::distinct(.data$AccID, .keep_all = T) %>%
              saveRDS(file = Precomp_tbl)
          }
        } else {
          print("Creating dataset")
          saveRDS(DocSum, file = Precomp_tbl_assign)
        }
      }
    }
    if (!is.na(Precomp_tbl)) {
      DocSum <- DocSum %>%
        dplyr::bind_rows(Precomp_tbl_sub)
    }
  } else {
    DocSum <- Precomp_tbl_sub
  }

  # Delete temp files
  unlink("temp_files", recursive = T)

  #  Clean and compile final dataset
  DocSum <- DocSum %>%
    {
      if ("TaxID" %in% names(input_table))
        dplyr::select(., -.data$TaxId)
      else
        dplyr::rename(., "TaxID" = .data$TaxId)
    } %>%
    dplyr::inner_join(input_table) %>%
    {
      if ("Species" %in% names(.))
        dplyr::select(., .data$ID, .data$Species, -.data$Organism, tidyselect::everything())
      else
        dplyr::rename(., "Species" = Organism)
    } %>%
    dplyr::arrange(.data$ID) %>%
    dplyr::select(.data$ID, .data$AccID, .data$Species, tidyselect::everything())




  if (do_filter) {

    ##### Host bank #####

    if (!is.na(filter_host[1])) {
      Host_word_bank <- purrr::map(filter_host, .f = collapse_multi) %>%
        purrr::set_names(filter_host)
    } else {
      Host_word_bank <- NA
    }


    ##### Site bank #####

    if (!is.na(filter_site[1])) {
      Site_word_bank <- purrr::map(filter_site, .f = collapse_multi) %>%
        purrr::set_names(filter_site)
    } else {
      Site_word_bank <- NA
    }


    ##### Negate bank #####

    if (!is.na(filter_negate[1])) {
      Negate_word_bank <- purrr::map(filter_negate, .f = collapse_multi) %>%
        purrr::set_names(filter_negate)
    } else {
      Negate_word_bank <- NA
    }


    ##### Filtration #####

    Data_extrt <- DocSum %>%
      dplyr::mutate_if(is.factor, as.character) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(host = trimws(stringr::str_extract(
        .data$meta,
        pattern = stringr::regex("(?<=host:).*?(?=\\s-)")
      ))) %>%
      dplyr::mutate(Isolation_source = trimws(
        stringr::str_extract(
          .data$meta,
          pattern = stringr::regex("(?<=isolation_source:).*?(?=\\s-)|(?<=isolation_source:).*?(?=$)")
        ))) %>%
      tidyr::unite("pooled_data", c(.data$MeshHeadingList, .data$Article, .data$meta),
                   na.rm = T)

    Data_to_filt <- Data_extrt

    ##### Host filtration #####
    if (!is.na(Host_word_bank[1])) {

      Host_data_tokeep <- c()
      i <- 1
      while (i <= length(Host_word_bank)) {
        print(paste("Searching through ", names(Host_word_bank[i]), " word banks", sep = ""))
        Host_data_filt <- Data_to_filt %>%
          dplyr::mutate_if(is.list, as.character) %>%
          dplyr::filter(stringr::str_detect(
            pooled_data, stringr::regex(as.character(unlist(Host_word_bank[i])),
                                        ignore_case = T
            ))) %>%
          dplyr::filter(is.na(host) | stringr::str_detect(
            host, pattern = stringr::regex(as.character(unlist(Host_word_bank[i])),
                                           ignore_case = T
            )
          ))

        Host_data_tokeep <- Host_data_tokeep %>%
          dplyr::bind_rows(Host_data_filt)

        Data_to_filt <- Data_to_filt %>%
          dplyr::filter(!.data$ID %in% Host_data_tokeep$ID) %>%
          dplyr::group_by(.data$ID, .data$AccID) %>%
          dplyr::distinct(.data$AccID, .keep_all = T) %>%
          dplyr::ungroup()

        i <- i + 1
      }
      Host_data_tokeep <- Host_data_tokeep %>%
        dplyr::arrange(.data$ID) %>%
        dplyr::select(.data$ID, .data$AccID, tidyselect::everything())

      Data_to_filt <- Host_data_tokeep %>%
        dplyr::group_by(.data$ID, .data$AccID) %>%
        dplyr::distinct(.data$AccID, .keep_all = T) %>%
        dplyr::ungroup()

      discarded_host <- Data_extrt %>%
        dplyr::filter(!AccID %in% Data_to_filt$AccID)
      base::assign("discarded_host", discarded_host, .GlobalEnv)
    }


    ##### Negate filtration #####
    if (!is.na(Negate_word_bank[1])) {
      Negate_data_tokeep <- c()
      i <- 1
      while (i <= length(Negate_word_bank)) {
        print(paste("Searching through ", names(Negate_word_bank[i]), " word banks", sep = ""))
        Negate_data_filt <- Data_to_filt %>%
          dplyr::mutate_if(is.list, as.character) %>%
          dplyr::filter(stringr::str_detect(
            pooled_data, stringr::regex(as.character(unlist(Negate_word_bank[i])),
                                        ignore_case = T), negate = T) |
              stringr::str_detect(pooled_data, "ISHAM") |
              stringr::str_detect(host, as.character(unlist(Host_word_bank)))
          )

        Negate_data_tokeep <- Negate_data_tokeep %>%
          dplyr::bind_rows(Negate_data_filt)

        Data_to_filt <- Data_to_filt %>%
          dplyr::filter(!.data$ID %in% Negate_data_tokeep$ID) %>%
          dplyr::group_by(.data$ID, .data$AccID) %>%
          dplyr::distinct(.data$AccID, .keep_all = T) %>%
          dplyr::ungroup()

        i <- i + 1
      }
      Negate_data_tokeep <- Negate_data_tokeep %>%
        dplyr::arrange(.data$ID) %>%
        dplyr::select(.data$ID, .data$AccID, tidyselect::everything())

      Data_to_filt <- Negate_data_tokeep %>%
        dplyr::group_by(.data$ID, .data$AccID) %>%
        dplyr::distinct(.data$AccID, .keep_all = T) %>%
        dplyr::ungroup()

      discarded_negate <- Data_extrt %>%
        {
          if (base::exists("discarded_host") == T)
            dplyr::filter(., !.$AccID %in% discarded_host$AccID)
          else .
        } %>%
        dplyr::filter(!AccID %in% Data_to_filt$AccID)
      base::assign("discarded_negate", discarded_negate, .GlobalEnv)
    }

    ##### Site filtration #####
    if (!is.na(Site_word_bank[1])) {

      Site_data_tokeep <- c()
      i <- 1
      while (i <= length(Site_word_bank)) {
        print(paste("Searching through ", names(Site_word_bank[i]), " word banks", sep = ""))
        Site_data_filt <- Data_to_filt %>%
          dplyr::mutate_if(is.list, as.character) %>%
          dplyr::filter(stringr::str_detect(
            pooled_data, stringr::regex(as.character(unlist(Site_word_bank[i])),
                                        ignore_case = T
            ))) %>%
          dplyr::filter(is.na(Isolation_source) | stringr::str_detect(
            Isolation_source, pattern = stringr::regex(as.character(unlist(Site_word_bank[i])),
                                                       ignore_case = T
            )
          ))

        Site_data_tokeep <- Site_data_tokeep %>%
          dplyr::bind_rows(Site_data_filt)

        Data_to_filt <- Data_to_filt %>%
          dplyr::filter(!.data$ID %in% Site_data_tokeep$ID) %>%
          dplyr::group_by(.data$ID, .data$AccID) %>%
          dplyr::distinct(.data$AccID, .keep_all = T) %>%
          dplyr::ungroup()

        i <- i + 1
      }
      Site_data_tokeep <- Site_data_tokeep %>%
        dplyr::arrange(.data$ID) %>%
        dplyr::select(.data$ID, .data$AccID, tidyselect::everything())

      discarded_site <- Data_extrt %>%
        {
          if (base::exists("discarded_host") == T)
            dplyr::filter(., !.$AccID %in% discarded_host$AccID)
          else .
        } %>%
        {
          if (base::exists("discarded_negate") == T)
            dplyr::filter(., !.$AccID %in% discarded_negate$AccID)
          else .
        } %>%
        dplyr::filter(!AccID %in% Site_data_tokeep$AccID) %>%
        dplyr::filter(!ID %in% Site_data_tokeep$ID)
      base::assign("discarded_site", discarded_site, .GlobalEnv)

      Data_to_filt <- Site_data_tokeep %>%
        dplyr::group_by(.data$ID, .data$AccID) %>%
        dplyr::distinct(.data$AccID, .keep_all = T) %>%
        dplyr::ungroup()
    }
  } else {
    DocSum <- DocSum
  }
}
