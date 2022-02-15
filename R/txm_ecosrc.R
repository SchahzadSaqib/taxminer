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
  "MeshHeadingList", 
  "ArticleTitle", 
  "Abstract"
))

#' Text mining and filtration
#'
#' Extract information from NCBI nucleotide and PubMed databases, attaching ecosystem specificity to each
#' accession ID. Different combinations of word banks are used to scan through this data and apply the filtration
#' criteria. The accession IDs are split into batches of 200, and \link[rentrez]{rentrez} is used to communicate with the NCBI databases.
#'
#'@param hit_tbl (Required) Default NULL. Output table obtained from \link[taxminer]{txm_align}. Alternatively,
#'                    A data.frame with at least 3 columns
#' \itemize{
#'  \item ID: sequence/hit ID numbers.
#'  \item AccID: Accession IDs.
#'}
#' @param filt_host (Optional) Default NA. Filter annotations by host
#' @param filt_site (Optional) Default NA. Filter annotations by body site or environment.
#' @param filt_negt (Optional) Default NA. Disregard annotations that contain these terms.
#' @param do_filt (Logical) Default TRUE. Perform filtration using the word banks. If FALSE the output will
#'                   contain all accession IDs and the extracted information associated to them.
#' @param asbd_tbl (Optional) Default NA. Specify the name of the pre-compiled database
#'                     present within the directory.
#' @param asgn_tbl (Optional) Default Dataset + system date. Name of a new compiled database
#'                            to be assigned
#' @param savedata (Logical) Default TRUE. Should a compiled database be saved to directory?
#' @param sys.break (integer) Default 1. Amount of time, in seconds, that the
#'                  system is paused between iterations. This is handy with larger
#'                  queries to reduce to load on the NCBI servers.         
#'
#' @export
txm_ecosrc <- function(
  hit_tbl,
  filt_host = NA,
  filt_site = NA,
  filt_negt = NA,
  do_filt = T,
  asbd_tbl = NA,
  asgn_tbl = paste("Dataset_", 
                   Sys.Date(), 
                   ".fst", sep = ""),
  savedata = T,
  sys.break = 1) {
  
  check_ecosrc(hit_tbl, 
               filt_host, 
               filt_site, 
               filt_negt)
  
  ##### Input preparation -----
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
                sum(nrow(asbd_tbl_sub), 
                    nrow(ids_src)),
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
                                 length.out = nrow(.))) %>%
      dplyr::group_by(.data$chunks) %>%
      tidyr::nest() %>%
      dplyr::pull()
    
    if (!dir.exists(
      here::here("temp_files"))) {
      dir.create(
        here::here("temp_files"), 
        recursive = T)
    }
    
    ##### Accession ID retrieval -----
    if (!file.exists(
      here::here("temp_files", 
                 "AccID_temp.fst"))) {
      acc_step <- 1
      doc_sum <- data.frame()
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
      print("Retrieving Accession ID data")
      while (acc_step <= base::ceiling(splits)) {
        chunk_ID <- try({
          ids_post <- rentrez::entrez_post(
            db = "nuccore",
            id = paste(
              purrr::flatten(ids_src[acc_step])$AccID,
              collapse = ","))
          ids_summ <- rentrez::entrez_summary(
            db = "nucleotide",
            web_history = ids_post,
            version = "2.0",
            retmode = "xml")
          ids_extr <- t(rentrez::extract_from_esummary(
            ids_summ,
            elements = c("AccessionVersion", 
                         "Organism",
                         "Strain", 
                         "Title", 
                         "TaxId",
                         "SubType", 
                         "SubName"))) %>%
            base::as.data.frame() %>% 
            dplyr::mutate(
              dplyr::across(
                tidyselect::vars_select_helpers$where(is.list),
                .fns = ~as.character(.x)))
        })
        if(!class(chunk_ID) == "try-error") {
          doc_sum <- doc_sum %>%
            dplyr::bind_rows(ids_extr)
          acc_step <- acc_step + 1
          pb_assign$tick()
          if (acc_step > base::ceiling(splits)) {
            message("Done!")
          }
          Sys.sleep(sys.break)
        } else {
          print("trying again")
          acc_step <- acc_step
          Sys.sleep(sys.break + 10)
        }
      }
      
      fst::write_fst(
        doc_sum, 
        here::here("temp_files",
                   "AccID_temp.fst"), 
        100)
    } else {
      doc_sum <- fst::read_fst(
        here::here("temp_files", 
                   "AccID_temp.fst"))
    }
    
    doc_sum <- doc_sum %>%
      dplyr::mutate(meta = Clean_list(doc_sum)) %>%
      dplyr::select(-c(.data$SubType, 
                       .data$SubName))
    
    ##### PMIDs retrieval -----
    if (file.exists(
      here::here("temp_files", 
                 "pmids_temp.fst"))) {
      pmids <- fst::read_fst(
        here::here("temp_files", 
                   "pmids_temp.fst"))
      pmid_step <- readr::read_rds(
        here::here("temp_files", 
                   "pmids_step.rds"))
      pb_assign <- readr::read_rds(
        here::here("temp_files", 
                   "pb_assign.rds"))
    } else {
      pmid_step <- 1
      pmids <- data.frame()
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
    }
    if (pmid_step > base::ceiling(splits)) {
      message("Done!")
    } else {
      print("Retrieving PubMed links")
      while (pmid_step <= base::ceiling(splits)) {
        chunk_link <- try({
          ids_elink <- rentrez::entrez_link(
            db = "pubmed", 
            dbfrom = "nuccore",
            id = purrr::flatten(ids_src[pmid_step])$AccID,
            cmd = "neighbor", 
            by_id = T,
            rettype = "native", 
            idtype = "acc"
          )
          
          if (length(ids_elink) > 0) {
            names(ids_elink) <- purrr::flatten(
              ids_src[pmid_step])$AccID
            ids_elink <- purrr::map(
              ids_elink,
              .f = function(x) 
                x$links[["nuccore_pubmed"]]) %>%
              purrr::compact() %>%
              purrr::map(.f = ~ paste(.x, 
                                      collapse = ",")) %>%
              purrr::map_df(.f = ~.x)
            
            if (nrow(ids_elink) > 0) {
              ids_el <- ids_elink %>%
                tidyr::pivot_longer(., cols = tidyselect::everything(),
                                    names_to = "AccID",
                                    values_to = "pmids") %>%
                dplyr::group_by(.data$AccID) %>%
                dplyr::distinct(.data$pmids, 
                                .keep_all = T)
            }
          }
        })
        if(!class(chunk_link)[1] == "try-error") {
          if (!is.null(chunk_link)) {
            pmids <- pmids %>%
              dplyr::bind_rows(chunk_link)
            
            fst::write_fst(
              pmids,
              here::here("temp_files", 
                         "pmids_temp.fst"), 
              100)
          }
          pmid_step <- pmid_step + 1
          readr::write_rds(
            pmid_step,
            here::here("temp_files", 
                       "pmids_step.rds"))
          pb_assign$tick()
          
          readr::write_rds(
            pb_assign,
            here::here("temp_files", 
                       "pb_assign.rds"))
          
          if (pmid_step > base::ceiling(splits)) {
            message("Done!")
          }
          Sys.sleep(sys.break)
        } else {
          print("reconnecting")
          pmid_step <- pmid_step
          Sys.sleep(sys.break + 10)
        }
      }
    }
    
    ##### PubMed data retrieval -----
    if (nrow(pmids) > 0) {
      pmids_src <- pmids %>%
        dplyr::distinct(pmids)
      
      splits <- round(nrow(pmids_src) / 200)
      if (splits < 1) splits <- 1
      
      pmids_src <- pmids_src %>%
        dplyr::mutate(chunks = rep(0:splits, 
                                   each = 200, 
                                   length.out = nrow(.))) %>%
        dplyr::group_by(.data$chunks) %>% 
        tidyr::nest() %>%
        dplyr::pull()
      
      pmid_step <- 1
      pmids_data <- data.frame()
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
      print("Retrieving PubMed data")
      while (pmid_step <= base::ceiling(splits)) {
        chunk_pub <- try({
          pmids_post <- rentrez::entrez_post(
            db = "pubmed",
            id = paste(purrr::flatten(pmids_src[pmid_step])$pmids,
                       collapse = ","
            )
          )
          pmids_fetch <- rentrez::entrez_fetch(
            db = "pubmed", 
            web_history = pmids_post,
            rettype = "native", 
            retmode = "xml",
            version = "2.0",
            parsed = T
          )
          xml_step <- 1
          while (xml_step <= nrow(XML::xmlToDataFrame(pmids_fetch))) {
            med_sub <- XML::xmlRoot(pmids_fetch) %>%
              .[[xml_step]] %>%
              .[["MedlineCitation"]]
            
            pmids_mesh <- med_sub %>%
              XML::xmlChildren() %>%
              purrr::set_names(
                base::as.character(
                  base::make.unique(names(.)))) %>%
              purrr::keep(.x = ., 
                          names(.) %in% c("PMID", 
                                          "MeshHeadingList")) %>%
              purrr::map_dfr(.x = ., 
                             .f = XML::xmlValue) %>%
              base::as.data.frame() %>%
              dplyr::rename("pmids" = .data$PMID)
            
            titl_abst <- med_sub %>%
              .[["Article"]] %>%
              XML::xmlChildren() %>%
              purrr::set_names(
                base::as.character(
                  base::make.unique(names(.)))) %>%
              purrr::keep(.x = ., 
                          names(.) %in% c("Journal", 
                                          "ArticleTitle", 
                                          "Abstract")) %>%
              purrr::map_dfr(.x = ., 
                             .f = XML::xmlValue) %>%
              base::as.data.frame() %>%
              dplyr::mutate(
                Abstract = ifelse("Abstract" %in% names(.), 
                                  Abstract, 
                                  NA), 
                Journal = ifelse("Journal" %in% names(.), 
                                 Journal, 
                                 NA), 
                ArticleTitle = ifelse("ArticleTitle" %in% names(.),
                                      ArticleTitle, 
                                      NA), 
                Article = paste(
                  "Journal: ", 
                  Journal, 
                  " - ArticleTitle: ", 
                  ArticleTitle, 
                  " - Abstract: ",
                  Abstract, 
                  sep = "")
              ) %>%
              dplyr::bind_cols(pmids_mesh) %>%
              dplyr::select(
                dplyr::any_of(c("pmids", 
                                "Article", 
                                "MeshHeadingList"))) %>%
              dplyr::mutate(MeshHeadingList = ifelse(
                "MeshHeadingList" %in% names(.), 
                MeshHeadingList, 
                NA),
                Article = ifelse(
                  "Article" %in% names(.), 
                  Article, 
                  NA))
            
            pmids_data <- pmids_data %>%
              dplyr::bind_rows(titl_abst) 
            
            xml_step <- xml_step + 1
          }
        })
        if(!class(chunk_pub) == "try-error") {
          pmid_step <- pmid_step + 1
          pb_assign$tick()
          if (pmid_step > base::ceiling(splits)) {
            message("Done!")
          }
          Sys.sleep(sys.break)
        } else {
          print("reconnecting")
          pmid_step <- pmid_step
          Sys.sleep(sys.break + 10)
        }
      }
      
      pmids_data <- pmids_data %>%
        dplyr::left_join(pmids, by = "pmids") %>%
        dplyr::select(.data$AccID, 
                      .data$pmids, 
                      .data$Article,
                      .data$MeshHeadingList)
      
      
      doc_sum <- doc_sum %>%
        dplyr::rename("AccID" = AccessionVersion) %>%
        dplyr::left_join(pmids_data, 
                         by = "AccID")
    } else {
      doc_sum <- doc_sum %>%
        dplyr::rename("AccID" = AccessionVersion) %>%
        dplyr::mutate(MeshHeadingList = NA) %>%
        dplyr::mutate(Article = NA)
    }
    
    if (savedata) {
      if (!is.null(doc_sum)) {
        if (!is.na(asbd_tbl)) {
          if (!is.null(ids_src)) {
            print("Writing new data to Pre-compiled dataset")
            asbd_full <- fst::read_fst(asbd_tbl) %>%
              dplyr::bind_rows(doc_sum) %>%
              dplyr::distinct(.data$AccID, 
                              .keep_all = T) %>%
            fst::write_fst(asbd_tbl, 
                           100)
          }
        } else {
          print("Creating dataset")
          fst::write_fst(doc_sum, 
                         asgn_tbl, 
                         100)
        }
      }
    }
    if (!is.na(asbd_tbl)) {
      doc_sum <- doc_sum %>%
        dplyr::bind_rows(asbd_tbl_sub)
    }
  } else {
    doc_sum <- asbd_tbl_sub
  }
  
  unlink("temp_files", 
         recursive = T)
  
  doc_sum <- doc_sum %>%
    {
      if ("TaxID" %in% names(hit_tbl))
        dplyr::select(., -.data$TaxId)
      else
        dplyr::rename(., "TaxID" = .data$TaxId)
    } %>%
    dplyr::inner_join(hit_tbl) %>%
    {
      if ("Species" %in% names(.))
        dplyr::select(., 
                      .data$ID, 
                      .data$Species, 
                      -.data$Organism, 
                      tidyselect::everything())
      else
        dplyr::rename(., 
                      "Species" = Organism)
    } %>%
    dplyr::arrange(.data$ID) %>%
    dplyr::select(.data$ID, 
                  .data$AccID, 
                  .data$Species, 
                  tidyselect::everything())
  
  
  
  ##### filter preparation ----  
  if (do_filt) {
    
    if (!is.na(filt_host[1])) {
      host_wbank <- purrr::map(
        filt_host, 
        .f = collapse_multi) %>%
        purrr::set_names(filt_host)
    } else {
      host_wbank <- NA
    }
    
    if (!is.na(filt_site[1])) {
      site_wbank <- purrr::map(
        filt_site, 
        .f = collapse_multi) %>%
        purrr::set_names(filt_site)
    } else {
      site_wbank <- NA
    }
    
    if (!is.na(filt_negt[1])) {
      negt_wbank <- purrr::map(
        filt_negt, 
        .f = collapse_multi) %>%
        purrr::set_names(filt_negt)
    } else {
      negt_wbank <- NA
    }
    
    df_extrt <- doc_sum %>%
      dplyr::mutate(
        dplyr::across(
          .cols = tidyselect::vars_select_helpers$where(is.factor), 
                           .fns = ~as.character(.x))) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(host = trimws(
        stringr::str_extract(
          .data$meta,
          pattern = stringr::regex("(?<=host:).*?(?=\\s-)")
        ))) %>%
      dplyr::mutate(Isolation_source = trimws(
        stringr::str_extract(
          .data$meta,
          pattern = stringr::regex(
            paste("(?<=isolation_source:).*?(?=\\s-)",
                  "(?<=isolation_source:).*?(?=$)", 
                  sep = "|"))))) %>%
      as.data.frame()
    
    df_extrt <- df_extrt %>%
      tidyr::unite("pooled_data", 
                   c(.data$MeshHeadingList, 
                     .data$Article, 
                     .data$meta),
                   na.rm = T)
    
    df_filt <- df_extrt
    
    ##### filtration -----
    if (!is.na(host_wbank[1])) {
      
      host_keep <- c()
      host_step <- 1
      while (host_step <= length(host_wbank)) {
        print(paste("Searching through ", 
                    names(host_wbank[host_step]), 
                    " word banks", sep = ""))
        host_filt <- df_filt %>%
          dplyr::mutate(
            dplyr::across(
              .cols = tidyselect::vars_select_helpers$where(is.factor), 
              .fns = ~as.character(.x))) %>%
          dplyr::filter(
            stringr::str_detect(
              pooled_data, 
              stringr::regex(
                as.character(
                  unlist(host_wbank[host_step])),
                ignore_case = T
              ))) %>%
          dplyr::filter(
            is.na(host) | stringr::str_detect(
              host, 
              pattern = stringr::regex(
                as.character(
                  unlist(host_wbank[host_step])),
                ignore_case = T
              )
            ))
        
        host_keep <- host_keep %>%
          dplyr::bind_rows(host_filt)
        
        df_filt <- df_filt %>%
          dplyr::filter(!.data$ID %in% host_keep$ID) %>%
          dplyr::group_by(.data$ID, .data$AccID) %>%
          dplyr::distinct(.data$AccID, .keep_all = T) %>%
          dplyr::ungroup()
        
        host_step <- host_step + 1
      }
      host_keep <- host_keep %>%
        dplyr::arrange(.data$ID) %>%
        dplyr::select(.data$ID, 
                      .data$AccID, 
                      tidyselect::everything())
      
      df_filt <- host_keep %>%
        dplyr::group_by(.data$ID, .data$AccID) %>%
        dplyr::distinct(.data$AccID, .keep_all = T) %>%
        dplyr::ungroup()
      
      dscrd_host <- df_extrt %>%
        dplyr::filter(!AccID %in% df_filt$AccID)
      base::assign("dscrd_host", 
                   dscrd_host, 
                   .GlobalEnv)
    }
    
    if (!is.na(negt_wbank[1])) {
      negt_keep <- c()
      negt_step <- 1
      while (negt_step <= length(negt_wbank)) {
        print(paste("Searching through ", 
                    names(negt_wbank[negt_step]), 
                    " word banks", 
                    sep = ""))
        negt_filt <- df_filt %>%
          dplyr::mutate(
            dplyr::across(
              .cols = tidyselect::vars_select_helpers$where(is.factor), 
              .fns = ~as.character(.x))) %>%
          dplyr::filter(
            stringr::str_detect(
              pooled_data, 
              stringr::regex(
                as.character(
                  unlist(negt_wbank[negt_step])),
                ignore_case = T), negate = T) |
              stringr::str_detect(pooled_data, "ISHAM") |
              stringr::str_detect(host, 
                                  as.character(
                                    unlist(host_wbank)))
          )
        
        negt_keep <- negt_keep %>%
          dplyr::bind_rows(negt_filt)
        
        df_filt <- df_filt %>%
          dplyr::filter(!.data$ID %in% negt_keep$ID) %>%
          dplyr::group_by(.data$ID, .data$AccID) %>%
          dplyr::distinct(.data$AccID, .keep_all = T) %>%
          dplyr::ungroup()
        
        negt_step <- negt_step + 1
      }
      negt_keep <- negt_keep %>%
        dplyr::arrange(.data$ID) %>%
        dplyr::select(.data$ID, 
                      .data$AccID, 
                      tidyselect::everything())
      
      df_filt <- negt_keep %>%
        dplyr::group_by(.data$ID, .data$AccID) %>%
        dplyr::distinct(.data$AccID, .keep_all = T) %>%
        dplyr::ungroup()
      
      dscrd_negt <- df_extrt %>%
        {
          if (base::exists("dscrd_host") == T)
            dplyr::filter(., !.$AccID %in% dscrd_host$AccID)
          else .
        } %>%
        dplyr::filter(!AccID %in% df_filt$AccID)
      base::assign("dscrd_negt", dscrd_negt, .GlobalEnv)
    }
    
    if (!is.na(site_wbank[1])) {
      
      site_keep <- c()
      site_step <- 1
      while (site_step <= length(site_wbank)) {
        print(paste("Searching through ", 
                    names(site_wbank[site_step]), 
                    " word banks", 
                    sep = ""))
        site_filt <- df_filt %>%
          dplyr::mutate(
            dplyr::across(
              .cols = tidyselect::vars_select_helpers$where(is.factor), 
              .fns = ~as.character(.x))) %>%
          dplyr::filter(
            stringr::str_detect(
              pooled_data, 
              stringr::regex(
                as.character(
                  unlist(site_wbank[site_step])),
                ignore_case = T
              ))) %>%
          dplyr::filter(
            is.na(Isolation_source) | stringr::str_detect(
              Isolation_source, 
              pattern = stringr::regex(
                as.character(
                  unlist(site_wbank[site_step])),
                ignore_case = T
              )
            ))
        
        site_keep <- site_keep %>%
          dplyr::bind_rows(site_filt)
        
        df_filt <- df_filt %>%
          dplyr::filter(!.data$ID %in% site_keep$ID) %>%
          dplyr::group_by(.data$ID, .data$AccID) %>%
          dplyr::distinct(.data$AccID, .keep_all = T) %>%
          dplyr::ungroup()
        
        site_step <- site_step + 1
      }
      site_keep <- site_keep %>%
        dplyr::arrange(.data$ID) %>%
        dplyr::select(.data$ID, 
                      .data$AccID, 
                      tidyselect::everything())
      
      discarded_site <- df_extrt %>%
        {
          if (base::exists("dscrd_host") == T)
            dplyr::filter(., !.$AccID %in% dscrd_host$AccID)
          else .
        } %>%
        {
          if (base::exists("dscrd_negt") == T)
            dplyr::filter(., !.$AccID %in% dscrd_negt$AccID)
          else .
        } %>%
        dplyr::filter(!AccID %in% site_keep$AccID) %>%
        dplyr::filter(!ID %in% site_keep$ID)
      
      base::assign("discarded_site", 
                   discarded_site, 
                   .GlobalEnv)
      
      df_filt <- site_keep %>%
        dplyr::group_by(.data$ID, .data$AccID) %>%
        dplyr::distinct(.data$AccID, .keep_all = T) %>%
        dplyr::ungroup()
    }
  } else {
    doc_sum <- doc_sum
  }
}
