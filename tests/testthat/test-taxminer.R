test_data <- data.frame(1:50) %>%
  dplyr::mutate(AccID = c(
    'CP046311.1','CP049223.1','CP049781.1','MK713563.1','KU726641.1','MH898663.1',
    'KC335149.1','CP046311.1','MN559429.1','MH898666.1','CP049223.1','MH898659.1',
    'CP049223.1','CP033426.1','KU726632.1','CP049223.1','KF738669.1','CP049223.1',
    'CP049226.1','KF280299.1','CP011280.1','JX104004.1','KU726642.1','KF007179.1',
    'CP049225.1','MH898660.1','KP996675.1','CP036376.1','KU726679.1','KP192298.1',
    'KP192306.1','KP192307.1','JX104020.1','KJ868804.1','MH898659.1','JX104009.1',
    'KP192303.1','KF007179.1','CP011280.1','CP011280.1','MN559429.1','KC999390.1',
    'JX104011.1','CP049223.1','CP049223.1','MK713563.1','AY959126.1','CP046311.1',
    'KC311734.1','AY078425.1'
  )) %>%
  purrr::set_names("ID", "AccID")

Search_test <- txm_ecosrc(
  input_table = test_data,
  filter_host = "human",
  filter_site = c("vagina", "FRS", "gut+oral+skin+clinical"),
  filter_negate = "non_human", savedata = F)

Lineage_test <- txm_lineage(Search_test)

unlink(paste("Dataset_", Sys.Date(), ".rds", sep = ""))
unlink("_snaps")

testthat::expect_length(Search_test, 11)
testthat::expect_equal(nrow(Search_test), 46)
testthat::expect_length(Lineage_test, 17)
