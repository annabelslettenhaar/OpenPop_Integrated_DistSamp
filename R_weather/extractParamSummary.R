
extractParamSummary <- function(modelOutput, parameterName, byArea = FALSE) {
  
  # 1. Convert to matrix
  samps <- as.matrix(modelOutput)
  
  # 2. Extract relevant columns
  param_cols <- grep(paste0("^", parameterName), colnames(samps), value = TRUE)
  
  # 3. Convert to data frame
  samps_df <- as.data.frame(samps[, param_cols])
  samps_df$iteration <- seq_len(nrow(samps_df))
  
  # 4. Pivot longer
  param_long <- samps_df %>%
    pivot_longer(
      cols = -iteration,
      names_to = "param",
      values_to = "value"
    ) %>%
    mutate(param = as.character(param))
  
  # 5. Extract indices
  if (byArea) {
    # Use gsub to extract area and year from two-index parameter
    param_long <- param_long %>%
      mutate(
        area = as.integer(gsub(paste0(parameterName, "\\[(\\d+),\\s*(\\d+)\\]"), "\\1", param)),
        year = as.integer(gsub(paste0(parameterName, "\\[(\\d+),\\s*(\\d+)\\]"), "\\2", param))
      )
    
    # 6. Summarize by area and year
    summary_df <- param_long %>%
      group_by(area, year) %>%
      summarise(
        mean = mean(value),
        sd = sd(value),
        .groups = "drop"
      )
    
  } else {
    # Extract year from single-index parameter
    param_long <- param_long %>%
      mutate(
        year = as.integer(str_extract(param, "\\d+(?=\\])"))
      )
    
    # 6. Summarize by year
    summary_df <- param_long %>%
      group_by(year) %>%
      summarise(
        mean = mean(value),
        sd = sd(value),
        .groups = "drop"
      )
  }
  
  return(summary_df)
}
