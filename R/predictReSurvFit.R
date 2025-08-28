#' Predict IBNR frequency
#'
#' This function predicts the results from the ReSurv fits.
#'
#' @param object \code{ResurvFit} object specifying start time, end time and status.
#' @param newdata \code{IndividualDataPP} object that contains new data to predict.
#' @param grouping_method \code{character}, use probability or exposure approach to group from input to output development factors. Choice between:
#' \itemize{
#' \item{\code{"exposure"}}
#' \item{\code{"probability"}}
#' }
#' Default is \code{"exposure"}.
#' @param check_value \code{numeric}, check hazard value on initial granularity, if above threshold we increase granularity to try and adjust the development factor.
#' @param lower_triangular_output \code{logical}, if set to \code{TRUE} we add the predicted lower triangle in input and output granularity to the \code{predict.ReSurvFit} output.
#' @param groups_encoding_output \code{logical}, if set to \code{TRUE} we add a \code{data.table} containing the groups encoding to the \code{predict.ReSurvFit} output.
#' @param ... Additional arguments to pass to the predict function.
#'
#'
#' @return Predictions for the \code{ReSurvFit} model. It includes
#' \itemize{
#' \item{\code{ReSurvFit}: Fitted \code{ReSurv} model.}
#' \item{\code{long_triangle_format_out}: \code{data.frame}. Predicted development factors and IBNR claim counts for each feature combination in long format.}
#' \itemize{
#' \item{\code{input_granularity}: \code{data.frame}. Predictions for each feature combination in long format for \code{input_time_granularity}.}
#' \itemize{
#' \item{\code{AP_i}: Accident period, \code{input_time_granularity}.}
#' \item{\code{DP_i}: Development period, \code{input_time_granularity}.}
#' \item{\code{f_i}: Predicted development factors, \code{input_time_granularity}.}
#' \item{\code{group_i}: Group code, \code{input_time_granularity}. This associates to each feature combination an identifier.}
#' \item{\code{expected_counts}: Expected counts, \code{input_time_granularity}.}
#' \item{\code{IBNR}: Predicted IBNR claim counts, \code{input_time_granularity}.}
#' }
#' \item{\code{output_granularity}: \code{data.frame}. Predictions for each feature combination in long format for \code{output_time_granularity}.}
#' \itemize{
#' \item{\code{AP_o}: Accident period, \code{output_time_granularity}.}
#' \item{\code{DP_o}: Development period, \code{output_time_granularity}.}
#' \item{\code{f_o}: Predicted development factors, \code{output_time_granularity}.}
#' \item{\code{group_o}: Group code, \code{output_time_granularity}. This associates to each feature combination an identifier.}
#' \item{\code{expected_counts}: Expected counts, \code{output_time_granularity}.}
#' \item{\code{IBNR}: Predicted IBNR claim counts, \code{output_time_granularity}.}
#' }
#' }
#' \item{\code{lower_triangle}: Predicted lower triangle.}
#' \itemize{
#' \item{\code{input_granularity}: \code{data.frame}. Predicted lower triangle for \code{input_time_granularity}.}
#' \item{\code{output_granularity}: \code{data.frame}. Predicted lower triangle for \code{output_time_granularity}.}
#' }
#' \item{\code{predicted_counts}: \code{numeric}. Predicted total frequencies.}
#' \item{\code{grouping_method}: \code{character}. Chosen grouping method.}
#'
#' }
#'
#' @importFrom dplyr bind_rows distinct relocate arrange
#' @export
#' @method predict ReSurvFit
predict.ReSurvFit <- function(object,
                              newdata = NULL,
                              grouping_method = "probability",
                              lower_triangular_output = TRUE,
                              groups_encoding_output =FALSE,
                              check_value = 1.85,
                              ...) {


  if (!is.null(newdata)) {
    pkg.env$check.newdata(newdata = newdata,
                          pastdata = object$IndividualDataPP)

    idata <- newdata
    # hazard_frame <- adjust.predictions(ResurvFit=object,
    #                                    hazard_model=object$hazard_model,
    #                                    idata=idata)


  } else{
    idata <- object$IndividualDataPP


  }



  is_baseline_model <- is.null(c(object$IndividualDataPP$categorical_features, object$IndividualDataPP$continuous_features))

  # Convert to data.table if not already
  hazard_frame <- as.data.table(object$hazard_frame)

  # Drop DP_i
  hazard_frame[, DP_i := NULL]

  # Rename columns
  setnames(hazard_frame, old = c("f_i", "cum_f_i"), new = c("dev_f_i", "cum_dev_f_i"))

  # hazard_frame <- object$hazard_frame %>%
  #   select(-DP_i) %>%
  #   rename(dev_f_i = f_i, cum_dev_f_i = cum_f_i)

  hazard_frame_grouped <- pkg.env$covariate_mapping(
    hazard_frame = hazard_frame,
    categorical_features = idata$categorical_features,
    continuous_features = idata$continuous_features,
    conversion_factor = idata$conversion_factor,
    calendar_period_extrapolation = idata$calendar_period_extrapolation
  )


  # browser()

  if(object$simplifier){

    missing.obsevations <- pkg.env$simplified_fill_data_frame(
      data = idata$full.data,
      continuous_features = idata$continuous_features,
      categorical_features =
        idata$categorical_features,
      years = idata$years,
      input_time_granularity =
        idata$input_time_granularity,
      conversion_factor = idata$conversion_factor
    )

  }else{
  missing.obsevations <- pkg.env$fill_data_frame(
    data = idata$full.data,
    continuous_features = idata$continuous_features,
    categorical_features =
      idata$categorical_features,
    years = idata$years,
    input_time_granularity =
      idata$input_time_granularity,
    conversion_factor = idata$conversion_factor
  )


  }

  # latest_observed computation ----
  ## Retrieve total amount of observed claims as of the evaluation date ----

  data_reserve = bind_rows(idata$training.data, missing.obsevations)

  max_DP_i <- data_reserve[, .(DP_max_rev = min(max(DP_rev_i) - DP_i) + 1), by = AP_i]

  # Step 1: select relevant columns
  cols_to_keep <- unique(c("AP_i", "AP_o", "DP_rev_i", "DP_i", idata$categorical_features, idata$continuous_features, "I"))
  data_reserve2 <- data_reserve[, ..cols_to_keep]  # .. prefix needed for a vector of names

  # Step 2: convert AP_i to numeric
  data_reserve2[, AP_i := as.numeric(AP_i)]

  # Step 3: left join with max_DP_i by AP_i
  data_reserve2 <- max_DP_i[data_reserve2, on = "AP_i"]

  if(is.null(idata$continuous_features)){

    # Define grouping columns
    group_cols <- c(idata$categorical_features, "AP_i", "AP_o", "DP_max_rev")

    # Compute sum of I by group
    observed_so_far <- data_reserve2[, .(latest_I = sum(I)), by = group_cols]

    # Define grouping columns
    group_cols <- c("AP_i", "AP_o", idata$categorical_features, "DP_rev_i", "DP_i")

    # Summarize I by group
    observed_dp_rev_i <- data_reserve2[, .(I = sum(I)), by = group_cols]

    ## add the group dimension

    # Define feature columns
    feature_cols <- c(idata$categorical_features, idata$continuous_features)

    # Create feature.id efficiently
    observed_so_far[, covariate := do.call(paste, c(.SD, sep = "_")), .SDcols = feature_cols]
    observed_dp_rev_i[, covariate := do.call(paste, c(.SD, sep = "_")), .SDcols = feature_cols]

    # Perform left join
    observed_so_far <- hazard_frame_grouped$groups[observed_so_far, on = "covariate"]

    # Select relevant columns
    observed_so_far <- observed_so_far[, .(AP_i, group_i, DP_max_rev, latest_I)]


    # Perform left join
    observed_dp_rev_i <- hazard_frame_grouped$groups[observed_dp_rev_i, on = "covariate"]

    # Select relevant columns
    observed_dp_rev_i <- observed_dp_rev_i[, .(AP_i, group_i, DP_max_rev, latest_I)]



  }else{

    if(length(idata$continuous_features)==1 & "AP_i" %in% idata$continuous_features){

      continuous_features_group<-NULL

    }else{

      continuous_features_group <- idata$continuous_features[!(idata$continuous_features %in% c("AP_i","RP_i") )]

    }

    if(is.null(continuous_features_group)){

      group_cols <- c("AP_i", "AP_o", idata$categorical_features, "DP_max_rev")
      observed_so_far <- data_reserve2[, .(latest_I = sum(I)), by = group_cols]

      # Define grouping columns
      group_cols <- c("AP_i", "AP_o", idata$categorical_features, "DP_rev_i", "DP_i")

      # Summarise
      observed_dp_rev_i <- data_reserve2[, .(I = sum(I)), by = group_cols]


    }else{

      # Define grouping columns
      group_cols <- c("AP_i", "AP_o", idata$categorical_features, continuous_features_group, "DP_max_rev")

      # Summarise
      observed_so_far <- data_reserve2[, .(latest_I = sum(I)), by = group_cols]


      # Define grouping columns
      group_cols <- c("AP_i", "AP_o", idata$categorical_features, continuous_features_group, "DP_rev_i", "DP_i")

      # Summarise
      observed_dp_rev_i <- data_reserve2[, .(I = sum(I)), by = group_cols]




    }



    # Define feature columns
    feature_cols <- c(idata$categorical_features, continuous_features_group)

    # Create feature.id efficiently
    observed_so_far[, covariate := do.call(paste, c(.SD, sep = "_")), .SDcols = feature_cols]
    observed_dp_rev_i[, covariate := do.call(paste, c(.SD, sep = "_")), .SDcols = feature_cols]


    time_features <- idata$continuous_features[idata$continuous_features %in% c("AP_i","RP_i")]


    # Left join
    latest_cumulative <- hazard_frame_grouped$groups[observed_so_far, on = c(time_features, "covariate")]

    # Select required columns
    latest_cumulative <- latest_cumulative[, c(unique("AP_i", time_features), "group_i", "DP_max_rev", "latest_I"), with = FALSE]


    # Step 1: left join
    observed_dp_rev_i <- hazard_frame_grouped$groups[
      observed_dp_rev_i,
      on = c(time_features, "covariate")
    ]

    # Step 2: keep only required columns
    observed_dp_rev_i <- observed_dp_rev_i[,
                                                   c(unique("AP_i", time_features), "group_i", "DP_rev_i", "DP_i", "I"),
                                                   with = FALSE
    ]

    # Step 3: inner join with (time_features, group_i) from hazard_frame_grouped$groups
    observed_dp_rev_i <- observed_dp_rev_i[
      hazard_frame_grouped$groups[, c(time_features, "group_i"), with = FALSE],
      on = c(time_features, "group_i"),
      nomatch = 0   # ensures inner join
    ]

    # latest_cumulative = observed_so_far_out, observed_pr_dp = observed_dp_rev_i_tmp


    }

    # browser()

  max_DP <- max(bind_rows(idata$training.data, missing.obsevations)$DP_rev_o)



  ## Computation of the expected IBNR ----

  grouped_hazard_0 <- latest_cumulative[
    hazard_frame_grouped$hazard_group,
    on = c("group_i", "AP_i")
  ]

  # Step 1: keep DP_max_rev explicitly
  expected_i <- grouped_hazard_0[, .(DP_rev_i, AP_i, group_i, S_i, S_i_lag,
                                   DP_max_rev, latest_I, DP_max_rev_keep = DP_max_rev)]


  # Step 2: hazard_tmp as before
  hazard_tmp <- hazard_frame_grouped$hazard_group[,
                                                  .(DP_rev_i_key = DP_rev_i + 1, AP_i, group_i, S_ultimate_i = S_i)]

  # Step 3: left join (DP_max_rev used only for matching)
  expected_i <- hazard_tmp[expected_i, on = .(DP_rev_i_key = DP_max_rev, AP_i, group_i)]

  # Step 4: fcase now can use DP_max_rev_keep
  ## in this case you only go with the probability approach
  expected_i[,U:=1]
  ## alternative is to use the exposure approach:
  # expected_i[, U := fcase(
  #   S_i_lag == 1, latest_I,
  #   DP_max_rev_keep == min(hazard_frame_grouped$hazard_group$DP_rev_i), latest_I,
  #   S_ultimate_i == 0, 0,
  #   AP_i != 1, latest_I / S_ultimate_i,
  #   default = latest_I
  # )]


  # Step 5: compute I_expected
  expected_i[, I_expected := U*(S_i_lag - S_i)]

  # Step 6: compute IBNR
  expected_i[, IBNR := fifelse(DP_rev_i < DP_max_rev_keep, I_expected, NA_real_)]

  # Step 7: keep only required columns
  expected_i <- expected_i[, .(AP_i, group_i, DP_rev_i, I_expected, IBNR)]





  browser()

######################################
  df_i <- pkg.env$retrieve_df_i(
    hazard_data_frame = hazard_frame_grouped$hazard_group,
    groups = hazard_frame_grouped$groups,
    is_baseline_model = is_baseline_model
  )
##########################################

  ## create df_i: this part assumes that adjust is false.

  df_i <- unique(hazard_frame_grouped$hazard_group[, .(group_i, DP_rev_i, dev_f_i)])[
    ,
    dcast(.SD, DP_rev_i ~ group_i, value.var = "dev_f_i")
  ][
    ,
    !"DP_rev_i"
  ]


  ###############################################
  if(!adjusted){

    df_i <- hazard_data_frame %>%
      select(group_i, DP_rev_i, dev_f_i) %>%
      distinct() %>%
      reshape2::dcast(DP_rev_i ~group_i, value.var="dev_f_i") %>%
      select(-DP_rev_i)


  }else{
    df_i <- hazard_data_frame %>%
      select(group_i, DP_rev_i, df_i_adjusted) %>%
      distinct() %>%
      reshape2::dcast(DP_rev_i ~group_i, value.var="df_i_adjusted") %>%
      select(-DP_rev_i)
  }

  #We only have 5 columns in the case of AP being included as covariate
  if(ncol(groups) == 5){
    colnames(df_i) <- c(paste0("AP_i_",groups$AP_i,",", groups$covariate ))
  }else{
    colnames(df_i) <- c(groups$covariate )
  }
  #

  if(is_baseline_model){

    df_i <- df_i %>%
      map_df(rev) %>%
      mutate(DP_i=row_number())

    return(df_i)

  }else{
    df_i <- as.data.frame(df_i[1:(nrow(df_i)-1),]) %>%
      map_df(rev) %>%
      mutate(DP_i=row_number())

  }

  return(df_i)
  ##########################################################



  # Creating the input granularity output ----

  hazard_frame_input <- pkg.env$input_hazard_frame(
    hazard_frame = hazard_frame_grouped$hazard_group,
    expected_i = expected_i ,
    categorical_features = idata$categorical_features,
    continuous_features = idata$continuous_features,
    df_i = df_i,
    groups = hazard_frame_grouped$groups,
    is_baseline_model = is_baseline_model
  )



  # browser()

  # In case of output granularity different from input granularity, creating the output -----

  if (idata$conversion_factor != 1) {
    development_periods <- distinct(select(data.frame(idata$training), AP_i, AP_o))

    # Calculate the minimum and maximum development periods for each row in development_periods for each DP_rev_o
    dp_ranges <- t(lapply(1:max_DP, function(DP_rev_o) {
      cbind(
        DP_rev_o,
        development_periods,
        min_dp = with(
          development_periods,
          AP_i + 1 / (idata$conversion_factor) * (DP_rev_o - AP_o)
        ),
        max_dp = with(
          development_periods,
          AP_i - 1 + 1 / (idata$conversion_factor) * (DP_rev_o - AP_o + 1)
        )
      )
    }))

    dp_ranges <- do.call(rbind, dp_ranges)


    #check_input_hazard <- pkg.env$check_input_hazard(hazard_frame_input,
    #                           check_value=check_value)

    #If we exeed the check value, we calculate on ouput granularity, predict on output granularity, and distribute evenly in the relevant input-periods.
    #From here we do simple chain-ladder to calculate new development factor.
    if (#check_input_hazard
      FALSE # I guess we  can remove this part if not used?
      ) {
      development_factor_o <- pkg.env$i_to_o_development_factor(
        hazard_data_frame = hazard_frame_grouped$hazard_group,
        expected_i = expected_i,
        dp_ranges = dp_ranges,
        groups = hazard_frame_grouped$groups,
        observed_pr_dp = latest_observed$observed_pr_dp,
        latest_cumulative = latest_observed$latest_cumulative,
        conversion_factor = idata$conversion_factor,
        grouping_method = "probability",
        min_DP_rev_i = min(hazard_frame_grouped$hazard_group$DP_rev_i)
      )

      if (ncol(hazard_frame_grouped$groups) == 5) {
        colnames(development_factor_o) <- unique(c(
          paste0(
            "AP_o_",
            hazard_frame_grouped$groups$AP_o,
            ",",
            hazard_frame_grouped$groups$covariate
          )
        ))
      }
      else{
        colnames(development_factor_o) <- c(hazard_frame_grouped$groups$covariate)
      }

      df_o <- as.data.frame(development_factor_o[1:(nrow(development_factor_o) -
                                                    1), ]) %>%
        map_df(rev) %>%
        mutate(DP_o = row_number())

      #We only update for relevant periods, hence for example for accident periods, where we have already seen the development, we just put to 1.
      hazard_frame_grouped$hazard_group <- pkg.env$update_hazard_frame(
        hazard_frame_input = hazard_frame_input,
        hazard_frame_grouped = hazard_frame_grouped$hazard_group,
        df_o = df_o,
        latest_observed_i = latest_observed$observed_pr_dp,
        groups = hazard_frame_grouped$groups,
        categorical_features = idata$categorical_features,
        continuous_features = idata$continuous_features,
        conversion_factor = idata$conversion_factor,
        check_value = check_value
      )

      expected_i <- pkg.env$predict_i(
        hazard_data_frame = hazard_frame_grouped$hazard_group,
        latest_cumulative = latest_observed$latest_cumulative,
        grouping_method = "exposure"
      )

      df_i <- pkg.env$retrieve_df_i(
        hazard_data_frame = hazard_frame_grouped$hazard_group,
        groups = hazard_frame_grouped$groups,
        adjusted = T
      )

      hazard_frame_input <- pkg.env$input_hazard_frame(
        hazard_frame = hazard_frame_grouped$hazard_group,
        expected_i = expected_i ,
        categorical_features = idata$categorical_features,
        continuous_features = idata$continuous_features,
        df_i = df_i,
        groups = hazard_frame_grouped$groups,
        adjusted = T
      )



  }

  expected_o <- pkg.env$predict_o(
    expected_i = expected_i,
    groups = hazard_frame_grouped$groups,
    conversion_factor = idata$conversion_factor,
    years = object$IndividualDataPP$years,
    input_time_granularity = object$IndividualDataPP$input_time_granularity
  )


  development_factor_o <- pkg.env$i_to_o_development_factor(
    hazard_data_frame = hazard_frame_grouped$hazard_group,
    expected_i = expected_i,
    dp_ranges = dp_ranges,
    groups = hazard_frame_grouped$groups,
    observed_pr_dp = latest_observed$observed_pr_dp,
    latest_cumulative = latest_observed$latest_cumulative,
    conversion_factor = idata$conversion_factor,
    grouping_method = grouping_method,
    min_DP_rev_i = min(hazard_frame_grouped$hazard_group$DP_rev_i),
    years = object$IndividualDataPP$years,
    input_time_granularity = object$IndividualDataPP$input_time_granularity
  )



  #We only have 5 groups if AP is included as a covariate
  if (ncol(hazard_frame_grouped$groups) == 5) {
    colnames(development_factor_o) <- unique(c(
      paste0(
        "AP_o_",
        hazard_frame_grouped$groups$AP_o,
        ",",
        hazard_frame_grouped$groups$covariate
      )
    ))
  }
  else{
    if (is_baseline_model) {
      colnames(development_factor_o) <- "0"
    } else{
      colnames(development_factor_o) <- c(hazard_frame_grouped$groups$covariate)
    }
  }

  df_o <- as.data.frame(development_factor_o[1:(nrow(development_factor_o) -
                                                  1), ]) %>%
    map_df(rev) %>%
    mutate(DP_o = row_number())


  hazard_frame_output <- pkg.env$output_hazard_frame(
    hazard_frame_input = hazard_frame_input,
    expected_o = expected_o,
    categorical_features = idata$categorical_features,
    continuous_features = idata$continuous_features,
    df_o = df_o,
    groups = hazard_frame_grouped$groups,
    is_baseline_model = is_baseline_model
  )

  # Final output formatting: no actual computations from here on

  max_dp_i = pkg.env$maximum.time(
    object$IndividualDataPP$years,
    object$IndividualDataPP$input_time_granularity
  )

  long_tr_input = hazard_frame_input %>%
    mutate(DP_i = max_dp_i - DP_rev_i + 1) %>%
    select(-expg, -baseline, -hazard, -DP_rev_i) %>%
    relocate(DP_i, .after =  AP_i) %>%
    rename(f_i = df_i, expected_counts = I_expected) %>%
    as.data.frame()

  long_tr_output = hazard_frame_output %>%
    mutate(DP_o = max(DP_rev_o) - DP_rev_o + 1) %>%
    select(-DP_rev_o) %>%
    arrange(DP_o) %>%
    relocate(DP_o, .after =  AP_o) %>%
    rename(f_o = df_o, expected_counts = I_expected) %>%
    as.data.frame()

  rm(list = c(
    "df_o",
    "df_i",
    "hazard_frame_input",
    "hazard_frame_output"
  ))





  out = list(
    ReSurvFit = object,
    # I removed these two (save memory)
    # df_output = as.data.frame(df_o),
    # df_input = as.data.frame(df_i),
    long_triangle_format_out = list(
      input_granularity = long_tr_input,
      output_granularity = long_tr_output
    ),
    predicted_counts = sum(long_tr_input$IBNR, na.rm = T),
    grouping_method = grouping_method
  )


  if(groups_encoding_output){

    d1 <- as.data.table(hazard_frame_grouped$hazard_group[,unique(c(idata$continuous_features,
                                                                    idata$categorical_features,
                                                                    "group_i",
                                                                    "DP_rev_i"))])[!duplicated(group_i)]

    d2 <- as.data.table(hazard_frame_grouped$groups[,c("group_i",
                                                       "group_o")])

    d3 <- merge(d1,d2,by=c("group_i"))

    out[['groups_encoding']] <- lower_triangle

  }


  if (lower_triangular_output) {
    ltr_input <- pkg.env$find_lt_input(long_tr_input, max_dp_i)

    ltr_output <- pkg.env$find_lt_output(
      long_tr_output,
      max(long_tr_output$DP_o),
      cut_point = ceiling(max_dp_i * idata$conversion_factor)
    )

    lower_triangle=list(input_granularity=ltr_input,
                        output_granularity=ltr_output)

    out[['lower_triangle']] <- lower_triangle
  }

  class(out) <- c('ReSurvPredict')

  return(out)
}


max_dp_i = pkg.env$maximum.time(object$IndividualDataPP$years,
                                object$IndividualDataPP$input_time_granularity)

long_tr_input = hazard_frame_input %>%
  mutate(DP_i = max_dp_i - DP_rev_i + 1) %>%
  select(-expg, -baseline, -hazard, -DP_rev_i) %>%
  as.data.frame()


# ltr_input <- pkg.env$find_lt_input(long_tr_input, max_dp_i)

out = list(
  ReSurvFit = object,
  # I removed these two (memory issues)
  # df_input = as.data.frame(df_i),
  long_triangle_format_out = list(input_granularity = long_tr_input),
  # lower_triangle = list(input_granularity = ltr_input),
  predicted_counts = sum(long_tr_input$IBNR, na.rm = T),
  grouping_method = grouping_method
)

if (lower_triangular_output) {
  ltr_input <- pkg.env$find_lt_input(long_tr_input, max_dp_i)

  ltr_output <- pkg.env$find_lt_output(
    long_tr_output,
    max(long_tr_output$DP_o),
    cut_point = ceiling(max_dp_i * idata$conversion_factor)
  )

  lower_triangle=list(input_granularity=ltr_input)

  out[['lower_triangle']] <- lower_triangle
}

class(out) <- c('ReSurvPredict')

return(out)

}
