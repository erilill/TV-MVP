## In this file, we define the methods of the R6 class

TVMVP$set("public", "set", function(...) {
  # for generic set function, the arguments cannot be general because they are private
  # by doing so, it makes no sense if the user input something really nonsense
  # the user can set the variables that have been defined in the class
  # we don't do type checking except for the data
  # but it may pop error message when the variables are misspecified and are being used.
  args <- list(...); arg_names <- names(args)

  # Special handling for 'data'
  if ("data" %in% arg_names) {
    if(is.matrix(args$data) && is.numeric(args$data)) {
      # set the data if it is numeric matrix
      private$data <- args$data
      private$iT <- nrow(args$data)
      private$ip <- ncol(args$data)

      tmp_size <- prettyunits::pretty_bytes(object.size(private$data))
      cli::cli_alert_info("data set {.val {tmp_size}} with {.val {private$iT}} rows and {.val {private$ip}} columns")
    } else {
      # data must be a numeric matrix!
      cli::cli_alert_warning("data must be a numeric matrix")
    }
  }

  # remove 'data' from args so it’s not double-assigned below
  # iT and ip are the variables that are not allowed to be changed manually
  arg_names <- setdiff(arg_names, private$important_fields)

  # nothing will happen if arg_names is empty
  for(name in arg_names){
    if (is.null(private[[name]])){
      cli::cli_alert_info("New variable set: {name}")
    } else {
      cli::cli_alert_info("Variable updated: {name}")
    }
    private[[name]] = args[[name]]
  }

  invisible(self)  # Enables chaining
})


TVMVP$set("public", "set_data", function(data) {
  # for set_data function, the argument data should not be missing
  # this is why I use missing(data) for checking without setting any default value of it
  # note that this is different from initialize function where data can be missing with default value NULL
  if(missing(data)){
    cli::cli_alert_warning("input is empty")
  } else{
    # set the data
    self$set(data=data)
  }

  invisible(self)  # Enables chaining
})

TVMVP$set("public", "get_data", function() {
  if(is.null(private$data)){
    cli::cli_alert_warning("The data is empty")
    return(invisible(NULL))
  }
  return(private$data)
  # make a copy of the data and return it
})


TVMVP$set("public", "set_max_m", function(max_m) {
  # for set_data function, the argument data should not be missing
  if(missing(max_m)){
    cli::cli_alert_warning("input is empty")
  } else{
    # set the data
    self$set(max_m=max_m)
  }

  invisible(self)  # Enables chaining
})

TVMVP$set("public", "set_bandwidth", function(bandwidth) {
  # for set_data function, the argument data should not be missing
  if(missing(bandwidth)){
    cli::cli_alert_warning("input is empty")
  } else{
    # set the data
    self$set(bandwidth=bandwidth)
  }

  invisible(self)  # Enables chaining
})


TVMVP$set("public", "get_optimal_m", function() {
  if(is.null(private$optimal_m))
    cli::cli_alert_warning("run {.code determine_factors(max_m = , bandwidth = )}")
  return(private$optimal_m)
})


TVMVP$set("public", "print", function(...) {
  # print function

  ### default priting

  cli::cli_alert_info("Object of {.strong TVMVP}")

  # Important fields
  important_fields = c("data","iT","ip",
                        "optimal_m","IC_values")

  # print data first
  if(!is.null(private$data)){
    if(is.matrix(private$data) && is.numeric(private$data)){
      tmp_size <- prettyunits::pretty_bytes(object.size(private$data))
      cli::cli_text("data set {.val {tmp_size}} with {.val {private$iT}} rows and {.val {private$ip}} columns")
    } else {
      # data must be a numeric matrix!
      cli::cli_alert_warning("data must be a numeric matrix")
    }
  } else {
    cli::cli_alert_warning("data is empty; use {.code set(data = )} or {.code set_data()}")
  }

  # Other fields
  other_fields <- setdiff(names(private), private$important_fields)
  if (length(other_fields) > 0) {
    for (name in other_fields) {
      if(!is.null(private[[name]]))
        cli::cli_text(" - {.field {name}} = {.val {private[[name]]}}")
    }
  }

  ### print results
  if(!is.null(private$optimal_m))
    cli::cli_text(" - {.field optimal_m} = {.val {private$optimal_m}}")
  if(!is.null(private$IC_values))
    cli::cli_text(" - {.field IC_values} = {.val {private$IC_values}}")
  if(!is.null(private$J_test))
    cli::cli_text(" - {.field J_test} = {.val {private$J_test}}")

  invisible(self)
})

TVMVP$set("public", "silverman", function() {
  # update the private bandwidth by Silverman
  # this can be run automatically in other functions
  if(is.null(private$data)) {
    cli::cli_alert_warning("data is empty")
    return(invisible(self)) # return
  }
  iT = private$iT; ip = private$ip

  private$bandwidth <- (2.35/sqrt(12)) * iT^(-0.2) * ip^(-0.1)

  invisible(self)
})

# if the object has already these arguments, then user need not specify them
# user can change them by inputting a new pair
TVMVP$set("public", "determine_factors", function(max_m=NULL, bandwidth=NULL) {
  if(!is.null(max_m)) private$max_m = max_m
  if(!is.null(bandwidth)) private$bandwidth = bandwidth
  max_m = private$max_m; bandwidth = private$bandwidth

  flag = TRUE
  if(is.null(private$data)) {
    cli::cli_alert_warning("data is empty")
    flag = FALSE
  }
  if(is.null(max_m)) {
    cli::cli_alert_warning("max_m is empty")
    flag = FALSE
  }
  if(is.null(bandwidth)) {
    if(flag){
      # we can use the default Silverman bandwidth
      self$silverman()
      bandwidth = private$bandwidth
      cli::cli_alert_warning("use default Silverman bandwidth")
    }
  }
  if(!flag) return(invisible(self)) # return

  # inform the user of what to use
  cli::cli_alert_info("using max_m = {max_m} and bandwidth = {bandwidth}")

  iT = private$iT; ip = private$ip

  # Initialize storage
  V <- numeric(max_m)
  penalty <- numeric(max_m)
  IC_values <- numeric(max_m)

  # Loop over possible number of factors (R)
  for (mi in 1:max_m) {
    residuals <- matrix(NA, nrow = iT, ncol = ip)
    prev_F = NULL
    for (r in 1:iT){
      # Step 1: Perform PCA with R factors
      pca_result <- try(local_pca(private$data, r = r, bandwidth = bandwidth,
                                  m = mi, kernel_func = epanechnikov_kernel,
                                  prev_F))
      if("try-error" %in% class(pca_result))
      {
        next
      }

      X_r <- matrix(0, nrow = iT, ncol = ip)
      X_r <- sweep(private$data, 1, sqrt(pca_result$w_r), `*`)
      scaled_loadings <- sqrt(ip) * sweep(pca_result$loadings, 2, sqrt(colSums(pca_result$loadings^2)), "/")
      Lambda_breve_R <- t((1/(iT*ip))*t(X_r)%*%X_r%*%scaled_loadings)
      F_breve_R <- solve((Lambda_breve_R)%*%t(Lambda_breve_R))%*%(Lambda_breve_R)%*%private$data[r,]

      # Step 2: Compute SSR (Sum of Squared Residuals)
      residuals[r,] <- private$data[r,] - t(F_breve_R) %*% (Lambda_breve_R)

      prev_F <- pca_result$F_hat_r
    }
    V[mi] <- sum(residuals^2) / (ip * iT)
    penalty[mi] <- mi * ((ip+iT*bandwidth)/(ip*iT*bandwidth))*log((ip*iT*bandwidth)/(ip+iT*bandwidth))
    IC_values[mi] <- log(V[mi]) + penalty[mi]
  }
  # Step 4: Determine optimal number of factors
  optimal_m <- which.min(IC_values)

  # results and return
  #message(sprintf("Optimal number of factors is %s.", optimal_m))
  private$optimal_m = optimal_m
  private$IC_values = IC_values
  # we return the optimal m and also save it in the object
  return(optimal_m)
})
