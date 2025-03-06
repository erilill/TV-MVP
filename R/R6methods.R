## In this file, we define the methods of the R6 class


TVMVP$set("public", "set_data", function(data) {
  # for set_data function, the argument data should not be missing
  # this is why I use missing(data) for checking without setting any default value of it
  # note that this is different from initialize function where data can be missing with default value NULL
  if(missing(data)){
    cli::cli_alert_warning("The data is empty.")
  } else{
    private$data <- data
    private$iT <- nrow(data)
    private$ip <- ncol(data)
    
    #tmp_size <- get_object_size(private$data)
    tmp_size <- prettyunits::pretty_bytes(object.size(private$data))
    cli::cli_alert_info("Tibble data set {tmp_size} with {private$iT} rows and {private$ip} columns successfully assigned.")
  }
  
  invisible(self)  # Enables chaining
})


TVMVP$set("public", "get_data", function() {
  if(is.null(private$data)){
    cli::cli_alert_warning("The data is empty.")
    return(invisible(NULL))
  }
  return(private$data)
  # make a copy of the data and return it
})