#' RNArray shiny-based graphical user interface for file selection and running RNArray.
#'
#'
#' @author Christina Bligaard Pedersen
#'
#'
#' @export
RNArray_UI <- function() {

  library("shiny")

  # Launching app
  res = tryCatch({
    runApp(appDir = file.path(system.file(package = "RNArray"), "shinyGUI"), launch.browser = T)
  }, warning = function(w){
    cat(paste(w,"\n"));
  },error = function(e){
    stop(paste("Unexpected Error:",e))
  })

  # Executing function
  if (RNArray_submit) {
    run_script <- file.path(out_dir, "RNArray_command.txt")
    cat(paste("Running RNArray pipeline from command in:", run_script, "\n"))
    log_file <- file.path(out_dir, "RNArray.log")
    cat(paste("Output log in:", log_file, "\n"))
    log = log_start(log_file)

    source(run_script)

    log_end(log)
  }
  return(paste("The output from the RNArray pipeline is in:", out_dir))

}

log_start <- function(filename) {
  file_out <- file(filename, open = 'wt')
  sink(file_out, split = T)

  return(file_out)
}

log_end <- function(file_out) {
  sink()
  close(file_out)
}


