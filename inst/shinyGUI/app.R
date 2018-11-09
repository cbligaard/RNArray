#### Shiny app for the RNArray pipeline

#### App code by Christina Bligaard Pedersen, 2018

library(shiny)
library(shinyFiles)
library(fs)

# Defining UI for RNArray app ----
ui <- fluidPage(

  # App title ----
  h1(id="big-heading", "The RNArray pipeline"),
  tags$style(HTML("#big-heading{color: #cb4154;}")),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(h3("File selection:"),

      # Input: Selection of RNA-seq fastq file(s) to process ----
      h5("RNA-Seq FASTQ file selection:"),
      shinyFilesButton("rnaseq_files", "Browse...", "RNA-Seq FASTQ file selection", multiple = T),
      br(),
      br(),
      h5("Microarray CEL file selection:"),
      shinyFilesButton("microarray_files", "Browse...", "Microarray CEL file selection", multiple = T),
      br(),
      br(),
      radioButtons("ref_specs", h3("Reference selection:"), selected = 'detect',
                   list("Use local file" = "file", "Pick from list" = "download", "Auto-detection" = "detect")),
      conditionalPanel(
        condition = "input.ref_specs=='file'",
        h5("Reference file selection:"),
        shinyFilesButton("ref_file", "Browse...", "Reference file selection", multiple = F)),
      conditionalPanel(
        condition = "input.ref_specs=='download'",
        selectInput("ref_platform", label = "Array platform:", choices = c("Affymetrix HG-U133 Plus 2.0", "Affymetrix HG-U133A", "Affymetrix HG-U133B", "Affymetrix HG-U133A Plus 2.0", "Affymetrix HG-U95A", "Affymetrix HG-U95A v. 2", "Affymetrix HG-U95B", "Affymetrix HG-U95C", "Affymetrix HG-U95D", "Affymetrix HG-U95E", "Affymetrix Mouse Genome 430 2.0", "Affymetrix Mouse Genome 430A 2.0"))),


      # Help text
      br(),
      br(),
      p("FASTQ file(s) must have one of the extensions: '.fq', '.fastq', '.fq.gz', '.fastq.gz' or '.fgz'."),
      p("CEL file(s) must have the extension: '.CEL'"),
      p("A reference file should contain the target sequences for your array - we also allow upload of pre-made indeces from kallisto.
        Allowed extensions are: '.idx', '.fasta', '.fsa', '.fas', or '.fa'.")

    ),

    # Main panel for displaying outputs ----
    mainPanel(

      # Make something about the pairing - so it can either choose for you or if it's hopeless names, you can specify the pairs

      br(),
      fluidRow(
        # Pipeline settings - including/excluding different part of the pipeline
        column(5, h3("Read specifications:"),
               checkboxInput("read_type", "FASTQ-files are paired-end", TRUE),
               # If paired-end
               conditionalPanel(
                 condition = "input.read_type",
                 checkboxInput('give_length', 'Specify fragment length', FALSE),
                 conditionalPanel(
                  condition = "input.give_length",
                  numericInput("frag_length_pe", "Fragment length", value = 100, min = 1)
               )),
               conditionalPanel(
                 condition = "input.read_type",
                 checkboxInput('give_sd', 'Specify fragment length standard deviation', FALSE),
                 conditionalPanel(
                   condition = "input.give_sd",
                   numericInput("frag_sd_pe", "Fragment length standard deviation", value = 10, min = 0)
                )),

               # If single-end
               conditionalPanel(
                 condition = "!input.read_type",
                  numericInput("frag_length_se", "Fragment length", value = 100, min = 1),
                 numericInput("frag_sd_se", "Fragment length standard deviation", value = 10, min = 0)
                 )),

        # Execution settings - run now or save script for running later
        column(5, radioButtons("exec_specs", h3("Execution options:"),
                    list("Quit GUI and run the pipeline in R" = "run", "Only save command for later" = "write")),
               numericInput("n_cores", "Number of cores", value = 1, min = 1, max = 28))
      ),

      # Output: Formatted text for caption ----
      br(),
      textOutput("n_files", container = span),

      br(),
      br(),
      uiOutput("rna_paths"),
      br(),
      uiOutput("microarray_paths"),
      br(),
      textOutput("ref_path"),

      # Submit button
      br(),
      br(),
      uiOutput("submit")
    )
  )
)



# Define server logic ----
server <- function(input, output, session) {

  # Defining a list of prettier names for Affymetrix arrays to match with the automatic detection
  platform_names <- list('HG-U133_Plus_2' = 'Affymetrix HG-U133 Plus 2.0',
                         'HG-U133A' = 'Affymetrix HG-U133A',
                         'HG-U133A_2' = 'Affymetrix HG-U133A Plus 2.0',
                         'HG-U133B' = 'Affymetrix HG-U133B',
                         'HG_U95A' = 'Affymetrix HG-U95A',
                         'HG_U95Av2' = 'Affymetrix HG-U95A v. 2',
                         'HG_U95B' = 'Affymetrix HG-U95B',
                         'HG_U95C' = 'Affymetrix HG-U95C',
                         'HG_U95D' = 'Affymetrix HG-U95D',
                         'HG_U95E' = 'Affymetrix HG-U95E',
                         'Mouse430A_2' = 'Affymetrix Mouse Genome 430A 2.0',
                         'Mouse430_2' = 'Affymetrix Mouse Genome 430 2.0')

  # Input files
  volumes <- c(Home = fs::path_home(), getVolumes()())
  shinyFileChoose(input, "rnaseq_files", roots = volumes, session = session, defaultPath = 'Documents', filetypes = c('gz', 'fastq', 'fq', 'fqz'))
  shinyFileChoose(input, "microarray_files", roots = volumes, session = session, defaultPath = 'Documents', filetypes = c('CEL'))
  shinyFileChoose(input, "ref_file", roots = volumes, session = session, defaultPath = 'Documents', filetypes = c("idx", "fasta", "fas", "fa", "fsa"))

  # Reactive inputs ----
  rna_input <- reactive({parseFilePaths(volumes, input$rnaseq_files)})
  array_input <- reactive({parseFilePaths(volumes, input$microarray_files)})
  ref_input <- reactive({parseFilePaths(volumes, input$ref_file)})
  ref_platform <- reactive({input$ref_platform})
  ref_specs <- reactive({input$ref_specs})
  ref_platform_auto <- reactiveValues(p = 0, file = 0)

  # Rendering information about input
  output$rna_paths <- renderText({
    if (length(input$rnaseq_files) > 1) {
      HTML(paste('Selected RNA-seq files: ', paste(rna_input()[,1][[1]], collapse = ', '), sep="<br/>"))
    }
  })
  output$microarray_paths <- renderText({
    if (length(input$microarray_files) > 1) {
      HTML(paste('Selected microarray files: ', paste(array_input()[,1][[1]], collapse = ', '), sep="<br/>"))

    }
  })
  output$ref_path <- renderText({
    if (length(input$ref_file) == 2 && input$ref_specs == 'file') {
      paste0('Selected reference file: ', ref_input()[,1][[1]])
    } else if (input$ref_specs == 'download') {

      # Detecting microarray platform
      ref_platform_det = affy::ReadAffy(filenames = array_input()[,'datapath'][[1]][1])@cdfName
      if (ref_platform_det %in% names(platform_names)) {
        ref_platform_auto$p <- platform_names[ref_platform_det]
      }

      paste0('Selected reference: ', ref_platform())

    } else if (input$ref_specs == 'detect' & nrow(array_input()) != 0) {

      # Detecting microarray platform
      ref_platform_det = affy::ReadAffy(filenames = array_input()[,'datapath'][[1]][1])@cdfName

      # Checking if there's a pre-made index or at least a pre-downloaded fasta in either the microarray of RNA-seq folder
      rna_seq_folder_f = list.files(dirname(rna_input()[,'datapath'][[1]][1]), pattern = '.(fa|fsa|fasta|fas|idx)$', full.names = T)
      array_folder_f = list.files(dirname(array_input()[,'datapath'][[1]][1]), pattern = '.(fa|fsa|fasta|fas|idx)$', full.names = T)

      if (length(grep('idx$', rna_seq_folder_f)) == 1) {
        ref_platform_auto$file = rna_seq_folder_f[grep('idx$', rna_seq_folder_f)]
      } else if (length(grep('idx$', array_folder_f)) == 1) {
        ref_platform_auto$file = array_folder_f[grep('idx$', array_folder_f)]
      } else if (length(grep('.(fa|fsa|fasta|fas)$', rna_seq_folder_f)) == 1) {
        ref_platform_auto$file = array_folder_f[grep('(fa|fsa|fasta|fas)$', array_folder_f)]
      } else if (length(grep('(fa|fsa|fasta|fas)$', array_folder_f)) == 1) {
        ref_platform_auto$file = array_folder_f[grep('(fa|fsa|fasta|fas)$', array_folder_f)]
      }

      if (ref_platform_det %in% names(platform_names)) {
        ref_platform_auto$p <- platform_names[ref_platform_det]

        if (ref_platform_auto$file == 0) {
          paste0('Reference is automatically detected to be ', ref_platform_auto$p, '. Could not determine a pre-downloaded reference file to use.')
        } else {
          paste0('Reference is automatically detected to be ', ref_platform_auto$p, '. Found pre-downloaded reference: ', ref_platform_auto$file, ' - if this is not the right file, please remove that file from the folder or choose another option from the menu.')
        }

      } else {
        ref_platform_auto$p <- 0
        paste0('There is no available reference file for the detected micraorray platform.')
      }
    }
  })


  # Monitoring input regarding length and standard deviation for reads
  frag_length_valid <- reactive({
    ((input$read_type && !is.null(input$frag_length_pe) && input$frag_length_pe > 0 && input$frag_length_pe %% 1 == 0 && is.numeric(input$frag_length_pe)) || (!input$give_length && input$read_type) || (!input$read_type && !is.null(input$frag_length_se) && input$frag_length_se > 0 && input$frag_length_se %% 1 == 0 && is.numeric(input$frag_length_se)))
  })

  frag_sd_valid <- reactive({
    ((input$read_type && !is.null(input$frag_sd_pe) && input$frag_sd_pe > 0 && is.numeric(input$frag_sd_pe)) || (!input$give_sd && input$read_type) || (!input$read_type && !is.null(input$frag_sd_se) && input$frag_sd_se > 0 && is.numeric(input$frag_sd_se)))
  })


  # Printing overall information about input
  output$n_files <- renderText({
    paste0("You have selected ", max(0, nrow(rna_input())), " " , ifelse(input$read_type,'paired','single'), "-end RNA-Seq file(s) and ", max(0, nrow(array_input())), " microarray file(s) for processing.")
  })


  # Printing information about output
  output$submit = renderUI({
    errors = error_checks(input, frag_length_valid, frag_sd_valid, rna_input, array_input, ref_input, ref_specs, ref_platform_auto, ref_platform)
    if ((!is.null(input$RNArray_submit)) && (input$RNArray_submit)) {
      if (length(errors) == 0) {
        write_exec_script(input, rna_input, array_input, ref_input, ref_platform, ref_specs, ref_platform_auto)
      } else {
        stop(simpleError("Cannot write command with errors."))
      }
      if (input$exec_specs %in% c("write", "run")){
        RNArray_submit <<- FALSE
        if (input$exec_specs == "run"){
          RNArray_submit <<- TRUE
          out_dir <<- file.path(dirname(rna_input()[,'datapath'][[1]][1]), "RNArray_output")
        }
        stop(simpleWarning("RNArray GUI Setup Complete."))
      }
    }
    else {
      if (length(errors) == 0){
        return(actionButton(inputId="RNArray_submit", label="Execute!"))
      } else {
        return(HTML(paste("The following problems must be corrected before running:", paste(errors, collapse = "<br/>"), sep="<br/>")))
      }
    }
  })
}


# Check for errors in input data
error_checks <- function(input, frag_length_valid, frag_sd_valid, rna_input, array_input, ref_input, ref_specs, ref_platform_auto, ref_platform) {
  errors <- c()

  if (nrow(rna_input()) == 0) {
    errors <- c(errors, "- No RNA-seq files selected.")
  }
  else {
    if (nrow(rna_input()) %% 2 != 0 & input$read_type) {
    errors <- c(errors, "- You selected an odd number of FASTQ files, but claim they are paired-end reads.")
    }
    if (any(tools::file_ext(rna_input()[,1]) == 'gz')) {
      first_ext <- tools::file_ext(tools::file_path_sans_ext(rna_input()[,1][tools::file_ext(rna_input()[,1]) == 'gz']))
      if (any(first_ext != 'fastq' & first_ext != 'fq')) {
        errors <- c(errors, "- You have selected invalid gzipped files for upload. They must either have the extension fastq.gz or fq.gz.")
      }
    }
  }
  if (nrow(array_input()) == 0) {
    errors <- c(errors, "- No microarray files selected.")
  }
  if (ref_specs() == 'file' & nrow(ref_input()) == 0) {
    errors <- c(errors, "- No reference file selected.")
  }
  if (ref_specs() == 'download' & nrow(array_input()) != 0 & ref_platform_auto$p[[1]] == 0 | (ref_platform() != ref_platform_auto$p[[1]])) {
    errors <- c(errors, "- The specified reference does not fit with the microarray files.")
  }
  if (ref_platform_auto$p[[1]] == 0 & ref_specs() == 'detect' & nrow(array_input()) != 0) {
    errors <- c(errors, "- Automatic detection of microarray platform not possible, because no pre-made index matches the detected platform.")
  }
  if (!frag_length_valid()) {
    errors <- c(errors, "- Invalid fragment length, it must be a positive integer.")
  }
  if (!frag_sd_valid()) {
    errors <- c(errors, "- Invalid fragment length standard deviation, it must be a positive number.")
  }
  if (is.na(input[['n_cores']]) | input[['n_cores']]  %% 1 != 0 | input[['n_cores']] < 1 | input[['n_cores']] > 28) {
    errors <- c(errors, "- Number of cores must be a positive integer between 1 and 28 (both inclusive).")
  }

  return(errors)
}

# Script for actual data processing ----
write_exec_script <- function(input, rna_input, array_input, ref_input, ref_platform, ref_specs, ref_platform_auto) {
  # Read input
  input_data = reactiveValuesToList(input)

  # Make folder for output - default is RNA-seq file location
  out_dir <- file.path(dirname(rna_input()[,'datapath'][[1]][1]), "RNArray_output")
  if (!file.exists(out_dir)){
    dir.create(out_dir, showWarnings = F)
  }

  # Generate command for processing to use in another R session and save that file

  # Getting final fragment length
  if (input_data[['read_type']]) {
    if (input_data[['give_length']]) {
      frag_length <- input_data[['frag_length_pe']]
    } else {
      frag_length <- NA
    }
    if (input_data[['give_sd']]) {
      frag_sd <- input_data[['frag_sd_pe']]
    } else {
      frag_sd <- NA
    }
  } else {
    frag_length <- input_data[['frag_length_se']]
    frag_sd <- input_data[['frag_sd_se']]
  }

  # Making command
  cmd_out <- file(file.path(out_dir, 'RNArray_command.txt'), 'w')

  # Make list of files to process
  if (length(rna_input()[,'datapath'][[1]]) == 1) {
    rna_list <- paste0('"', rna_input()[,'datapath'][[1]], '"')
  } else {
    rna_list <- paste(rna_input()[,'datapath'], collapse=", ")
  }

  if (length(array_input()[,'datapath'][[1]]) == 1) {
    array_list <- paste0('"', array_input()[,'datapath'][[1]], '"')
  } else {
    array_list <- paste(array_input()[,'datapath'], collapse=", ")
  }

  if (ref_specs() == 'file') {
    ref_file_out <- ref_input()[,'datapath'][[1]]
  } else if (ref_specs() == 'download') {
    ref_file_out <- ref_platform()
  } else if (ref_specs() == 'detect' && ref_platform_auto$file != 0) {
    ref_file_out <- ref_platform_auto$file
  } else {
    ref_file_out <- ref_platform_auto$p[[1]]
  }

  # Write to file
  writeLines(paste0("### Command generated for running RNArray by the RNArray GUI. Command generated ", date(), "."), cmd_out)
  writeLines(paste0("RNArray_output <- RNArray::RNArray(rna_seq_files = ", rna_list,
         ", array_files = ", array_list,
         ", ref_file = '", ref_file_out, "', paired_end = ", input_data[['read_type']],
         ", frag_length = ", frag_length, ", frag_sd = ", frag_sd,
         ", cores = ", input_data[['n_cores']], ", kallisto_cmd = 'kallisto', tech_prefix = 'AFFX', out_dir = '", out_dir, "')"), cmd_out)
  close(cmd_out)

}


shinyApp(ui, server)
