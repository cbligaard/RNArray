#' RNArray pipeline for aligning RNA-seq reads to microarray probesets and making distributions comparable
#'
#' @param rna_seq_files VECTOR The RNA-seq FASTQ-files of interest (they must have one of the extensions fastq.gz, fq.gz, fqz, fq, or fastq). If your files are paired-end, they should follow the naming scheme: PairA_1.fq, PairA_2.fq, PairB_1.fq, PairB_2.fq etc. or PairA_R1/2.fq etc.
#' @param array_files VECTOR The microarray CEL-files of interest (they must have the extension CEL).
#' @param ref_file STRING The name of the fasta file (or kallisto .idx file) containing the reference sequences for you array of interest. ALternatively, the name of the platform for you array of interest.
#' @param paired_end LOGICAL, either true (RNA-seq reads are paired-end, default) or false (reads are single-end).
#' @param frag_length INTEGER Fragment length of the RNA-seq. Required when paired_end = F, optional when paired_end = T.
#' @param frag_sd NUMERIC Estimated standard deviation of read fragment length.  Required when paired_end = F, optional when paired_end = T.
#' @param cores (Optional) Number of cores to use - multithreading significantly increases the speed (default = 1).
#' @param kallisto_cmd (Optional) STRING giving full command to use to call kallisto, if simply typing "kallisto" at the command line does not give the required version of kallisto or does not work. Default is simply "kallisto". If used, this argument should give the full path to the desired kallisto binary.
#' @param tech_prefix (Optional) STRING giving the prefix of technical rows on the applied platform - for Affymetrix it's 'AFFX'.
#' @param out_dir (Optional) STRING giving an output directory. If not specified, the results will be saved in a sub-directory (RNArray_output) in the directory containing the reference file.
#'
#' @return Quantile normalized and batch corrected RNA-seq and microarray data
#'
#' @author Christina Bligaard Pedersen and Lars Rønn Olsen
#'
#' @examples
#' RNArray(c('Desktop/FASTQ/SampleA_1.fq.gz', 'Desktop/FASTQ/SampleA_2.fq.gz', 'Desktop/FASTQ/SampleB_1.fq.gz', 'Desktop/FASTQ/SampleB_2.fq.gz'), c('Desktop/CEL/FileX.CEL', 'Desktop/CEL/FileY.CEL'), 'Desktop/my_targets.fa', paired_end = T, cores = 4)
#' RNArray(c('Desktop/FASTQ/FileA.fq.gz', 'Desktop/FASTQ/FileB.fq.gz'), c('Desktop/CEL/FileX.CEL', 'Desktop/CEL/FileY.CEL'), 'Desktop/my_targets.fa', paired_end = F, frag_length = 89, frag_sd = 7, cores = 2)
#'
#' @export
RNArray <- function(rna_seq_files,
                    array_files,
                    ref_file,
                    paired_end = T,
                    frag_length = NA,
                    frag_sd = NA,
                    cores = 1,
                    kallisto_cmd = 'kallisto',
                    tech_prefix = 'AFFX',
                    out_dir = NA
) {

  cat("RNArray run started", date(), "\n")

  # Premade index location
  platform_references <- list('Affymetrix HG-U133 Plus 2.0' = 'https://bitbucket.org/cbligaard/rnarray-files/downloads/HG-U133_Plus_2.idx',
                              'Affymetrix HG-U133A Plus 2.0' = 'https://bitbucket.org/cbligaard/rnarray-files/downloads/HG-U133A_2.idx',
                              'Affymetrix HG-U133A' = 'https://bitbucket.org/cbligaard/rnarray-files/downloads/HG-U133A.idx',
                              'Affymetrix HG-U133B' = 'https://bitbucket.org/cbligaard/rnarray-files/downloads/HG-U133B.idx',
                              'Affymetrix HG-U95A' = 'https://bitbucket.org/cbligaard/rnarray-files/downloads/HG-U95A.idx',
                              'Affymetrix HG-U95A v. 2' = 'https://bitbucket.org/cbligaard/rnarray-files/downloads/HG-U95A_2.idx',
                              'Affymetrix HG-U95B' = 'https://bitbucket.org/cbligaard/rnarray-files/downloads/HG-U95B.idx',
                              'Affymetrix HG-U95C' = 'https://bitbucket.org/cbligaard/rnarray-files/downloads/HG-U95C.idx',
                              'Affymetrix HG-U95D' = 'https://bitbucket.org/cbligaard/rnarray-files/downloads/HG-U95D.idx',
                              'Affymetrix HG-U95E' = 'https://bitbucket.org/cbligaard/rnarray-files/downloads/HG-U95E.idx',
                              'Affymetrix Mouse Genome 430A 2.0' = 'https://bitbucket.org/cbligaard/rnarray-files/downloads/Mouse430A_2.idx',
                              'Affymetrix Mouse Genome 430 2.0' = 'https://bitbucket.org/cbligaard/rnarray-files/downloads/Mouse430_2.idx')

  # Minimum number of files allowed
  n_min_a = 2
  n_min_r = 2

  # Checking if input is okay
  if (!paired_end && (is.na(frag_length) || is.na(frag_sd))) {
    stop('When using single-end RNA-seq reads, you must specify the fragment length and standard deviation.')
  } else if (cores %% 1 != 0 || cores <= 0) {
    stop('The number of cores must be a positive integer.')
  } else if (!is.character(ref_file) || ((tools::file_ext(ref_file) != 'fa') && (tools::file_ext(ref_file) != 'fasta') && (tools::file_ext(ref_file) != 'idx') && !(ref_file %in% names(platform_references)))) {
    stop('Your reference file must be a string with the extension .fa, .fasta, or .idx or it must be in the list of valid platforms.')
  } else if (!is.vector(rna_seq_files) || length(rna_seq_files) == 0) {
    stop('No FASTQ files were specified. This must be a vector.')
  } else if (length(rna_seq_files) != length(grep('(*.fastq.gz)|(*.fq.gz)|(*.fqz)|(*.fq)|(*.fastq)', rna_seq_files))) {
    stop('At least one of your RNA-seq files does not have a valid extension.')
  } else if (!is.vector(array_files) || length(array_files) == 0) {
    stop('No CEL files were specified. This must be a vector.')
  } else if (length(array_files) < n_min_a || ((paired_end) && length(rna_seq_files) < n_min_r*2) || (!(paired_end) && length(rna_seq_files) < n_min_r)) {
    stop(paste0('The number of microarray samples should be at least ', n_min_a, ' and the number of RNA-seq samples should be at least ', n_min_r, '.'))
  } else if (!all(tools::file_ext(array_files)=='CEL')) {
    stop('At least one of your microarray files does not have the correct extension.')
  } else if ((paired_end) && (length(rna_seq_files) %% 2 != 0)) {
    stop('You specified paired-end reads, but the number of input files is uneven.')
  } else if (!is.character(kallisto_cmd) || nchar(kallisto_cmd) == 0) {
    stop('kallisto command must be a string.')
  } else if (!is.character(tech_prefix)) {
    stop('Tech prefix must be a string.')
  } else if (!is.na(out_dir) && !is.character(out_dir)) {
    stop('If you specify an output directory, it must be a character.')
  }


  # Checking that it is possible to run kallisto on command line and that its version is >= 0.42.2
  cat("\nChecking kallisto installation.\n\n")
  if (package_version(strsplit(system(paste(kallisto_cmd, 'version'), intern = T), ' ')[[1]][3]) < '0.42.2') {
    stop('kallisto version < 0.42.2. Please upgrade your installation!')
  }

  # Download reference if only platform is provided
  if (ref_file %in% names(platform_references)) {
    dest <- file.path(dirname(rna_seq_files)[1], paste0(gsub(' ', '_', ref_file), '.idx'))

    # Get the technical prefix
    if (startsWith(ref_file, 'Affymetrix')) {
      tech_prefix = 'AFFX'
    }

    download.file(url = platform_references[[ref_file]], destfile = dest)
    ref_file <- dest # Path to downloaded index - the indeces are rather big compared to raw fastas...
  }

  # Setting an output directory if not provided
  if (is.na(out_dir)) {
    out_dir <- file.path(dirname(rna_seq_files[1]), 'RNArray_output')
    dir.create(out_dir, showWarnings = F)
  } else {
    if (!dir.exists(out_dir)) {
      stop('Invalid output directory specified.')
    }
  }


  ### Run microarray data pre-processing ----

  # Reading data and RMA normalization
  cat("Reading and processing microarray data.")
  array_data <- affy::ReadAffy(filenames = array_files); array_type <- array_data@cdfName
  array_data <- affy::exprs(affy::rma(array_data))

  # Get non-technical rows
  include_rows <- rownames(array_data)[!startsWith(rownames(array_data), tech_prefix)]


  ### Run kallisto pseudo-alignment for all RNA-seq files ----
  cat("\nReading and processing RNA-seq data.\n")

  # Check if the uploaded file is already an index - if not, make one'
  if (tools::file_ext(ref_file) != 'idx') {
    cat("Making kallisto index.\n")
    idx_file <- paste0(dirname(rna_seq_files)[1], '/', basename(tools::file_path_sans_ext(ref_file)), '.idx')
    system(paste('kallisto index -i', idx_file, ref_file))
    ref_file <- idx_file
  }

  # Stop if index was not made
  stopifnot(file.exists(ref_file))

  # Make sure kallisto index contains as many probe sets as the microarray
  n_probes_idx = strsplit(system(paste('kallisto inspect', ref_file), intern = T)[3], ' ')[[1]][6]
  if (n_probes_idx != nrow(array_data)) {
    stop('The number of probe sets in the index does not match the number of probe sets on the microarray platform.')
  }


  # Making the targets_file for kallisto and running the program - two options, paired-end or single-end reads
  cat("Pseudo-aligning reads with kallisto.\n")
  rna_seq_files_table = file.path(dirname(rna_seq_files)[1], 'tmp_rna_seq_files_table.txt')
  file_bases <- basename(unique(gsub('_(1|2|R1|R2).((fastq.gz)|(fq.gz)|(fqz)|(fq)|(fastq))', '', rna_seq_files)))

  if (paired_end) {
    pairs <- lapply(file_bases, function(x) rna_seq_files[startsWith(basename(rna_seq_files), paste0(x, '_'))])

    # Check that it makes sense
    if (any(sapply(pairs, function(x) {length(x) != 2}))) {
      stop('Your paired-end RNA-seq files do not conform to the expected naming standard.')
    }

    write.table(x = cbind(file_bases, matrix(basename(unlist(pairs)), byrow = T, ncol = 2)), file = rna_seq_files_table, quote = F, sep = '\t', row.names = F, col.names = c('name', 'file1', 'file2'))

    # Actually running kallisto
    if (is.na(frag_length) && is.na(frag_sd)) {
      kallisto_log <- scater::runKallisto(targets_file = rna_seq_files_table, transcript_index = ref_file, single_end = F, n_cores = cores, output_prefix = file.path(out_dir, 'kallisto'), verbose = T, kallisto_cmd = kallisto_cmd)
    } else if (is.na(frag_sd)) {
      kallisto_log <- scater::runKallisto(targets_file = rna_seq_files_table, transcript_index = ref_file, single_end = F, n_cores = cores, output_prefix = 'kallisto', fragment_length = frag_length, verbose = T, kallisto_cmd = kallisto_cmd)
    } else if (is.na(frag_length)) {
      kallisto_log <- scater::runKallisto(targets_file = rna_seq_files_table, transcript_index = ref_file, single_end = F, n_cores = cores, output_prefix = 'kallisto', fragment_standard_deviation = frag_sd, verbose = T, kallisto_cmd = kallisto_cmd)
    } else {
      kallisto_log <- scater::runKallisto(targets_file = rna_seq_files_table, transcript_index = ref_file, single_end = F, n_cores = cores, output_prefix = 'kallisto', fragment_length = frag_length, fragment_standard_deviation = frag_sd, verbose = T, kallisto_cmd = kallisto_cmd)
    }

  } else {
    write.table(x = cbind(file_bases, rna_seq_files), file = rna_seq_files_table, quote = F, sep = '\t', row.names = F, col.names = c('name', 'file'))

    # Actually running kallisto
    kallisto_log <- scater::runKallisto(targets_file = rna_seq_files_table, transcript_index = ref_file, single_end = T, n_cores = cores, output_prefix = file.path(out_dir, 'kallisto'), fragment_length = frag_length, fragment_standard_deviation = frag_sd, verbose = T, kallisto_cmd = kallisto_cmd)

  }

  if (length(grep('had status 1', toString(kallisto_log[[1]]$kallisto_log)))==1) {
    cat('kallisto started with command(s):', paste(sapply(kallisto_log, function(x){x$kallisto_call}), '\n'), '\n')
    stop('kallisto did not run properly.')
  }

  # Removing kallisto file table after use to avoid clutter
  system(paste0('rm ', rna_seq_files_table))

  # Provide statistics about mapping (make prettier)
  cat(toString(kallisto_log[[1]]$kallisto_log[4]), "\n")
  cat(toString(lapply(kallisto_log, function(x) x$kallisto_log[12])), "\n")


  # Read alignment data into R for downstream processing ----
  mapping_folders <- paste0(out_dir, '/kallisto_', file_bases, '/abundance.tsv')
  mapping_data <- plyr::ldply(mapping_folders, function(fn) data.frame(Filename=fn, read.table(fn, header = T)))
  RNAseq_TPM <- matrix(mapping_data$tpm, ncol = length(file_bases))
  colnames(RNAseq_TPM) <- file_bases
  rownames(RNAseq_TPM) <- gsub(paste0('.+:(', tech_prefix, '-.+?|[0-9]+.+?)(_at).+'), '\\1\\2', levels(mapping_data$target_id)) # regex may fail in some cases - check it and fix it so it works for homemade indeces and my pre-made versions.

  # Checking that rownames correspond to those of the microarray
  if (!all(include_rows %in% rownames(RNAseq_TPM))) {
    stop(paste0('It seems like your microarray data is not from the same platform as the selected reference file. ', table(include_rows %in% rownames(RNAseq_TPM))[[1]], ' of the micraorray probes was not found in the resulting RNA-seq data.'))
  }


  #### Quantile normalization of the two data types ----
  cat("\nPerforming quantile normalization and batch correction.\n")
  RNAseq_qnorm <- preprocessCore::normalize.quantiles.use.target(RNAseq_TPM[include_rows,], target = rowMeans(array_data[include_rows,]))
  rownames(RNAseq_qnorm) <- rownames(RNAseq_TPM[include_rows,]); colnames(RNAseq_qnorm) <- colnames(RNAseq_TPM)
  cat("- Quantile normalization done.\n")


  ### Batch correction of the data ----
  batch <- sva::ComBat(cbind(RNAseq_qnorm, array_data[include_rows,]),c(rep(1,ncol(RNAseq_qnorm)),rep(2,ncol(array_data))))
  RNAseq_qnorm_bc <- batch[,1:ncol(RNAseq_qnorm)]
  array_bc <- batch[,(ncol(RNAseq_qnorm)+1):(ncol(RNAseq_qnorm)+ncol(array_data))]
  cat("- Batch correction done.\n\n")


  ### Comparing distributions by making a density plot ---
  png(filename = file.path(out_dir, 'Density_plot.png'), res = 300, height = 2000, width = 3000)
  density_plot(array_bc, RNAseq_qnorm_bc)
  dev.off()
  cat('A plot showing the data distributions has been saved in the output directory.')

  ### Outputting results to files ---
  write.table(array_bc, file = file.path(out_dir, 'Microarray.txt'), sep = '\t', quote = F, row.names = T, col.names = T)
  write.table(RNAseq_qnorm_bc, file = file.path(out_dir, 'RNAseq.txt'), sep = '\t', quote = F, row.names = T, col.names = T)

  cat('Tables with the expression data have been saved in the output directory.')

  ### Return resulting dataframes ----
  return(list('rna_seq' = RNAseq_qnorm_bc, 'microarray' = array_bc))

}


### Density plot function ---
density_plot <- function(dataset1, dataset2) {
  dens1 <- density(dataset1)
  dens2 <- density(dataset2)

  par(mar=c(5.1, 5.1, 4.1, 2.1))
  plot(1, type="n", xlim = c(min(min(dens1$x),min(dens2$x)), max(max(dens1$x),max(dens2$x))),
       ylim = c(min(min(dens1$y),min(dens2$y)), max(max(dens1$y),max(dens2$y))),
       xlab = 'Value', ylab = 'Density', cex.lab=2.5, cex.axis = 2.0, yaxt = 'n')
  axis(2, cex.axis = 2.0, mgp=c(4,1,0))
  lines(dens1, col = 'blue', lty = 1, lwd = 4.5)
  lines(dens2, col = 'red', lty = 2, lwd = 4.5)
  legend('topright', c('Microarray data', 'RNA-Seq data'), col=c("blue", "red"), lwd=3.5, lty=c(1,2), bty='n', cex = 1.5)

}
