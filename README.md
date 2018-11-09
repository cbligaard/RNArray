# Welcome to the RNArray repository!

The package code was written by [Christina Bligaard Pedersen](https://www.dtu.dk/english/service/phonebook/person?id=88697&cpid=176923) and [Lars Rønn Olsen](https://www.dtu.dk/english/service/phonebook/person?id=26586).

The RNArray pipeline is a generalized version of the method from “Using microarray-based subtyping methods for breast cancer in the era of high-throughput RNA sequencing” by Pedersen C. B., Nielsen, F. C., Rossing, M., and Olsen, L. R., 2018 ([DOI:10.1002/1878-0261.12389](https://doi.org/10.1002/1878-0261.12389)). The stand-alone implementation of the tool is currently in submission.

To install the package, first clone the package to a local path and run the following commands in R:
```
install.packages("devtools")
library(devtools)
devtools::install_github("cbligaard/RNArray")
```

You also need to install kallisto (v. >= 0.42.2 required). It be installed from [The Pachter Lab](https://pachterlab.github.io/kallisto/download). 

To use the graphical user interface for RNArray, run the following command:
```
RNArray::RNArray_UI()
```

The interface will open in a browser window upon command execution.
This will give you the primary functionalities of the RNArray function, but with a user-friendly interface to aid file selection and running for those not familiar/comfortable with scripting. 

There are two options for you run 'Quit GUI and run the pipeline in R' and 'Only save command for later'. The first will save the resulting command for running RNArray in a file in the output directory (RNArray_command.txt) and execute it immidiately. The latter will only write the file with the command. When using the GUI, the output directory is the location of the RNA-seq file(s).


#### User guide
Paired-end RNA-seq files should be named such that a pair is specified by `<prefix>_1.fq` and `<prefix>_2.fq` OR `<prefix>_R1.fq` and `<prefix>_R2.fq`. Valid file extensions are .fastq.gz, .fq.gz, .fgz .fastq, and .fq.

Microarray files should be from an Affymetrix platform and must have the extension .CEL and the reference file should have one of the extensions .idx, .fasta, .fa, .fsa, or .fas. If .idx is used, the file is expected to be a kallisto index. As an alternative to providing the reference file, the RNArray program can also auto-detect the platform, or you can choose from a drop-down menu.

If you choose to upload a reference file, the options are endless, but the other selection options only support the following platforms:

* [Affymetrix HGU-133 Plus 2.0](http://www.affymetrix.com/support/technical/byproduct.affx?product=hg-u133-plus)
* [Affymetrix HGU-133A](http://www.affymetrix.com/support/technical/byproduct.affx?product=hgu133-20)
* [Affymetrix HGU-133A Plus 2.0](http://www.affymetrix.com/support/technical/byproduct.affx?product=hgu133-20) 
* [Affymetrix HGU-133B](http://www.affymetrix.com/support/technical/byproduct.affx?product=hgu133-20)
* [Affymetrix HGU-95A-E](http://www.affymetrix.com/support/technical/byproduct.affx?product=hgu95)
* [Affymetrix HGU-95A v2](http://www.affymetrix.com/support/technical/byproduct.affx?product=hgu95) 
* [Affymetrix Mouse Genome 430 2.0](http://www.affymetrix.com/support/technical/byproduct.affx?product=moe430-20)
* [Affymetrix Mouse Genome 430A 2.0](http://www.affymetrix.com/support/technical/byproduct.affx?product=moe430A-20)

If you choose 'Auto-detection' of the reference, the program will also look in the folders containing the microarray and RNA-seq files to identify pre-downloaded reference files - if these are not found, the file will be downloaded during the run.
If you choose 'Pick from list' for the reference, the reference file will be downloaded in any case.

*Last updated on November 7, 2018.*