### Cleaning up Lacto Project

This assumes you already have SRA fastq files downloaded for multiple strains within the same species (i.e. Lactobacillus iners UMB4066, UMB4068).
This code is adapted from [LactoCompare](https://github.com/putonti/LactoCompare). I cleaned up this code and present it in the form of a pipeline to automate all of the steps for user-friendly access. 


### Software

* [SRA toolkit](https://github.com/ncbi/sra-tools) 
* [Fastq-dump](https://rnnh.github.io/bioinfo-notebook/docs/fastq-dump.html)
* [SPAdes](https://github.com/ablab/spades)
* [BBDUk](https://github.com/BioInfoTools/BBMap/blob/master/sh/bbduk.sh)
* [Seqtk](https://github.com/lh3/seqtk)
* [FastANI](https://github.com/ParBLiSS/FastANI)


### Packages
* os
* argsparse
* random
* glob
* subprocess
* pandas

### Running Main Pipeline
* Running with default parameters - RUN IN SCREEN
 ```
 python3 lacto_pipeline.py --i 'path/to/fastq/files' --o 'path/to/output'
 ```
Options include --subsamples (number of subsamples per strain - default 50) and --reads_in_subsample (number of reads in each subsample - default 15000, you can adjust based on coverage).

It is recommended to specify an output folder that is empty, or to specify an output folder you would like to create.

### Output Folders of Main Pipeline
Subfolders are outputted for each strain. Within each subfolder (assuming you are using default parameters):
* ani
  * Contains random pairwise ANI calculations (from fastani) from the 50 subsampled assemblies 
* assemble
  * Contains contigs.fasta file for each of the 50 subsampled assemblies
* ss
  * Contains the subsampled paired end fastqs (n = 100) 
* trim
  * Contains the trimmed subsampled paired end fastqs (n = 100)

Additionally, outside of the subfolders:
* coallate_ani.csv
  * Contains the ani values for each strain (col: strain, row: ani)
 
### Running Rscript
* simu_anis
  * Contains the output from the lacto_pipeline.py script (coallated anis for each strain)
* test_anis
  * Contains the fastani output from your own set of strains (must be the same species as your fastqs) 
 ```
 Rscript stats.r --simu_anis '/path/to/coallate_ani.csv' --test_anis '/path/to/test_ani.csv' --out '/path/to/output'
 ```

 ### Rscript output
* density.png
  * Contains the kernel density plot of the bootstrapped ANI values for your input fastq files
* stats.txt
  * Basic statistics on your ANI values across your strains
* empirical_ps.csv
  * Contains the empirical p-values (last column) imposed on the density plot for each of your test fastANI inputs
 
