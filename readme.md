# CRESCENT

Comprehensive RNA-seq Expression, Splicing, and Coding/non-coding Element Network Tool

## Summary

CRESCENT is a Snakemake workflow following a typical bioinformatic RNA-Seq analysis providing multiple tools for each analysis step and offering flexibility in the use of their parameters. It enables differential expression analysis (DEA) and in-depth splicing analysis through differentially alternative splicing (DAS) and differential transcript usage (DTU), functional enrichment analysis through gene ontology (GO), and provides a MultiQC comprehensive report for many of its steps. CRESCENT incorporates multiple Snakemake wrappers (Mölder et al. 2021) whenever suitable, allowing for quick integration of popular tools, ensuring version compatibility and enhancing reproducibility.

## Features

* Raw sequencing files QC
* Reads cleaning
* Mapping / salmon quantification
* BAM sorting and indexing
* Read counting
* Differential Expression Analysis (DEA) on all genomic features available in annotation
* Differential Alternative Splicing (DAS) analysis 
* Differential Transcript Usage (DTU) analysis
* Gene Ontology (GO) enrichment analysis for DEA, DAS and DTU results
* Comprehensive outputs reporting

## Dependencies

CRESCENT is built to work on a Linux system. Many bioinformatic tools it uses are not available at all on Windows, and some libraries are not up to date enough to work on MacOS.

This project uses Conda for package management. All dependencies are listed in `env_crescent.yaml`. To set up the environment, download the yaml and run:

```bash
conda env create -f env_crescent.yaml
conda activate crescent
```

Multiple tools are also deployed via snakemake wrappers as CRESCENT is launched, requiring a working conda installation but no further step by the user than to include '--use-conda' in the snakemake command.
Some older systems might run into trouble with multiqc's latest versions. If multiqc specifically fails when running the pipeline then specifying that version 1.13 should be installed by replacing the line containing it in env_crescent.yaml with the following should solve the problem:
`  - multiqc`

If you intend to use CRESCENT via executors such as an HPC scheduler like slurm or lsf, you need to install the corresponding plugin as well:

	Slurm: `conda install bioconda::snakemake-executor-plugin-slurm`
	LSF: `conda install bioconda::snakemake-executor-plugin-lsf`
	Kubernetes: `conda install bioconda::snakemake-executor-plugin-kubernetes`
	Flux: `conda install bioconda::snakemake-executor-plugin-flux`
	DRMAA: `conda install bioconda::snakemake-executor-plugin-drmaa`

Others can be found in the [snakemake documentation](https://snakemake.github.io/snakemake-plugin-catalog/index.html).

## How to install and prepare

### Dataset and reference genome

In this repository you will find the files required to run the pipeline (snakefile, rule modules, scripts, config file) under snakemake_crescent as well as a small test dataset under Test_data (along with its reference genome, transcriptome, and annotation which are all compressed as .gz and **need to be unzipped before use**).

The config file provided, snakemake_crescent/config/config_crescent.yaml, is set up to work with the test data. Only the **WORKING_DIR and OUTPUT_DIR need to be changed**, providing the directory where both `RNAseq` and `refgen` folders are found as described below. The expected results are provided in Test_data_results.

The dataset needs to be placed in a working folder that you can define in the config file, and which needs to contain a specific structure to work shown below.

Inside an `RNAseq/` folder and either `bam/` if your data is in bam format, or `fastq/` if it is in fastq format, your dataset.

Under `refgen/` your reference genome / transcriptome / annotation.

```
├── RNAseq/
│ ├── bam/  # include .bam files if starting with bam files
│ ├── fastq/ # include .fastq files if starting with fastq files
├── refgen/ # include reference genomes/transcriptomes (fasta) & annotations (.gff or .gtf)
```

After the analysis a new folder will appears in `RNAseq/` called results, though the output can be redirected to another folder than the working folder in which case a new RNAseq/results/ folder will be created there.

Under results will be all result files of the CRESCENT pipeline. See [outputs](#Outputs) for details.
### Config file

The config file under snakemake_crescent/config/ is where all parameters of the analysis are to be set.

The main paramaters to consider before launch are :
- the working directory (pointing to the folder containing `RNAseq/` and `refgen/` described above),
- the output directory where an `RNAseq/results/` folder will be created to contain results. The output directory can be the same as the working directory, or another folder,
- the parameters defining which steps are to be executed (do_rawqc, do_trimming, do_DEA, do_DAS, do_DTU),
- the parameters specific to the dataset being used such as input_mapped, paired_end, stranded, and everything under genome_and_annotation.

Read through all the parameters (except advanced parameters at the end) at least once to be sure you are not missing anything before launching the pipeline.

## How to launch it

Once you have set up the conda environment, have all the pipeline's files, have set up your raw data as required, and have set up your config file, the pipeline is read to be launched.

It is always best to try a dry run first with snakemake pipelines, using -n. As such, first try:
```bash
snakemake -s workflow/Snakefile --use-conda --conda-frontend conda -n
```

If a problem arises, check your config file for potential errors.

When ready to start the analysis, use:
```bash
snakemake -s workflow/Snakefile --use-conda --conda-frontend conda
```

The `--use-conda` option is used to get the wrappers from [snakemake's wrappers repository](https://github.com/snakemake/snakemake-wrappers) for the appropriate steps. The `--conda-frontend conda` option is required to ensure that recent versions of mamba (>2.0) do not conflict with snakemake's calls as they include the now deprecated `--no-default-packages` option.
If using a version of mamba that still has this option available, feel free to remove this part of the above command.
If using conda>=23.10.0 libmamba is not only built into conda but the default solver, so using conda will not slow down dependency resolution.

When using executors such as an HPC scheduler, you will need to use its plugin and specify information such as the number of jobs that can be executed at once or the partition to use. I recommend looking through the documentation of your specific executor [here](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html).

## Outputs

This lists all the outputs produced if the entire workflow is executed.
### Sequencing QC
The sequencing quality control results can be found in `results/sequencing_qc/`.
Html and zip files will be produced for each raw data file provided.
### Trimming
The trimming results can be found in `results/trimmed/`, this includes the trimmed fastq files and the report of the result of trimming produced by the tool used (fastp or cutadapt).
### Mapping and salmon quantification
Mapping results (hisat2 or STAR) will be placed in the `results/mapped/` folder, with text reports on the mapping.

Within `results/mapped/` will also be a `sorted/` folder containing sorted bam files and index bai files.
### Read counting
Read counting via htseq-counts or featureCounts, done for DEA or by itself, will provide outputs in `results/counts/`.
These results include both counts and summaries of the counts.

If DEA is also executed then a `colData.tsv` file containing the match between sample and condition names will be provided, along with a `count_matrix.tsv` file collating the expression data of all the samples.
### DEA outputs
DEA outputs are found in the `results/DEA/` folder, with a folder for edgeR results and a folder for DESeq2 results.

Each folder contains normalized counts, a PCA plot with the two principal components that reflect the most variance in the normalized counts, and a file with all the DEFs (Differentially Expressed Features) containing first gene ids, then the fields provided by the respective tool used.
Most importantly, all contain feature ids, log2FC, and adjusted p-value (under padj in DESeq2 and FDR in edgeR).

Finally a `background` file containing the ids of all features that are at least somewhat expressed in the data (>10 reads for at least one sample) is provided.
### DTU outputs
DTU outputs contain abundance files for the transcripts separated by condition, gene-wise results, and transcript-wise results, all of which are detailed in the output.html file [here](https://github.com/bartongroup/RATS/tree/master/doc) (download  output.html and open it, then click on Abundances, Genes or Transcripts in the summary).

Also found here is `DTU_ids.txt` which lists all the ids of differentially transcribed genes according to DTU analysis.
### DAS outputs
DAS outputs are all found in `results/DAS/` which itself contains a folder of all rMATS outputs, information about which can be found [here](https://rmats.sourceforge.io/user_guide.htm#output).

It also contains the results of filtering on read counts (minimum average of 10 reads across all samples), adjusted p-value and differential in level of inclusion of exon (as defined by the user in config file) on JCEC files and a `DAS_ids.txt` file containing the ids of all genes found differentially spliced in some way (from all types of splicing events) by DAS analysis.
### GO analysis
All GO analysis results are found in `results/GO_analysis`, separated by the type of analysis (DEA, DAS, DTU) and type of GO (BP=Biological Process, MF= Molecular Function, CC=Cellular Compartment). In every case, only significant results are used for GO analysis, dependent on analysis-specific thresholds set by the user in the config file.

For each type of analysis and each type of GO, four output files are produced:

- Bar plots depict enrichment score (p-value) as color and gene count along the X axis for the most significant results.
- Dot plots depict the same information with dot size representing count and the X axis is used for gene ratio (defined as the ratio of input genes that are annotated in a term).
- The cnetplot or gene-concept network depicts the linkages of genes and GO terms as a network.
- The upset plot provided depicts the overlapping versus non-overlapping genes of the most significant GO terms.
### Multiqc
The multiqc results summarizing whichever of sequencing qc, trimming, mapping, salmon quantification and read counting were done will be found in `results/global_report/allqc/multiqc_report.html`.

## Support

You can get support by contacting the project collaborators, in this order :

* gilles.sireta@uca.fr (main developer)
* gwendal.cueff@inrae.fr
* darbotvincent@gmail.com

Corresponding authors (less bioinformatically inclined):
* christophe.tatout@uca.fr
* aline.probst@uca.fr

## Contributing

This project is hosted on gitlab (clean repository to be made). Please feel free to contribute in
any manner, but especially on code optimization, testing and debugging. Additional features could also be added into the pipeline, like additional statistical or functional (biological) analysis. Supplemental visualisations could also improve our tool, at any step.

## Release notes

v1.0
This is a first release of CRESCENT, waiting for new improvements.
 
## License
This project is opensource under a GPLv3 license

