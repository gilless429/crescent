# Defining the inputs for the three 'mapping' possibilities star or hisat2 or salmon depending on : paired_end or not / trimming done or not
if do_trimming:
    map_reads1 = [OUTPUT_DIR+"RNAseq/results/trimmedfq/{filename}_trim.fq"] if not paired_end else [OUTPUT_DIR+"RNAseq/results/trimmedfq/{filename}_1_trim.fq"]
    map_reads2 = [OUTPUT_DIR+"RNAseq/results/trimmedfq/{filename}_2_trim.fq"]
    salmon_reads1 = [OUTPUT_DIR+"RNAseq/results/trimmedfq/{filename}_trim.fq"] if not paired_end else [OUTPUT_DIR+"RNAseq/results/trimmedfq/{filename}_1_trim.fq"]
    salmon_reads2 = [OUTPUT_DIR+"RNAseq/results/trimmedfq/{filename}_2_trim.fq"]
else:
    map_reads1 = [WORKING_DIR+"RNAseq/fastq/{filename}" + fastq_ext] if not paired_end else [WORKING_DIR+"RNAseq/fastq/{filename}_1" + fastq_ext]
    map_reads2 = [WORKING_DIR+"RNAseq/fastq/{filename}_2" + fastq_ext]
    salmon_reads1 = [WORKING_DIR+"RNAseq/fastq/{filename}" + fastq_ext] if not paired_end else [WORKING_DIR+"RNAseq/fastq/{filename}_1" + fastq_ext]
    salmon_reads2 = [WORKING_DIR+"RNAseq/fastq/{filename}_2" + fastq_ext]

## Salmon mapping/quantification is done if salmon is used for DEA or if DTU will be done
# Salmon indexing of reference genome
rule salmon_index:
    input:
        WORKING_DIR+refgen_folder+reftr
    output:
        WORKING_DIR+refgen_folder + "salmon_index_" + config["genome_and_annotation"]["reference_index_transcriptome"] + "/complete_ref_lens.bin"
    params:
        output_directory = WORKING_DIR+refgen_folder + "salmon_index_" + config["genome_and_annotation"]["reference_index_transcriptome"] + "/"
    threads: threads_per_rule
    benchmark:
        OUTPUT_DIR+"benchmark/salmon_index/"+reftr+".tsv"
    run:
        shell("salmon index -p {threads} -t {input} -i {params.output_directory}")

# Create output directory
if (salmon_map or do_DTU) and not os.path.exists(OUTPUT_DIR+"RNAseq/results/salmon_quant/"):
    os.makedirs(OUTPUT_DIR+"RNAseq/results/salmon_quant/")
# Quantify reads via salmon
rule salmon_mapping:
    input:
        fastq = salmon_reads1 if not paired_end else [],  # trimmed single-end output
        fastq1 = salmon_reads1 if paired_end else [],  # trimmed paired-end output 1
        fastq2 = salmon_reads2 if paired_end else [],  # trimmed paired-end output 2
        outindex = WORKING_DIR+refgen_folder + "salmon_index_" + config["genome_and_annotation"]["reference_index_transcriptome"]  + "/complete_ref_lens.bin"
    output:
        quant = OUTPUT_DIR+"RNAseq/results/salmon_quant/{filename}/quant.sf",
        folder = directory(OUTPUT_DIR+"RNAseq/results/salmon_quant/{filename}/")
    params:
        library_type = config["library_type"],
        index_folder = WORKING_DIR+refgen_folder + "salmon_index_" + config["genome_and_annotation"]["reference_index_transcriptome"],
        outdir = OUTPUT_DIR+"RNAseq/results/salmon_quant/{filename}",
        bootstraps = config["bootstraps"],
        extra = config["Mapping_extra_parameters"]
    priority: 1
    benchmark:
        OUTPUT_DIR+"benchmark/salmon_mapping/{filename}.tsv"
    run:
        if paired_end:
            shell("salmon quant -l A --index {params.index_folder} -1 {input.fastq1} -2 {input.fastq2} -o {params.outdir} --numBootstraps {params.bootstraps}")
        else:
            shell("salmon quant -l A --index {params.index_folder} -r {input.fastq} -o {params.outdir} --numBootstraps {params.bootstraps}")

# Classic mapping (hisat2/STAR) is done if a tool other than salmon is specified for DEA or if DAS will be done
# Reference genome index building
rule index_refgen:
    input:
        fasta = WORKING_DIR+refgen_folder+refgen  # reference genome fasta file
    output:
        directory(WORKING_DIR+refgen_folder +"index_"+reference_index)  # reference genome output folder name, using the reference index to separate various reference genomes' indexes
    params:
        prefix = WORKING_DIR+refgen_folder+"index_"+reference_index+"/"+reference_index,  # only used if the indexer is hisat2 index, provides the prefix such that files will be named prefix.1.ht2, prefix.2.ht2, etc...
        extra = "--sjdbOverhang " + str(indexation_sjdbOverhang) + " --sjdbGTFfile " + WORKING_DIR+refgen_folder+annotation + " --genomeSAindexNbases 12" if "star" in map_tool else []  # Only gets used if the indexer is STAR index
    log:
        OUTPUT_DIR + refgen_folder + "index_"+reference_index+"/logs/hisat2_index_"+refgen+".log" if "hisat2" in map_tool else OUTPUT_DIR + refgen_folder+"index_"+reference_index+"/logs/star_index_"+refgen+".log"
    threads: config["threads_genomeindex"]
    resources: mem_mb = config["memory_genomeindex"]
    priority:1
    benchmark:
        OUTPUT_DIR+"benchmark/index_refgen/"+refgen+".tsv"
    wrapper:
        index_wrapper

# Mapping reads to reference genome, changing output to bam file immediately, and removing trimmed fq files once they are unnecessary
rule mapping:
    input:
        fq1 = map_reads1,
        fq2 = map_reads2 if paired_end else [],
        idx = WORKING_DIR+refgen_folder+"index_"+reference_index
    output:
        aln = temp(OUTPUT_DIR+'RNAseq/results/mapped/{filename}_mapped.bam' if "hisat2" in map_tool else OUTPUT_DIR+"RNAseq/results/mapped/{filename}_Aligned.out.bam"),
        log_final = OUTPUT_DIR+'RNAseq/results/mapped/{filename}_Log.final.out' if "star" in map_tool else OUTPUT_DIR+"RNAseq/results/mapped/{filename}.txt"
    log:
        OUTPUT_DIR+"RNAseq/results/mapped/{filename}.log" if "hisat2" in map_tool else OUTPUT_DIR+"RNAseq/results/mapped/{filename}.log"
    params:
        idx = WORKING_DIR+refgen_folder+"index_"+reference_index+"/"+reference_index if "hisat2" in map_tool else WORKING_DIR+refgen_folder+"index_"+reference_index,  # reference genome index, hisat2 requires the prefix of the files where star requires the folder
        extra = (('--alignIntronMin '+alignment_intron_min+' --alignIntronMax '+alignment_intron_max+' --outSAMprimaryFlag AllBestScore --outSAMtype BAM Unsorted --outStd BAM_Unsorted') if "star" in map_tool else "") + config["Mapping_extra_parameters"] 
    threads: threads_per_rule
    priority:1
    benchmark:
        OUTPUT_DIR+"benchmark/mapping/{filename}.tsv"
    run:
        if map_tool == "hisat2":
            if not paired_end:
                shell("hisat2 -p {threads} {params.extra} -x {params.idx} -U {input.fq1} --new-summary --summary-file {output.log_final} | samtools view -bSh - > {output.aln}")
            else:
                shell("hisat2 -p {threads} {params.extra} -x {params.idx} -1 {input.fq1} -2 {input.fq2} --new-summary --summary-file {output.log_final} | samtools view -bSh - > {output.aln}")
        elif map_tool == "star":
            if not paired_end:
                shell("STAR --runMode alignReads --alignEndsType Local {params.extra} --genomeDir {params.idx} --readFilesIn {input.fq1} > {output.aln} {output.log_final}")
            else:
                shell("STAR --runMode alignReads --alignEndsType Local {params.extra} --genomeDir {params.idx} --readFilesIn {input.fq1} {input.fq2} > {output.aln} {output.log_final}")



# The input of the sorting step depends upon whether the input to the pipeline is mapped already (in which case we go right into this step without any trimming/mapping)
# But also upon the mapper used if the mapping is done by the pipeline instead
sorting_input = OUTPUT_DIR+"RNAseq/bam/{filename}.bam" if input_mapped else (OUTPUT_DIR+"RNAseq/results/mapped/{filename}_mapped.bam" if "hisat2" in map_tool else OUTPUT_DIR+"RNAseq/results/mapped/{filename}_Aligned.out.bam") 

# Sorting bam by coordinates and removing original mapped bam ; coordinate-sorted bams are better for visualization among others
rule sort:
    input:
        sorting_input  # Output name-sorted bam file
    output:
        OUTPUT_DIR+"RNAseq/results/mapped/sorted/{filename}.bam"  # Output coordinate-sorted bam file
    params:
        extra = config["BamSorting_extra_parameters"]
    threads: threads_per_rule
    priority:1
    benchmark:
        OUTPUT_DIR+"benchmark/sort/{filename}.tsv"
    run:
        if sort_tool == "samtools":
            shell("samtools sort {input} -o {output} -@ {threads} {params.extra}")
        elif sort_tool == "sambamba":
            shell("sambamba sort -o {output} -t {threads} {params.extra} {input}")


# Sorting by names : if in paired-end, it's better for featurecounts (which namesorts itself if necessary) and htseq-count (which doesn't but handles other sorting poorly at times)
# The name-sorted bams are dropped once the quantification has occurred as they are not useful beyond that whereas position-sorted bams can be used for visualization
# Position-sorted bam are used as input so that if the workflow is done again without count files, it is unnecessary to re-execute mapping to get namesorted bams again.
rule sortname:
    input:
        OUTPUT_DIR+"RNAseq/results/mapped/sorted/{filename}.bam"
    output:
        OUTPUT_DIR+"RNAseq/results/mapped/namesorted/{filename}.bam"  # Output name-sorted bam file
    params:
        extra = config["BamSorting_extra_parameters"]
    threads: threads_per_rule
    priority:1
    benchmark:
        OUTPUT_DIR+"benchmark/sortname/{filename}.tsv"
    run:
        # -n is the required option to sort by name instead of coordinates
        if sort_tool == "samtools":
            shell("samtools sort -n {input} -o {output} -@ {threads} {params.extra}")
        elif sort_tool == "sambamba":
            shell("sambamba sort -n -o {output} -t {threads} {params.extra} {input}")

# Indexing sorted bam files
rule index_bam:
    input:
        OUTPUT_DIR+"RNAseq/results/mapped/sorted/{filename}.bam"  # Input sorted bam file
    output:
        OUTPUT_DIR+"RNAseq/results/mapped/sorted/{filename}.bai"  # Output bai index file of the above bam file
    params:
        extra = config["BamIndexing_extra_parameters"]
    threads: threads_per_rule
    priority:1
    benchmark:
        OUTPUT_DIR+"benchmark/index_bam/{filename}.tsv"
    wrapper:
        bamindex_wrapper
