# Read counting is only done if hisat2 or STAR mapping was used
# Counting of reads and assignation to gene / exon / transposable element / etc. via gff annotation file
rule count_reads:
    input:
        bam = OUTPUT_DIR+"RNAseq/results/mapped/namesorted/{filename}.bam" if paired_end else OUTPUT_DIR+"RNAseq/results/mapped/sorted/{filename}.bam",  # mapped and sorted bam file
        bai = OUTPUT_DIR+"RNAseq/results/mapped/sorted/{filename}.bai",  # index for bam files
        annotation = WORKING_DIR+refgen_folder + annotation  # the annotation file to be used (gtf / gff)
    output:
        readcount = OUTPUT_DIR+"RNAseq/results/counts/{filename}_count.txt",  # read count result file
        count_summary = OUTPUT_DIR+"RNAseq/results/counts/{filename}_count.txt.summary" if quantif_tool == "featurecounts" else []  # htseq-count does not output any summary files, as the summary is found a the bottom of the count file
    params:
        featureType = feature,  # exon / transposable_element / mRNA / etc
        idType = id_type,  # gene_id / ID / etc
        stranded = stranded,
        fractional = "-O -M --fraction" if config["DEA"]["overlapping_fractional_counts"]=="yes" and config["DEA"]["multimapped_fractional_counts"]=="yes" else "-O --fraction" if config["DEA"]["overlapping_fractional_counts"]=="yes" else "-M --fraction" if config["DEA"]["multimapped_fractional_counts"]=="yes" else "",
        extra = config["ReadCounting_extra_parameters"]
    threads: threads_per_rule
    priority:1
    benchmark:
        OUTPUT_DIR+"benchmark/count_reads/{filename}.tsv"
    run:
        if quantif_tool == 'featurecounts':
            if not paired_end:
                shell("featureCounts {params.extra} -s {params.stranded} -t {params.featureType} -g {params.idType} -d 20 {params.fractional} -T {threads} -a {input.annotation} -o {output.readcount} {input.bam}")
            else:  # the difference between paired and non-paired in featurecounts is the -p and --countReadPairs parameter : -p states the data is paired-end, and --countReadPairs says to count not individual reads, but the pairs
                shell("featureCounts {params.extra} -s {params.stranded} -t {params.featureType} -p --countReadPairs -g {params.idType} -d 20 {params.fractional} -T {threads} -a {input.annotation} -o {output.readcount} {input.bam}")
        elif quantif_tool == 'htseq':  # htseq-count has no paired-end specific parameters
                shell("htseq-count {params.extra} -s {params.stranded} -i {params.idType} -t {params.featureType} -f bam {input.bam} {input.annotation} > {output.readcount}")

rule count_matrix:
    input:
        expand(OUTPUT_DIR+"RNAseq/results/counts/{filename}_count.txt" if not salmon_map else OUTPUT_DIR+"RNAseq/results/salmon_quant/{filename}/quant.sf", filename = filenames)  # count files
    output:
        count_matrix = DEA_quant_folder + "count_matrix.tsv",  # merged count matrix
        colData = DEA_quant_folder + "colData.tsv",  # samples to conditions match for the DESeq2 implementation
        rds = temp(DEA_quant_folder + "tximp.RDS" if salmon_map else [])
    threads: threads_per_rule
    params:
        OUTPUT_DIR = OUTPUT_DIR,
        path = DEA_quant_folder,
        quantif = quantif_tool,  # method for merging read counts = quantification tool dependant
        sep_file = "_", # Separator of information in file name
        index_sample = 2,  # index (based on separator) of sample in file name
        sample_names = config["samples_and_conditions"]["sample_names"],  # sample names if provided directly
        conditions = config["samples_and_conditions"]["conditions"]  # condition names if provided separately
    benchmark:
        OUTPUT_DIR+"benchmark/count_matrix/count_matrix.tsv"
    script:
        "../scripts/script_count_matrix2.R"

# Building the name of the file which will contain a list of DE elements for use as output of DESeq2 and input of GO_analysis
# First check whether the logFC threshold is an integer
# If so make it type int first so as to not have a .0 at the end of its string
logFC_threshold = config["DEA"]["cutoffs"]["log2FC"]
if logFC_threshold.is_integer():
    logFC_threshold = int(logFC_threshold)
DE_list = "_".join([not_control, 'VS', control, "".join(["DEGstable_padj", str(config["DEA"]["cutoffs"]["pADJ"]),"_logFC",  str(logFC_threshold), ".tsv"])])

rule DESeq2:
    input:
        rds = DEA_quant_folder + "tximp.RDS" if salmon_map else [],
        count_matrix = DEA_quant_folder + "count_matrix.tsv",  # from rule count_matrix
        colData = DEA_quant_folder + "colData.tsv",  # from rule count_matrix
    output:
        norm_read_counts = OUTPUT_DIR + "RNAseq/results/DEA/DESeq2/normalized_counts_matrix.tsv",  # normalized read counts
        DE_list_output = OUTPUT_DIR + "RNAseq/results/DEA/DESeq2/" + DE_list,  # file containing list of DE elements with assorted log2FC
        PCA_plot_pdf = OUTPUT_DIR + "RNAseq/results/DEA/DESeq2/PCAplot.pdf",  # PCA plot pdf
        PCA_plot_png = OUTPUT_DIR + "RNAseq/results/DEA/DESeq2/PCAplot.png",
        background = OUTPUT_DIR + "RNAseq/results/DEA/DESeq2/background"
    threads: threads_per_rule
    params:        
        input_path = DEA_quant_folder,  # count or quantification files path
        output_path = OUTPUT_DIR + "RNAseq/results/DEA/DESeq2/",  # path to place outputs
        quantif = quantif_tool,  # method for merging read counts = quantification tool dependant
        compare = config["samples_and_conditions"]["compare"],  # comparison inputted by user in config file
        padj = config["DEA"]["cutoffs"]["pADJ"],
        log2FC = config["DEA"]["cutoffs"]["log2FC"],
        PCA_labels = config["DEA"]["PCA_labels"]
    benchmark:
        OUTPUT_DIR+"benchmark/DESeq2/DESeq2.tsv"
    script:
        "../scripts/script_DESeq2_2.R"

rule edgeR:
    input:
        rds = DEA_quant_folder + "tximp.RDS" if salmon_map else [],
        count_matrix = DEA_quant_folder + "count_matrix.tsv",  # from rule count_matrix
        colData = DEA_quant_folder + "colData.tsv",  # from rule count_matrix
    output:
        norm_read_counts = OUTPUT_DIR + "RNAseq/results/DEA/edgeR/normalized_counts_matrix.tsv",  # normalized read counts
        DE_list_output = OUTPUT_DIR + "RNAseq/results/DEA/edgeR/" + DE_list,  # file containing list of DE elements with assorted log2FC
        PCA_plot_pdf = OUTPUT_DIR + "RNAseq/results/DEA/edgeR/PCAplot.pdf",  # PCA plot pdf
        PCA_plot_png = OUTPUT_DIR + "RNAseq/results/DEA/edgeR/PCAplot.png",
        background = OUTPUT_DIR + "RNAseq/results/DEA/edgeR/background"
    threads: threads_per_rule
    params:        
        input_path = DEA_quant_folder,  # count or quantification files path
        output_path = OUTPUT_DIR + "RNAseq/results/DEA/edgeR/",  # path to place outputs
        quantif = quantif_tool,  # method for merging read counts = quantification tool dependant
        compare = config["samples_and_conditions"]["compare"],  # comparison inputted by user in config file
        padj = config["DEA"]["cutoffs"]["pADJ"],
        log2FC = config["DEA"]["cutoffs"]["log2FC"],
        PCA_labels = config["DEA"]["PCA_labels"]
    benchmark:
        OUTPUT_DIR+"benchmark/edgeR/edgeR.tsv"
    script:
        "../scripts/script_edgeR.R"

# Setting a variable for the common_DEGs rule's output txt file such that it includes the comparison name
# This allows multiple comparisons to be done without erasing results
common_DEGS_output = OUTPUT_DIR + "RNAseq/results/DEA/" + "common_ids_" + config["samples_and_conditions"]["compare"][0] + "vs" + config["samples_and_conditions"]["compare"][1] + ".txt"
rule common_DEGS:
    input:
        DE_list_edgeR = OUTPUT_DIR + "RNAseq/results/DEA/edgeR/" + DE_list,
        DE_list_DESeq2 = OUTPUT_DIR + "RNAseq/results/DEA/DESeq2/" + DE_list
    output:
        common_DEGs = common_DEGS_output
    shell:
        """
        comm -12 \
            <(tail -n+2 {input.DE_list_edgeR} | cut -f1 | sort) \
            <(tail -n+2 {input.DE_list_DESeq2} | cut -f1 | sort) \
            > {output.common_DEGs}
        """

rule GO_analysis_DEA:
    input:
        gene_list_input = common_DEGS_output if do_DESeq2 and do_edgeR else OUTPUT_DIR + "RNAseq/results/DEA/" + "DESeq2/" + DE_list if do_DESeq2 else OUTPUT_DIR + "RNAseq/results/DEA/" + "edgeR/" + DE_list if do_edgeR else [],  # file containing list of DE elements with assorted log2FC
        background = OUTPUT_DIR + "RNAseq/results/DEA/" + "DESeq2/background" if do_DESeq2 else OUTPUT_DIR + "RNAseq/results/DEA/" + "edgeR/background",
        DESeq2_output = OUTPUT_DIR + "RNAseq/results/DEA/" + "DESeq2/" + DE_list if do_DESeq2 else []
    output:
        barplotgo = OUTPUT_DIR + "RNAseq/results/GO_analysis/DEA/"+config['GO']['ontology_types'][0]+"/barplotGO_" + ("_".join([not_control, "VS", control, "all_genes", config['GO']['ontology_types'][0]]) if config['GO']['separate_up_down'] == "no" else "_".join([not_control, "VS", control, "upregulated", config['GO']['ontology_types'][0]])) + ".pdf"
    threads: threads_per_rule
    params:
        path = OUTPUT_DIR+"RNAseq/results/",
        chosen_species = config['GO']['species'],
        species_index = config['GO']['species_index_among_results'],
        key_type = config['GO']['key_type'],
        ontology_types = config['GO']['ontology_types'],
        compare = config["samples_and_conditions"]["compare"],
        data_type = "DEA",
        min_genes = 10,
        separate_up_down = config['GO']['separate_up_down']
    benchmark:
        OUTPUT_DIR+"benchmark/GO_analysis/DEA_GO.tsv"
    script:
        "../scripts/script_GO.R"
