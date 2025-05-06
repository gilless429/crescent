import glob  # imported to find bam files to do rMATs on and put their names into files to give it


# Create the two files detailing which groups the sorted bam belong to
rule file_condition_match:
    input:
        expand(OUTPUT_DIR+"RNAseq/results/mapped/sorted/{filename}.bam", filename=filenames)  # checks that all of the sorting and indexing is done post-mapping
    output:
        b1 = OUTPUT_DIR+"RNAseq/results/DAS/b1.txt",
        b2 = OUTPUT_DIR+"RNAseq/results/DAS/b2.txt"
    benchmark:
        OUTPUT_DIR+"benchmark/file_condition_match/file_condition_match.tsv"
    run:
        # The two parameters below use list comprehension to grab the list of files matching each condition by sequentially comparing both of the conditions to use (compare parameter in config file) with the ordered list of conditions in the conditions parameter
        # glob is required, thus its importing above
        rMATS_files1 = [cond1_bam for index1, cond1_bam in enumerate(sorted(glob.glob(OUTPUT_DIR+"RNAseq/results/mapped/sorted/*bam"))) if config['samples_and_conditions']['conditions'][index1] == config['samples_and_conditions']['compare'][0]]
        rMATS_files2 = [cond2_bam for index2, cond2_bam in enumerate(sorted(glob.glob(OUTPUT_DIR+"RNAseq/results/mapped/sorted/*bam"))) if config['samples_and_conditions']['conditions'][index2] == config['samples_and_conditions']['compare'][1]]
        with open(output.b1,'w') as f1:
            f1.write(",".join(rMATS_files1))
        with open(output.b2,'w') as f2:
            f2.write(",".join(rMATS_files2))

# Do the DAS analysis
rule rMATS:
    input:
        annotation = WORKING_DIR+refgen_folder+annotation_DTU,  # the annotation file to be used (gtf)
        mapped = expand(OUTPUT_DIR+"RNAseq/results/mapped/sorted/{filename}.bam", filename=filenames),
        b1 = OUTPUT_DIR+"RNAseq/results/DAS/b1.txt",
        b2 = OUTPUT_DIR+"RNAseq/results/DAS/b2.txt"
    output:
        OUTPUT_DIR+"RNAseq/results/DAS/rMATS_output/SE.MATS.JCEC.txt"  # testing for completion via one of the output files
    params:
        read_length = config["DAS"]["read_length"],
        bam_dir = OUTPUT_DIR+"RNAseq/results/mapped/sorted/",
        output_dir = OUTPUT_DIR+"RNAseq/results/DAS/rMATS_output/",
        tmp_output_dir = OUTPUT_DIR+"RNAseq/results/DAS/tmp/",
        stranded = 'fr-firststrand' if stranded == 1 or stranded == "yes" else 'fr-unstranded' if stranded == 0 or stranded == "no" else 'fr-secondstrand' if stranded == 2 or stranded == "reverse" else [],
        paired = 'paired' if paired_end else 'single',
        variable_length = '--variable-read-length' if config['DAS']['variable_read_length'] == 'yes' else []
    threads: threads_per_rule
    benchmark:
        OUTPUT_DIR+"benchmark/rMATS/rMATS.tsv"
    run:
        shell("rmats.py --b1 {input.b1} --b2 {input.b2} --gtf {input.annotation} -t {params.paired}  --libType {params.stranded} --readLength {params.read_length} {params.variable_length} --nthread {threads} --od {params.output_dir} --tmp {params.tmp_output_dir}")

# Filter the results
rule rMATS_sig:
    input:
        OUTPUT_DIR+"RNAseq/results/DAS/rMATS_output/SE.MATS.JCEC.txt"
    output:
        OUTPUT_DIR+"RNAseq/results/DAS/DAS_ids.txt"
    params:
        rmats_output = OUTPUT_DIR+"RNAseq/results/DAS/rMATS_output/",
        script_output_folder = OUTPUT_DIR+"RNAseq/results/DAS/",
        padj = config["DAS"]["pADJ"],
        psi = config["DAS"]["cutoff"]
    benchmark:
        OUTPUT_DIR+"benchmark/rMATS_sig/rMATS_sig.tsv"
    script:
        "../scripts/rMATS_filt.py"

rule GO_analysis_DAS:
    input:
        gene_list_input = OUTPUT_DIR+"RNAseq/results/DAS/DAS_ids.txt"
    output:
        barplotgo = OUTPUT_DIR + "RNAseq/results/GO_analysis/DAS/"+config['GO']['ontology_types'][0]+"/barplotGO_" + "_".join([not_control, "VS", control, "all_genes", config['GO']['ontology_types'][0]]) + ".pdf",
    threads: threads_per_rule
    params:
        path = OUTPUT_DIR+"RNAseq/results/",
        chosen_species = config['GO']['species'],
        species_index = config['GO']['species_index_among_results'],
        key_type = config['GO']['key_type'],
        ontology_types = config['GO']['ontology_types'],
        compare = config["samples_and_conditions"]["compare"],
        data_type = "DAS",
        min_genes = 10
    benchmark:
        OUTPUT_DIR+"benchmark/GO_analysis/DAS_GO.tsv"
    script:
        "../scripts/script_GO.R"