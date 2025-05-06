
# Grab the list of all conditions, in order
conditions_list = config['samples_and_conditions']['conditions']

rule run_RATS:
    input:
        expand(OUTPUT_DIR+"RNAseq/results/salmon_quant/{filename}/quant.sf", filename = filenames)
    output:
        DTU_ids_tsv = OUTPUT_DIR+"RNAseq/results/DTU/DTU_ids.txt",
        DTU_genes_tsv = OUTPUT_DIR+"RNAseq/results/DTU/DTU_results_gene_level.tsv",
        DTU_transcripts_tsv = OUTPUT_DIR+"RNAseq/results/DTU/DTU_results_transcript_level.tsv",
        abundance_cond1 = OUTPUT_DIR+'RNAseq/results/DTU/abundance_'+config['samples_and_conditions']['compare'][0]+'.txt',
        abundance_cond2 = OUTPUT_DIR+'RNAseq/results/DTU/abundance_'+config['samples_and_conditions']['compare'][1]+'.txt'
    params:
        quant_folder = OUTPUT_DIR+"RNAseq/results/salmon_quant/",
        # Provide the list of conditions in order
        conditions = config['samples_and_conditions']['conditions'],
        # Grab the two comparison conditions in their own variable each to 
        comparison_condition_1 = config['samples_and_conditions']['compare'][0],
        comparison_condition_2 = config['samples_and_conditions']['compare'][1],
        padj = config["DTU"]["pADJ"],
        cutoff = config["DTU"]["cutoff"],
        annotation = OUTPUT_DIR + refgen_folder + annotation_DTU
    benchmark:
        OUTPUT_DIR+"benchmark/run_RATS/run_RATS.tsv"
    script:
        "../scripts/DTU.R"

rule GO_analysis_DTU:
    input:
        gene_list_input = OUTPUT_DIR+"RNAseq/results/DTU/DTU_ids.txt"
    output:
        barplotgo = OUTPUT_DIR + "RNAseq/results/GO_analysis/DTU/"+config['GO']['ontology_types'][0]+"/barplotGO_" + "_".join([not_control, "VS", control, "all_genes", config['GO']['ontology_types'][0]]) + ".pdf",
    threads: threads_per_rule
    params:
        path = OUTPUT_DIR+"RNAseq/results/",
        chosen_species = config['GO']['species'],
        species_index = config['GO']['species_index_among_results'],
        key_type = config['GO']['key_type'],
        ontology_types = config['GO']['ontology_types'],
        compare = config["samples_and_conditions"]["compare"],
        data_type = "DTU",
        min_genes = 10
    benchmark:
        OUTPUT_DIR+"benchmark/GO_analysis/DTU_GO.tsv"
    script:
        "../scripts/script_GO.R"