#######################################################
# Verification of sufficient gene number for analysis #
#######################################################

# Input and expected output
genes_file <- snakemake@input[["gene_list_input"]]
barplot_pdf <- snakemake@output[["barplotgo"]]
dotplot_pdf <- snakemake@output[["dotplotgo"]]
cnetplot_pdf <- snakemake@output[["cnetplotgo"]]
upsetplot_pdf <- snakemake@output[["upsetplotgo"]]

# Minimum number of genes needed to run GO analysis (default 10)
min_genes <- as.integer(snakemake@params[["min_genes"]])
# Count the number of genes in the input file
num_genes <- length(readLines(genes_file)) - 1

# If there aren't enough genes a message is displayed and placeholder pdfs created
if (num_genes < min_genes) {
  message(paste0("GO analysis skipped for ", snakemake@params[['data_type']], ": Not enough genes (minimum ", min_genes, " required)."))
  # Create a placeholder PDF
  pdf(barplot_pdf, width=8, height=6)
  plot.new()
  text(0.5, 0.5, paste("GO analysis skipped:", num_genes, "genes found, but", min_genes, "required"), cex=1.5)
  dev.off()
  pdf(dotplot_pdf, width=8, height=6)
  plot.new()
  text(0.5, 0.5, paste("GO analysis skipped:", num_genes, "genes found, but", min_genes, "required"), cex=1.5)
  dev.off()
  pdf(cnetplot_pdf, width=8, height=6)
  plot.new()
  text(0.5, 0.5, paste("GO analysis skipped:", num_genes, "genes found, but", min_genes, "required"), cex=1.5)
  dev.off()
  pdf(upsetplot_pdf, width=8, height=6)
  plot.new()
  text(0.5, 0.5, paste("GO analysis skipped:", num_genes, "genes found, but", min_genes, "required"), cex=1.5)
  dev.off()
} else {
  # Else we run the GO analysis
  library(tidyverse, quietly = T)
  library(clusterProfiler, quietly = T)
  library(enrichplot, quietly = T)

  ###############################
  # Database and key type setup #
  ###############################
  # All available orgdb pre-builts (as of 03/10/2024)
  orglist <- c("org.Hs.eg.db",      # Homo sapiens
              "org.Mm.eg.db",      # Mus Musculus
              "org.Ce.eg.db",      # Caenorhabditis elegans
              "org.Sc.sgd.db",     # Saccharomyces cerevisiae
              "org.Dm.eg.db",      # Drosophila melanogaster
              "org.At.tair.db",    # Arabidopsis thaliana
              "org.Dr.eg.db",      # Danio rerio
              "org.Ag.eg.db",      # Anopheles gambiae
              "org.Bt.eg.db",      # Bos taurus
              "org.Cf.eg.db",      # Canis familiaris
              "org.EcK12.eg.db",   # Escherichia coli strain K12
              "org.EcSakai.eg.db", # Escherichia coli strain Sakai
              "org.Gg.eg.db",      # Gallus gallus
              "org.Mmu.eg.db",     # Macaca mulatta
              "org.Mxanthus.db",   # Myxococcus xanthus
              "org.Pf.plasmo.db",  # Plasmodium falciparum
              "org.Pt.eg.db",      # Pan troglodytes
              "org.Rn.eg.db",      # Rattus norvegicus
              "org.Ss.eg.db",      # Sus scrofa
              "org.Xl.eg.db"       # Xenopus laevis
              )
  # All available orgdb pre-builts --> species names as list names for packages
  names(orglist) <- c("Homo sapiens",
                      "Mus musculus",
                      "Caenorhabditis elegans",
                      "Saccharomyces cerevisiae",
                      "Drosophila melanogaster",
                      "Arabidopsis thaliana",
                      "Danio rerio",
                      "Anopheles gambiae",
                      "Bos taurus",
                      "Canis familiaris",
                      "Escherichia coli strain K12",
                      "Escherichia coli strain Sakai",
                      "Gallus gallus",
                      "Macaca mulatta",
                      "Myxococcus xanthus",
                      "Plasmodium falciparum",
                      "Pan troglodytes",
                      "Rattus norvegicus",
                      "Sus scrofa",
                      "Xenopus laevis") 

  # Species chosen by user
  chosen_species <- snakemake@params[["chosen_species"]]
  # If the chosen species is among those with a pre-built package, then use that
  if (chosen_species %in% names(orglist)) {
    chosen_annotation <- orglist[[chosen_species]]
    if (!require(chosen_annotation, character.only = T, quietly = T)) {
      BiocManager::install(chosen_annotation, update = F)
      invisible(library(chosen_annotation, character.only = T, quietly = T))
    }
    org <- eval(as.name(chosen_annotation))
    } else {  # If the chosen species does not have a pre-built package then draw from AnnotationHub's list of orgdb
      if (!require('AnnotationHub', character.only = T)) {
        BiocManager::install('AnnotationHub')
        library(AnnotationHub, quietly = T)
      }
      ah <- AnnotationHub::AnnotationHub()
      query_result <- query(ah, chosen_species)
      if (length(query_result) > 1) {
        print(query_result)
        cat("Multiple results for the specified species.")
        cat("Which of the above results would you like to keep (give the number in order) ?")
        if (snakemake@params[["species_index"]] != "") {
          index_org = snakemake@params[["species_index"]]
        } else {
          index_org <- readLines("stdin", n = 1)
        }
        org <- query_result[[as.numeric(index_org)]]
        } else if (length(query_result) == 0) {
          stop("No species in the list of available species matches your input. Check that it is spelled correctly (in latin).")
          } else {
            org <- query_result
          }
    }

  # Ask user for key type to use
  if (snakemake@params[["key_type"]] != "") {
      key_type = snakemake@params[["key_type"]]
      } else {
        print(keytypes(org))
        cat("Multiple key types are available in the annotation of this species.")
        not_chosen <- T
        while (not_chosen) {
          cat("Type or paste here whichever of the above matches your annotation's gene IDs (key type + ? shows examples, ex: REFSEQ?):")
          
            key_type <- readLines("stdin", n = 1)
          }
          if (grepl("?", key_type, fixed = T)) {
            cat(head(keys(org, gsub("?", "", key_type, fixed = T))))
          } else {
            not_chosen <- F
          }
        }

  #######################
  # Enrichment analysis #
  #######################

  # Grab table with gene_ids, then the gene ids in the first column
  gene_ids_table <- read_delim(genes_file, delim = "\t")
  gene_ids <- gene_ids_table[[1]]
  
  # Separation of positive VS negative logFC / IncLevelDiff / 
  if (snakemake@params[["data_type"]] == "DEA" && snakemake@params[["separate_up_down"]] == "yes") {
    if ("padj" %in% colnames(gene_ids_table)) {  # DESeq2 results
      up_gene_ids <- gene_ids[gene_ids_table$log2FoldChange > 0]
      down_gene_ids <- gene_ids[gene_ids_table$log2FoldChange < 0]
    } else if ("FDR" %in% colnames(gene_ids_table)) {  # edgeR results
      up_gene_ids <- gene_ids[gene_ids_table$logFC > 0]
      down_gene_ids <- gene_ids[gene_ids_table$logFC < 0]
    } else {  # both edgeR and DESeq2 results used, DESeq2 results used to define if gene up/down-regulated
      deseq2_results <- read_delim(snakemake@input[["DESeq2_output"]], delim="\t")
      # Filter deseq2_results to keep only genes in gene_ids
      filtered_results <- deseq2_results[deseq2_results[[1]] %in% gene_ids, ]
      # Create up and down gene lists
      up_gene_ids <- filtered_results[[1]][filtered_results$log2FoldChange > 0]
      down_gene_ids <- filtered_results[[1]][filtered_results$log2FoldChange < 0]
    }
  }

  no_background <- F
  # Grab background (all genes expressed at least somewhat in the dataset)
  if (is.null(snakemake@input[["background"]])) {
    no_background <- T
  } else { background <- scan(snakemake@input[["background"]], character()) }

  # Get compare value to create plots with names matching them
  compare <- snakemake@params[["compare"]]
  not_control <- compare[1]
  control <- compare[2]

  # Determine which gene list(s) to use for GO analysis
  if (exists("up_gene_ids") && exists("down_gene_ids") && (length(up_gene_ids) > 0 | length(down_gene_ids) > 0)) {
    gene_sets <- list(
      upregulated = up_gene_ids,
      downregulated = down_gene_ids
    )
  } else {
    gene_sets <- list(all_genes = gene_ids)
    }

  # GO analysis
  ontology_types <- snakemake@params[["ontology_types"]]
  for (ontology_type in ontology_types) {
    for (set_name in names(gene_sets)) {
      current_genes <- gene_sets[[set_name]]
      
      if (length(current_genes) == 0) next  # Skip if the list is empty

      # Create and go to appropriate directory to put results in
      results_path <- file.path(snakemake@params[["path"]], "GO_analysis", snakemake@params[['data_type']], ontology_type)
      dir.create(results_path, showWarnings = FALSE, recursive = TRUE)
      # Enrichment analysis
      if (no_background) {
        ego <- try(enrichGO(gene = current_genes,
                            OrgDb = org,
                            ont = ontology_type,
                            keyType = key_type,
                            pvalueCutoff = 1,
                            qvalueCutoff = 1), silent = TRUE)
      } else {
        ego <- try(enrichGO(gene = current_genes,
                            OrgDb = org,
                            universe = background,
                            ont = ontology_type,
                            keyType = key_type,
                            pvalueCutoff = 1,
                            qvalueCutoff = 1), silent = TRUE)
      }
      if (inherits(ego, "try-error") || is.null(ego)) {
        stop("There was a problem with the GO analysis. Perhaps the key type used in your annotation is not the same as that specified for use here, or perhaps there were too few valid elements for GO analysis.")
      }
      # Simplify and generate plots
      simplego <- clusterProfiler::simplify(ego)
      # Create base file name to be reused for all plots
      if
      base_filename <- paste(not_control, "VS", control, set_name, ontology_type, sep = "_")
      # Output plots
      barplot(simplego)
      ggsave(file.path(results_path, paste0("barplotGO_", base_filename, ".pdf")))
      dotplot(simplego)
      ggsave(file.path(results_path, paste0("dotplotGO_", base_filename, ".pdf")))
      cnetplot(simplego)
      ggsave(file.path(results_path, paste0("cnetplotGO_", base_filename, ".pdf")))
      enrichplot::upsetplot(simplego)
      ggsave(file.path(results_path, paste0("upsetplotGO_", base_filename, ".pdf")))
    }
  }
}