Phase 2
conda install -c bioconda bedtools -y

# Filter the annotation file to only look at whole genes (ignores individual exons to prevent duplicates)
awk '$3 == "gene"' genes.gff3 > genes_only.gff3

# Find which genes physically overlap with the Gene X binding site
bedtools intersect -a genes_only.gff3 -b gene_x_binding_sites.bed -wa -u > 201950506_target_genes.gff3
