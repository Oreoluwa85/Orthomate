# OrthoMate

This is a collection of scripts in python3 and R, that generate heatmaps corresponding to the numbers of proteins present in [OrthoDB v12](https://www.ezlab.org/orthodb.html). OrthoDB has a mapping function, available with a subscription. This collection of scripts retrieves NCBI taxonomy trees and numbers of proteins per taxa represented in odb12 for a subset of user selected orthogroups.

# Requirements
- Bash
    - grep
- Python3
- R 4.5.0
    - ggtree, ggplot2, treeio, rentrez, tidyr, tidyverse, scales, dplyr, plyr, data.table, UniProt.ws, ggtext, janitor, ggrepel, purrr, patchwork, gt, ggnewscale, stringr, ggExtra, gtExtras, egg, servr, DiagrammeR, webshot2

# Usage
- download using Github
- download odb12 odb12v0_OG2genes.tab
- generate your taxon specific tabfile

> the script generate_odb12_tabfile.sh takes a list of taxon numbers of interest as an input query and a prospective NCBI taxonomy tree filename as an output query of the form: <./generate_odb12_tabfile.sh TaxonLists/querylist.txt querytree>

> the output of the script, if the taxa are in odb12, is a tabfile corresponding to every gene from every taxon represented in odb12 and an NCBI taxonomy tree, corresponding to all taxa in the taxonlist, retrieved by automatically running the script, Get_taxid_newick.py.

 - To generate a linear repertoire of proteins per orthogroup per taxa, run generate_compOG_repertoire.R as follows:
> <Rscript generate_compOG_repertoire.R TabFiles/Tabfile.tab TreeFiles/Treefile.nwk>

- To generate a histogram, describing the number of identified orthogroups per species, run generate_compOG_totals.R as follows:
> <Rscript generate_compOG_totals.R TabFiles/Tabfile.tab TreeFiles/Treefile.nwk>

- To generate an all-against-all, pairwise heatmap comparison of numbers of proteins per pair of compOGs, related by heterodimeric interactions for all taxa of a given list, run generate_compOG_pair_repertoire.R as follows:
> <Rscript generate_compOG_pair_repertoire.R TabFiles/Tabfile.tab>

- To generate an all-against-all, pairwise heatmap comparison of numbers of proteins per pair of compOGs, related by heterodimeric interactions for an individual taxo, run generate_compOG_pair_score.R as follows:
> <Rscript generate_compOG_pair_score.R TabFiles/Tabfile.tab>



