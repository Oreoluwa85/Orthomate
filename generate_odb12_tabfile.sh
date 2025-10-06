#!/bin/bash
clear
echo 'usage: generate_odb12_tabfile.sh <TaxonLists/TaxonFile.txt> <OutputTabFile.tab>'
read -p "Press enter to continue"

{
if [ ! -f "TabFiles/$2.tab" ]; then
    echo "building tabfile"
    # build grep command without LC_ALL
    cmd=(grep -F -i)
    while IFS= read -r p; do
      cmd+=("-e" $'\t'"${p}_")
    done < "$1"
    # print the full command as it will run
    printf 'LC_ALL=C '
    printf '%q ' "${cmd[@]}"
    printf 'odb12v0_OG2genes.tab > %q\n' "TabFiles/${2}.tab"
    # execute the command with LC_ALL set
    LC_ALL=C "${cmd[@]}" odb12v0_OG2genes.tab > TabFiles/"${2}.tab"
    echo "done"
fi
}

echo "building treefile"
python3 Get_taxid_newick.py $1 TreeFiles/$2.nwk
echo "Done. Ready to run Orthosearch"

echo "type Rscript generate_compOG_repertoire.R TabFiles/$(basename "$1" .txt).tab TreeFiles/$2.nwk to generate compOG repertoire"
echo "type Rscript generate_compOG_totals.R TabFiles/$(basename "$1" .txt).tab TreeFiles/$2.nwk to generate histogram of total orthologous groups identified per taxon"
echo "type Rscript generate_compOG_pair_repertoire.R TabFiles/$(basename "$1" .txt).tab to generate an all-against-all, pairwise heatmap comparison of numbers of proteins across a group of taxa"
echo "type Rscript generate_compOG_pair_score.R TabFiles/$(basename "$1" .txt).tab taxid to generate an all-against-all, pairwise heatmap comparison of numbers of proteins per pair of compOGs related by heterodimeric interactions, using any taxid from the original TaxonList file to select a taxon of interest"
