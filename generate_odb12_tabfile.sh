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
    LC_ALL=C "${cmd[@]}" odb12v0_OG2genes.tab > "TabFiles/${2}.tab"
    echo "done"
fi
}

echo "building treefile"
python3 Get_taxid_newick.py $1 $2.nwk
echo "Done. Ready to run Orthosearch"

echo "type Rscript generate_compOG_repertoire.R TabFiles/$(basename "$1" .txt).tab TreeFiles/$2.nwk to generate compOG repertoire"

echo "type Rscript generate_compOG_totals.R TabFiles/$(basename "$1" .txt).tab TreeFiles/$2.nwk to generate histogram and NCBI taxonomy tree, describing numbers of orthogroups identified per species from the submitted TaxonList"

echo "type Rscript generate_compOG_pair_repertoire.R TabFiles/$(basename "$1" .txt).tab to generate an all-against-all heatmap corresponding to the numbers of proteins shared between orthogroups linked by proteins that heterodimerically interact, summed for all taxa in the submitted TaxonList"

echo "type Rscript generate_compOG_pair_score.R TabFiles/$(basename "$1" .txt).tab to generate an all-against-all heatmap corresponding to the numbers of proteins shared between orthogroups linked by proteins that heterodimerically interact, for a particular taxon from the submitted TaxonList, compared to all other taxa"