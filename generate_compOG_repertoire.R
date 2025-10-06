library(ggtree)
library(ggplot2)
library(treeio)
library(rentrez)
library(tidyr)
library(tidyverse)
library(scales)
library(dplyr)
library(data.table)
library(UniProt.ws)
library(ggtext)
library(janitor)
library(ggrepel)
library(purrr)
library(patchwork)
library(gt)
library(ggnewscale)
library(stringr)
library(ggExtra)
library(gtExtras)
library(egg)
library(servr)
library(DiagrammeR)
library(webshot2)
args <- commandArgs(TRUE)
i <- read.table(args[1], sep="\t", header=F)
j <- read.nhx(args[2])
filename <- gsub("TabFiles/","", args[1])

# i <- read.table("/Users/of2/Documents/Complex_Portal/ComPred/Orthomate/TabFiles/microsporidia.tab", sep="\t", header=F)
# j <- read.nhx("/Users/of2/Documents/Complex_Portal/ComPred/Orthomate/TreeFiles/microsporidia.nwk")
# setwd("/Users/of2/Documents/Complex_Portal/ComPred/Orthomate/")

Uniprot_annotated_complexes<- as.data.frame(read_csv("Resources/Uniprot_annotated_complexes.csv"))
Complexes2 <- read_csv("Resources/Complexes2.csv")
Uniprot_annotated_complexes$OrthoDB <- gsub('.{1}$', '', Uniprot_annotated_complexes$OrthoDB)



compOGrepertoire <- c("1937912at2759", "413460at2759", "10261556at2759", "5575062at2759", 
                      "10249535at2759", "10060499at2759", "5132116at2759", "308383at2759", 
                      "1918432at2759", "427795at2759", "5857104at2759", "1608002at2759", 
                      "1926878at2759", "21449at2759", "534348at2759", "3763at2759", 
                      "613763at2759", "2414538at2759", "336885at2759", "29755at2759", 
                      "4199794at2759", "3000483at2759", "430051at2759", "446168at2759", 
                      "10250002at2759", "419317at2759", "206088at2759", "2143914at2759", 
                      "10258825at2759", "272141at2759", "5800476at2759", "1732493at2759", 
                      "7851174at2759", "448448at2759", "8962942at2759", "2415936at2759", 
                      "1716531at2759", "9984419at2759", "6600758at2759", "411372at2759", 
                      "9332038at2759", "193931at2759", "63267at2759", "10254671at2759", 
                      "537915at2759", "417078at2759", "108365at2759", "10262475at2759", 
                      "10264728at2759", "340346at2759", "1930084at2759", "16120at2759", 
                      "10252509at2759", "10255768at2759", "10256771at2759", "1093at2759", 
                      "268428at2759", "268763at2759", "429533at2759", "431557at2759", 
                      "2187549at2759", "2877at2759", "1183224at2759", "1933107at2759", 
                      "412748at2759", "64353at2759", "372487at2759", "10265243at2759", 
                      "28737at2759", "193499at2759", "10248812at2759", "10253204at2759", 
                      "6513042at2759", "10256122at2759", "248779at2759", "10248617at2759", 
                      "270173at2759", "282152at2759", "270392at2759", "10253254at2759", 
                      "289038at2759", "64767at2759", "10251154at2759", "269804at2759", 
                      "273340at2759", "441223at2759", "511529at2759", "2127950at2759", 
                      "10262986at2759", "272481at2759", "9991317at2759", "8068875at2759", 
                      "431715at2759", "3142434at2759", "6780543at2759", "10250970at2759", 
                      "10260794at2759", "196131at2759", "364892at2759", "5575at2759", 
                      "258143at2759", "7875889at2759", "1882297at2759", "428577at2759", 
                      "2401965at2759", "342024at2759", "10255414at2759", "1685042at2759", 
                      "424572at2759", "24966at2759", "10265785at2759", "590761at2759", 
                      "10254527at2759", "2110130at2759", "68056at2759", "1706657at2759", 
                      "197206at2759", "21243at2759", "10250478at2759", "10264910at2759", 
                      "5571054at2759", "1698572at2759", "238316at2759", "2250022at2759", 
                      "27435at2759", "422728at2759", "6500128at2759", "28053at2759", 
                      "10256289at2759", "49016at2759", "542917at2759", "660555at2759", 
                      "10264038at2759", "1650at2759", "372421at2759", "5061070at2759", 
                      "192611at2759", "2423701at2759", "3176171at2759", "10036721at2759", 
                      "10251574at2759", "1882346at2759", "844at2759", "6108017at2759", 
                      "6353017at2759", "2357150at2759", "256303at2759", "340608at2759", 
                      "377733at2759", "248923at2759", "5839at2759", "294251at2759", 
                      "10250354at2759", "550424at2759", "10255013at2759", "1719357at2759", 
                      "27923at2759", "405996at2759", "10250817at2759", "1727884at2759", 
                      "10248520at2759", "1748577at2759", "240298at2759", "10262857at2759", 
                      "246406at2759", "1662883at2759", "10263554at2759", "10264220at2759", 
                      "2011769at2759", "1706066at2759", "10253113at2759", "284782at2759", 
                      "26525at2759", "10260625at2759", "364224at2759", "415696at2759", 
                      "1649088at2759", "10266385at2759", "10250117at2759", "1741334at2759", 
                      "1856718at2759", "427280at2759", "407658at2759", "10265628at2759"
)

  colnames(i) <- c("OrthoDB","OrthoDB_genes")
  # ChatGPT pivot solution
  
  # Extract taxon ID (everything before the first underscore)
  query_orthodb_stripped <- i %>%
    mutate(TaxonID = sub("_.*", "", OrthoDB_genes))
  
  # Restrict OrthoDB_stripped to eukaryote-wide only
  #query_orthodb_stripped <- query_orthodb_stripped[query_orthodb_stripped$OrthoDB %like% "at2759", ]
  query_orthodb_stripped <- query_orthodb_stripped[query_orthodb_stripped$OrthoDB %in% compOGrepertoire, ] 
  
  # Count occurrences of each taxon ID per OrthoDB - alternative below
  query_orthodb_pivot <- query_orthodb_stripped %>%
    group_by(OrthoDB, TaxonID) %>%
    summarise(Count = n(), .groups = "drop") %>%
    pivot_wider(names_from = TaxonID, values_from = Count, values_fill = 0)
  
  ## Merge?
  query_orthodb_merge <- merge(Uniprot_annotated_complexes,query_orthodb_pivot, .keep_all=T)
  query_orthodb_merge <- query_orthodb_merge %>% distinct(OrthoDB, .keep_all = T)
  
  ## Pivot merged table to columns per taxa
  orthodb_merge_pivot <- as.tibble(query_orthodb_merge[c(1,7:ncol(query_orthodb_merge))]) %>% pivot_longer(cols=!OrthoDB, names_to = "taxid", values_to = "count") %>% pivot_wider(names_from = OrthoDB, values_from = count)
  
  unloadNamespace("plyr")
  # Count ocurrences of each taxonID per OrthoDB
  query_orthodb_pivot2 <- merge(query_orthodb_stripped,Uniprot_annotated_complexes) %>%
    distinct(OrthoDB_genes, .keep_all=T) %>%
    group_by(OrthoDB, TaxonID) %>%
    summarise(Count = n(), .groups = "drop") %>%
    pivot_wider(names_from = OrthoDB, values_from = Count, values_fill = 0) %>% remove_rownames %>% column_to_rownames(var="TaxonID")

  library(ggnewscale)
  # orthodb DCC-USCO (singular complex components) uniplot - stop inheriting color from other splits
  numhomologs = 1
  
  # Filter and build dataframe - negates need for "filter_if(Orthodb_merge[c(1,7:36)], is.numeric, all_vars((.) >= 2))"
  query_orthodb_filter <- (filter_if(merge(query_orthodb_stripped,Uniprot_annotated_complexes) %>%
                                       distinct(OrthoDB_genes, .keep_all=T) %>%
                                       group_by(OrthoDB, TaxonID) %>%
                                       summarise(Count = n(), .groups = "drop") %>%
                                       pivot_wider(names_from = TaxonID, values_from = Count, values_fill = 0), is.numeric, all_vars((.) >= numhomologs)))
  
  
  query_orthodb_seprow <- separate_rows(merge(
    (filter_if(merge(query_orthodb_stripped,Uniprot_annotated_complexes) %>%
                 distinct(OrthoDB_genes, .keep_all=T) %>%
                 group_by(OrthoDB, TaxonID) %>%
                 summarise(Count = n(), .groups = "drop") %>%
                 pivot_wider(names_from = TaxonID, values_from = Count, values_fill = 0), is.numeric, all_vars((.) >= numhomologs)))
    ,Uniprot_annotated_complexes, all.x =TRUE) %>% distinct(OrthoDB, .keep_all=T) %>% mutate(across(everything(),~ gsub(';$', '', .))) %>% mutate(across(everything(),~ gsub("[[Q].*$", "", .))),ncol(query_orthodb_filter)+5,sep = ";")
  colnames(query_orthodb_seprow)[ncol(query_orthodb_filter)+5] <- "#Complex ac"
  queryorthodb_complexes2_merge <- merge(query_orthodb_seprow,Complexes2, all.x=T)
  
  # Fix heatmap
  ## See applicable labels
  
  # Load plot and get heatmap label positions
  widthno = 9
  firstoffset = 1
  querytreeplot <- ggtree(j)
  querytreeplot <- ggtree(j) + 
    geom_tiplab(size = 3, vjust = 0, aes(label = querytreeplot$data$sci_name))# +
  #xlim(NA, 60)
  orthodb_heatmap <- gheatmap(querytreeplot, query_orthodb_pivot2 %>% dplyr::select(where(~!any(. <= numhomologs-1))), offset = firstoffset, colnames=F, width=widthno) + ylim(-10, NA) + theme(legend.position = "none")
  orthodb_heatmappositions = get_heatmap_column_position(orthodb_heatmap, by="bottom")
  
  # Get protein descriptions from Uniprot list
  orthodb_heatmappositions2 <- merge((merge(query_orthodb_filter,Uniprot_annotated_complexes, all.x =TRUE) %>% distinct(OrthoDB, .keep_all=T))[c(1,ncol((merge(query_orthodb_filter,Uniprot_annotated_complexes, all.x =TRUE) %>% distinct(OrthoDB, .keep_all=T)))-2)], orthodb_heatmappositions, by.x='OrthoDB', by.y='label')
  
  colnames(orthodb_heatmappositions2) <- c("label","Protein.names","x","y")
  
  orthodb_heatmappositionsarranged <-cbind((merge((merge(query_orthodb_filter,Uniprot_annotated_complexes, all.x =TRUE) %>% distinct(OrthoDB, .keep_all=T))[c(1,ncol(query_orthodb_filter)+4)], queryorthodb_complexes2_merge, by.x='OrthoDB', by.y='OrthoDB') %>% group_by(OrthoDB) %>% summarise(Complex = paste(sort(unique(`Protein.names`)), collapse = ";"), across()) %>%
                                               distinct(OrthoDB, .keep_all=T) %>% mutate(Complex=case_when(
                                                 str_detect(OrthoDB, regex("10036721at2759",ignore_case=F)) ~ "MCM complex",
                                                 str_detect(OrthoDB, regex("10060499at2759",ignore_case=F)) ~ "INO80 complex",
                                                 str_detect(OrthoDB, regex("10248520at2759",ignore_case=F)) ~ "T-complex",
                                                 str_detect(OrthoDB, regex("10248617at2759",ignore_case=F)) ~ "Pol (II) complex",
                                                 str_detect(OrthoDB, regex("10248812at2759",ignore_case=F)) ~ "MPP10 complex",
                                                 str_detect(OrthoDB, regex("10249535at2759",ignore_case=F)) ~ "Histone cleavage complex",
                                                 str_detect(OrthoDB, regex("10250002at2759",ignore_case=F)) ~ "Telomerase holoenzyme",
                                                 str_detect(OrthoDB, regex("10250117at2759",ignore_case=F)) ~ "NIAUFX (FeS) assembly",
                                                 str_detect(OrthoDB, regex("10250354at2759",ignore_case=F)) ~ "R2TP co-chaperone",
                                                 str_detect(OrthoDB, regex("10250478at2759",ignore_case=F)) ~ "MetRS",
                                                 str_detect(OrthoDB, regex("10250817at2759",ignore_case=F)) ~ "SRP",
                                                 str_detect(OrthoDB, regex("10250970at2759",ignore_case=F)) ~ "snRNP U2",
                                                 str_detect(OrthoDB, regex("10251154at2759",ignore_case=F)) ~ "SSU processome",
                                                 str_detect(OrthoDB, regex("10251574at2759",ignore_case=F)) ~ "MCM complex",
                                                 str_detect(OrthoDB, regex("10252509at2759",ignore_case=F)) ~ "Proteasome",
                                                 str_detect(OrthoDB, regex("10253113at2759",ignore_case=F)) ~ "CIA (FeS) assembly",
                                                 str_detect(OrthoDB, regex("10253204at2759",ignore_case=F)) ~ "MPP10 complex",
                                                 str_detect(OrthoDB, regex("10253254at2759",ignore_case=F)) ~ "Post-spliceosomal",
                                                 str_detect(OrthoDB, regex("10254527at2759",ignore_case=F)) ~ "ERF1-ERF3 complex",
                                                 str_detect(OrthoDB, regex("10254671at2759",ignore_case=F)) ~ "Kinase Casein CKII",
                                                 str_detect(OrthoDB, regex("10255013at2759",ignore_case=F)) ~ "SNARE",
                                                 str_detect(OrthoDB, regex("10255414at2759",ignore_case=F)) ~ "EIF2A",
                                                 str_detect(OrthoDB, regex("10255768at2759",ignore_case=F)) ~ "Proteasome",
                                                 str_detect(OrthoDB, regex("10256122at2759",ignore_case=F)) ~ "PRP19-associated complex",
                                                 str_detect(OrthoDB, regex("10256289at2759",ignore_case=F)) ~ "COPII",
                                                 str_detect(OrthoDB, regex("10256771at2759",ignore_case=F)) ~ "Proteasome",
                                                 str_detect(OrthoDB, regex("10258825at2759",ignore_case=F)) ~ "eNoSc complex",
                                                 str_detect(OrthoDB, regex("10260625at2759",ignore_case=F)) ~ "Foxo3-Ywhaz complex",
                                                 str_detect(OrthoDB, regex("10260794at2759",ignore_case=F)) ~ "snRNP U2",
                                                 str_detect(OrthoDB, regex("10261556at2759",ignore_case=F)) ~ "BTR complex",
                                                 str_detect(OrthoDB, regex("10262475at2759",ignore_case=F)) ~ "Mitotic checkpoint",
                                                 str_detect(OrthoDB, regex("10262857at2759",ignore_case=F)) ~ "TRAPP II/III",
                                                 str_detect(OrthoDB, regex("10262986at2759",ignore_case=F)) ~ "TFIIH complex",
                                                 str_detect(OrthoDB, regex("10263554at2759",ignore_case=F)) ~ "V-ATPase",
                                                 str_detect(OrthoDB, regex("10264038at2759",ignore_case=F)) ~ "Exosome",
                                                 str_detect(OrthoDB, regex("10264220at2759",ignore_case=F)) ~ "V-ATPase",
                                                 str_detect(OrthoDB, regex("10264728at2759",ignore_case=F)) ~ "NatA complex",
                                                 str_detect(OrthoDB, regex("10264910at2759",ignore_case=F)) ~ "PeBoW complex",
                                                 str_detect(OrthoDB, regex("10265243at2759",ignore_case=F)) ~ "Elongator holoenzyme complex",
                                                 str_detect(OrthoDB, regex("10265628at2759",ignore_case=F)) ~ "serine hydroxymethyltransferase ",
                                                 str_detect(OrthoDB, regex("10265785at2759",ignore_case=F)) ~ "EIF4",
                                                 str_detect(OrthoDB, regex("10266385at2759",ignore_case=F)) ~ "MitoPDH",
                                                 str_detect(OrthoDB, regex("108365at2759",ignore_case=F)) ~ "Kinase Pyruvate",
                                                 str_detect(OrthoDB, regex("1093at2759",ignore_case=F)) ~ "Proteasome",
                                                 str_detect(OrthoDB, regex("1183224at2759",ignore_case=F)) ~ "CCR4-NOT deadenylase",
                                                 str_detect(OrthoDB, regex("1608002at2759",ignore_case=F)) ~ "NuRF complex",
                                                 str_detect(OrthoDB, regex("16120at2759",ignore_case=F)) ~ "Phosphatase TAP42",
                                                 str_detect(OrthoDB, regex("1649088at2759",ignore_case=F)) ~ "Histidine synthase",
                                                 str_detect(OrthoDB, regex("1650at2759",ignore_case=F)) ~ "Exosome",
                                                 str_detect(OrthoDB, regex("1662883at2759",ignore_case=F)) ~ "Tubulin dimer",
                                                 str_detect(OrthoDB, regex("1685042at2759",ignore_case=F)) ~ "EIF2A",
                                                 str_detect(OrthoDB, regex("1698572at2759",ignore_case=F)) ~ "PheRS",
                                                 str_detect(OrthoDB, regex("1706066at2759",ignore_case=F)) ~ "Acetyl CoA-synthase",
                                                 str_detect(OrthoDB, regex("1706657at2759",ignore_case=F)) ~ "MSC",
                                                 str_detect(OrthoDB, regex("1716531at2759",ignore_case=F)) ~ "E3 UL complex (HRD1)",
                                                 str_detect(OrthoDB, regex("1719357at2759",ignore_case=F)) ~ "SNARE",
                                                 str_detect(OrthoDB, regex("1727884at2759",ignore_case=F)) ~ "SRP",
                                                 str_detect(OrthoDB, regex("1732493at2759",ignore_case=F)) ~ "Cyclin C-CDK3 complex",
                                                 str_detect(OrthoDB, regex("1741334at2759",ignore_case=F)) ~ "NUBP1 (FeS) assembly",
                                                 str_detect(OrthoDB, regex("1744952at2759",ignore_case=F)) ~ "MCM complex",
                                                 str_detect(OrthoDB, regex("1748577at2759",ignore_case=F)) ~ "T-complex",
                                                 str_detect(OrthoDB, regex("1856718at2759",ignore_case=F)) ~ "SiR complex",
                                                 str_detect(OrthoDB, regex("1882297at2759",ignore_case=F)) ~ "60S subunit",
                                                 str_detect(OrthoDB, regex("1882346at2759",ignore_case=F)) ~ "MCM complex",
                                                 str_detect(OrthoDB, regex("1918432at2759",ignore_case=F)) ~ "MBD2/NuRD HDAC",
                                                 str_detect(OrthoDB, regex("192611at2759",ignore_case=F)) ~ "GPI-anchor transamidase",
                                                 str_detect(OrthoDB, regex("1926878at2759",ignore_case=F)) ~ "Nuclear ORC",
                                                 str_detect(OrthoDB, regex("1930084at2759",ignore_case=F)) ~ "Phosphatase SIT4/SAP155",
                                                 str_detect(OrthoDB, regex("1933107at2759",ignore_case=F)) ~ "CCR4-NOT deadenylase",
                                                 str_detect(OrthoDB, regex("193499at2759",ignore_case=F)) ~ "Intron-binding complex",
                                                 str_detect(OrthoDB, regex("1937912at2759",ignore_case=F)) ~ "ATAC complex",
                                                 str_detect(OrthoDB, regex("193931at2759",ignore_case=F)) ~ "Kinase AMPK",
                                                 str_detect(OrthoDB, regex("196131at2759",ignore_case=F)) ~ "snRNP U5",
                                                 str_detect(OrthoDB, regex("197206at2759",ignore_case=F)) ~ "MSC",
                                                 str_detect(OrthoDB, regex("2011769at2759",ignore_case=F)) ~ "membrane remodeling complex",
                                                 str_detect(OrthoDB, regex("206088at2759",ignore_case=F)) ~ "XRCC1 complex",
                                                 str_detect(OrthoDB, regex("2110130at2759",ignore_case=F)) ~ "GCN1-GCN20 complex",
                                                 str_detect(OrthoDB, regex("21243at2759",ignore_case=F)) ~ "MSC",
                                                 str_detect(OrthoDB, regex("2127950at2759",ignore_case=F)) ~ "TFIIB/D",
                                                 str_detect(OrthoDB, regex("2143914at2759",ignore_case=F)) ~ "c-Myb-C/EBPbeta complex",
                                                 str_detect(OrthoDB, regex("21449at2759",ignore_case=F)) ~ "Nucleolar repressor",
                                                 str_detect(OrthoDB, regex("2187549at2759",ignore_case=F)) ~ "BUD23/TRM112 methyltransferase",
                                                 str_detect(OrthoDB, regex("2250022at2759",ignore_case=F)) ~ "mTORC I/II",
                                                 str_detect(OrthoDB, regex("2357150at2759",ignore_case=F)) ~ "Nuclear pore complex",
                                                 str_detect(OrthoDB, regex("238316at2759",ignore_case=F)) ~ "PheRS",
                                                 str_detect(OrthoDB, regex("2401965at2759",ignore_case=F)) ~ "COX1 complex",#mitochondrial
                                                 str_detect(OrthoDB, regex("240298at2759",ignore_case=F)) ~ "TIM23 translocase",
                                                 str_detect(OrthoDB, regex("2414538at2759",ignore_case=F)) ~ "Pol (zeta) complex",
                                                 str_detect(OrthoDB, regex("2415936at2759",ignore_case=F)) ~ "E3 UL complex (GID)",
                                                 str_detect(OrthoDB, regex("2423701at2759",ignore_case=F)) ~ "Hsp90-sti1 chaperone",
                                                 str_detect(OrthoDB, regex("246406at2759",ignore_case=F)) ~ "TRAPP II/III",
                                                 str_detect(OrthoDB, regex("248779at2759",ignore_case=F)) ~ "Pol (I/II) complex",
                                                 str_detect(OrthoDB, regex("248923at2759",ignore_case=F)) ~ "PAR polarity complex",
                                                 str_detect(OrthoDB, regex("24966at2759",ignore_case=F)) ~ "EIF3",
                                                 str_detect(OrthoDB, regex("256303at2759",ignore_case=F)) ~ "Nuclear pore complex",
                                                 str_detect(OrthoDB, regex("258143at2759",ignore_case=F)) ~ "tRNA-intron endonuclease",
                                                 str_detect(OrthoDB, regex("26525at2759",ignore_case=F)) ~ "Calprotectin & S100A8",
                                                 str_detect(OrthoDB, regex("268428at2759",ignore_case=F)) ~ "Proteasome",
                                                 str_detect(OrthoDB, regex("268763at2759",ignore_case=F)) ~ "Proteasome",
                                                 str_detect(OrthoDB, regex("269804at2759",ignore_case=F)) ~ "SSU processome",
                                                 str_detect(OrthoDB, regex("270173at2759",ignore_case=F)) ~ "Pol (II) complex",
                                                 str_detect(OrthoDB, regex("270392at2759",ignore_case=F)) ~ "Pol I",
                                                 str_detect(OrthoDB, regex("272141at2759",ignore_case=F)) ~ "B-catenin destruction complex",
                                                 str_detect(OrthoDB, regex("272481at2759",ignore_case=F)) ~ "TFIIH complex",
                                                 str_detect(OrthoDB, regex("273340at2759",ignore_case=F)) ~ "SSU processome",
                                                 str_detect(OrthoDB, regex("27435at2759",ignore_case=F)) ~ "AAA ATPase",
                                                 str_detect(OrthoDB, regex("27923at2759",ignore_case=F)) ~ "SNARE",
                                                 str_detect(OrthoDB, regex("28053at2759",ignore_case=F)) ~ "Adaptor complex",
                                                 str_detect(OrthoDB, regex("282152at2759",ignore_case=F)) ~ "Pol (III) complex",
                                                 str_detect(OrthoDB, regex("284782at2759",ignore_case=F)) ~ "CIA (FeS) assembly",
                                                 str_detect(OrthoDB, regex("28737at2759",ignore_case=F)) ~ "HMC complex",
                                                 str_detect(OrthoDB, regex("2877at2759",ignore_case=F)) ~ "BUD23/TRM112 methyltransferase",
                                                 str_detect(OrthoDB, regex("289038at2759",ignore_case=F)) ~ "SAGA complex",
                                                 str_detect(OrthoDB, regex("294251at2759",ignore_case=F)) ~ "Polarisome",
                                                 str_detect(OrthoDB, regex("29755at2759",ignore_case=F)) ~ "RAD53-ASF1 complex",
                                                 str_detect(OrthoDB, regex("3000483at2759",ignore_case=F)) ~ "RR1 complex",
                                                 str_detect(OrthoDB, regex("308383at2759",ignore_case=F)) ~ "KMT complex",
                                                 str_detect(OrthoDB, regex("3142434at2759",ignore_case=F)) ~ "UTP-B complex",
                                                 str_detect(OrthoDB, regex("3176171at2759",ignore_case=F)) ~ "KIF3 complex",
                                                 str_detect(OrthoDB, regex("336885at2759",ignore_case=F)) ~ "Pol-prim complex",
                                                 str_detect(OrthoDB, regex("340346at2759",ignore_case=F)) ~ "Phosphatase PP2A complex",
                                                 str_detect(OrthoDB, regex("340608at2759",ignore_case=F)) ~ "P4 ATPase",
                                                 str_detect(OrthoDB, regex("342024at2759",ignore_case=F)) ~ "EF1A",
                                                 str_detect(OrthoDB, regex("364224at2759",ignore_case=F)) ~ "GATOR2 complex",
                                                 str_detect(OrthoDB, regex("364892at2759",ignore_case=F)) ~ "snRNP U5",
                                                 str_detect(OrthoDB, regex("372421at2759",ignore_case=F)) ~ "Exosome",
                                                 str_detect(OrthoDB, regex("372487at2759",ignore_case=F)) ~ "Decapping complex",
                                                 str_detect(OrthoDB, regex("3763at2759",ignore_case=F)) ~ "Pol (delta) complex",
                                                 str_detect(OrthoDB, regex("377733at2759",ignore_case=F)) ~ "P4 ATPase",
                                                 str_detect(OrthoDB, regex("405996at2759",ignore_case=F)) ~ "SPOTS complex",
                                                 str_detect(OrthoDB, regex("407658at2759",ignore_case=F)) ~ "a16mannosyltransferase complex",
                                                 str_detect(OrthoDB, regex("411372at2759",ignore_case=F)) ~ "E3 UL complex (lep-2)",
                                                 str_detect(OrthoDB, regex("412748at2759",ignore_case=F)) ~ "CPSF",
                                                 str_detect(OrthoDB, regex("413460at2759",ignore_case=F)) ~ "B-wich complex",
                                                 str_detect(OrthoDB, regex("415696at2759",ignore_case=F)) ~ "GRX4 (FeS) assembly",
                                                 str_detect(OrthoDB, regex("417078at2759",ignore_case=F)) ~ "Kinase PKA",
                                                 str_detect(OrthoDB, regex("419317at2759",ignore_case=F)) ~ "XPC complex",
                                                 str_detect(OrthoDB, regex("4199794at2759",ignore_case=F)) ~ "RFC",
                                                 str_detect(OrthoDB, regex("422728at2759",ignore_case=F)) ~ "AAA ATPase",
                                                 str_detect(OrthoDB, regex("424572at2759",ignore_case=F)) ~ "EIF2B",
                                                 str_detect(OrthoDB, regex("427280at2759",ignore_case=F)) ~ "Tapasin-ERp57 complex",
                                                 str_detect(OrthoDB, regex("427795at2759",ignore_case=F)) ~ "MBD2/NuRD HDAC",
                                                 str_detect(OrthoDB, regex("428577at2759",ignore_case=F)) ~ "60S subunit",
                                                 str_detect(OrthoDB, regex("429533at2759",ignore_case=F)) ~ "Proteasome",
                                                 str_detect(OrthoDB, regex("430051at2759",ignore_case=F)) ~ "RecQ helicase",
                                                 str_detect(OrthoDB, regex("431557at2759",ignore_case=F)) ~ "Proteasome",
                                                 str_detect(OrthoDB, regex("431715at2759",ignore_case=F)) ~ "UTP-A complex",
                                                 str_detect(OrthoDB, regex("441223at2759",ignore_case=F)) ~ "SSU processome",
                                                 str_detect(OrthoDB, regex("446168at2759",ignore_case=F)) ~ "Replication factor C",
                                                 str_detect(OrthoDB, regex("448448at2759",ignore_case=F)) ~ "E3 UL complex (CUL3)",
                                                 str_detect(OrthoDB, regex("49016at2759",ignore_case=F)) ~ "COPII",
                                                 str_detect(OrthoDB, regex("5061070at2759",ignore_case=F)) ~ "FZO1-MGM1-UGO1 complex",
                                                 str_detect(OrthoDB, regex("511529at2759",ignore_case=F)) ~ "TFIIB complex",
                                                 str_detect(OrthoDB, regex("5132116at2759",ignore_case=F)) ~ "INO80 complex",
                                                 str_detect(OrthoDB, regex("534348at2759",ignore_case=F)) ~ "PCNA homotrimer",
                                                 str_detect(OrthoDB, regex("537915at2759",ignore_case=F)) ~ "Kinase PFK-1",
                                                 str_detect(OrthoDB, regex("542917at2759",ignore_case=F)) ~ "COPII",
                                                 str_detect(OrthoDB, regex("550424at2759",ignore_case=F)) ~ "Radial spoke complex",
                                                 str_detect(OrthoDB, regex("5571054at2759",ignore_case=F)) ~ "PeBoW complex",
                                                 str_detect(OrthoDB, regex("5575062at2759",ignore_case=F)) ~ "Condensin I/II",
                                                 str_detect(OrthoDB, regex("5575at2759",ignore_case=F)) ~ "snRNP U5",
                                                 str_detect(OrthoDB, regex("5800476at2759",ignore_case=F)) ~ "B-catenin destruction complex",
                                                 str_detect(OrthoDB, regex("5839at2759",ignore_case=F)) ~ "PSMG3/4 proteasomal chaperone",
                                                 str_detect(OrthoDB, regex("5857104at2759",ignore_case=F)) ~ "MBD2/NuRD HDAC",
                                                 str_detect(OrthoDB, regex("590761at2759",ignore_case=F)) ~ "EIF4",
                                                 str_detect(OrthoDB, regex("6108017at2759",ignore_case=F)) ~ "Myosin complex",
                                                 str_detect(OrthoDB, regex("613763at2759",ignore_case=F)) ~ "Pol (III) complex",
                                                 str_detect(OrthoDB, regex("63267at2759",ignore_case=F)) ~ "Kinase Atk-1/2",
                                                 str_detect(OrthoDB, regex("6353017at2759",ignore_case=F)) ~ "Ndc80 complex",
                                                 str_detect(OrthoDB, regex("64353at2759",ignore_case=F)) ~ "CPSF complex",
                                                 str_detect(OrthoDB, regex("64767at2759",ignore_case=F)) ~ "SKI complex",
                                                 str_detect(OrthoDB, regex("64767at2759",ignore_case=F)) ~ "SKI complex",
                                                 str_detect(OrthoDB, regex("6500128at2759",ignore_case=F)) ~ "ABC transporter",
                                                 str_detect(OrthoDB, regex("6513042at2759",ignore_case=F)) ~ "NMD complex",
                                                 str_detect(OrthoDB, regex("6600758at2759",ignore_case=F)) ~ "E3 UL complex (RANBP2)",
                                                 str_detect(OrthoDB, regex("660555at2759",ignore_case=F)) ~ "DISP complex",
                                                 str_detect(OrthoDB, regex("6780543at2759",ignore_case=F)) ~ "snRNP C/D",
                                                 str_detect(OrthoDB, regex("68056at2759",ignore_case=F)) ~ "MARS",
                                                 str_detect(OrthoDB, regex("7851174at2759",ignore_case=F)) ~ "E2 UBC (UBC13)",
                                                 str_detect(OrthoDB, regex("7875889at2759",ignore_case=F)) ~ "40S subunit",
                                                 str_detect(OrthoDB, regex("8068875at2759",ignore_case=F)) ~ "TRAMP complex",
                                                 str_detect(OrthoDB, regex("844at2759",ignore_case=F)) ~ "MCM complex",
                                                 str_detect(OrthoDB, regex("8962942at2759",ignore_case=F)) ~ "E3 UL complex (CUL4B)",
                                                 str_detect(OrthoDB, regex("9332038at2759",ignore_case=F)) ~ "Kinase (Egg-3/4/5 MBK-2",
                                                 str_detect(OrthoDB, regex("9984419at2759",ignore_case=F)) ~ "E3 UL complex (RAD6-18)",
                                                 str_detect(OrthoDB, regex("9991317at2759",ignore_case=F)) ~ "TFIIIC complex",)) %>%
                                               mutate(result=case_when(
                                                 str_detect(OrthoDB, regex("10036721at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("10060499at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("10248520at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("10248617at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("10248812at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("10249535at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("10250002at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("10250117at2759",ignore_case=F)) ~ "xOther",
                                                 str_detect(OrthoDB, regex("10250354at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("10250478at2759",ignore_case=F)) ~ "Translation",
                                                 str_detect(OrthoDB, regex("10250817at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("10250970at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("10251154at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("10251574at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("10252509at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("10253113at2759",ignore_case=F)) ~ "xOther",
                                                 str_detect(OrthoDB, regex("10253204at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("10253254at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("10254527at2759",ignore_case=F)) ~ "Translation",
                                                 str_detect(OrthoDB, regex("10254671at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("10255013at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("10255414at2759",ignore_case=F)) ~ "Translation",
                                                 str_detect(OrthoDB, regex("10255768at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("10256122at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("10256289at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("10256771at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("10258825at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("10260625at2759",ignore_case=F)) ~ "xOther",
                                                 str_detect(OrthoDB, regex("10260794at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("10261556at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("10262475at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("10262857at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("10262986at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("10263554at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("10264038at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("10264220at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("10264728at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("10264910at2759",ignore_case=F)) ~ "Translation",
                                                 str_detect(OrthoDB, regex("10265243at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("10265628at2759",ignore_case=F)) ~ "xOther",
                                                 str_detect(OrthoDB, regex("10265785at2759",ignore_case=F)) ~ "Translation",
                                                 str_detect(OrthoDB, regex("10266385at2759",ignore_case=F)) ~ "xOther",
                                                 str_detect(OrthoDB, regex("108365at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("1093at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("1183224at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("1608002at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("16120at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("1649088at2759",ignore_case=F)) ~ "xOther",
                                                 str_detect(OrthoDB, regex("1650at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("1662883at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("1685042at2759",ignore_case=F)) ~ "Translation",
                                                 str_detect(OrthoDB, regex("1698572at2759",ignore_case=F)) ~ "Translation",
                                                 str_detect(OrthoDB, regex("1706066at2759",ignore_case=F)) ~ "xOther",
                                                 str_detect(OrthoDB, regex("1706657at2759",ignore_case=F)) ~ "Translation",
                                                 str_detect(OrthoDB, regex("1716531at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("1719357at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("1727884at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("1732493at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("1741334at2759",ignore_case=F)) ~ "xOther",
                                                 str_detect(OrthoDB, regex("1744952at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("1748577at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("1856718at2759",ignore_case=F)) ~ "xOther",
                                                 str_detect(OrthoDB, regex("1882297at2759",ignore_case=F)) ~ "Translation",
                                                 str_detect(OrthoDB, regex("1882346at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("1918432at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("192611at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("1926878at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("1930084at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("1933107at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("193499at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("1937912at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("193931at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("196131at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("197206at2759",ignore_case=F)) ~ "Translation",
                                                 str_detect(OrthoDB, regex("2011769at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("206088at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("2110130at2759",ignore_case=F)) ~ "Translation",
                                                 str_detect(OrthoDB, regex("21243at2759",ignore_case=F)) ~ "Translation",
                                                 str_detect(OrthoDB, regex("2127950at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("2143914at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("21449at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("2187549at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("2250022at2759",ignore_case=F)) ~ "Translation",
                                                 str_detect(OrthoDB, regex("2357150at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("238316at2759",ignore_case=F)) ~ "Translation",
                                                 str_detect(OrthoDB, regex("2401965at2759",ignore_case=F)) ~ "Translation",#mitochondrial
                                                 str_detect(OrthoDB, regex("240298at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("2414538at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("2415936at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("2423701at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("246406at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("248779at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("248923at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("24966at2759",ignore_case=F)) ~ "Translation",
                                                 str_detect(OrthoDB, regex("256303at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("258143at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("26525at2759",ignore_case=F)) ~ "xOther",
                                                 str_detect(OrthoDB, regex("268428at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("268763at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("269804at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("270173at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("270392at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("272141at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("272481at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("273340at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("27435at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("27923at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("28053at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("282152at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("284782at2759",ignore_case=F)) ~ "xOther",
                                                 str_detect(OrthoDB, regex("28737at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("2877at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("289038at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("294251at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("29755at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("3000483at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("308383at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("3142434at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("3176171at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("336885at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("340346at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("340608at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("342024at2759",ignore_case=F)) ~ "Translation",
                                                 str_detect(OrthoDB, regex("364224at2759",ignore_case=F)) ~ "xOther",
                                                 str_detect(OrthoDB, regex("364892at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("372421at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("372487at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("3763at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("377733at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("405996at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("407658at2759",ignore_case=F)) ~ "xOther",
                                                 str_detect(OrthoDB, regex("411372at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("412748at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("413460at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("415696at2759",ignore_case=F)) ~ "xOther",
                                                 str_detect(OrthoDB, regex("417078at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("419317at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("4199794at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("422728at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("424572at2759",ignore_case=F)) ~ "Translation",
                                                 str_detect(OrthoDB, regex("427280at2759",ignore_case=F)) ~ "xOther",
                                                 str_detect(OrthoDB, regex("427795at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("428577at2759",ignore_case=F)) ~ "Translation",
                                                 str_detect(OrthoDB, regex("429533at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("430051at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("431557at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("431715at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("441223at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("446168at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("448448at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("49016at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("5061070at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("511529at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("5132116at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("534348at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("537915at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("542917at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("550424at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("5571054at2759",ignore_case=F)) ~ "Translation",
                                                 str_detect(OrthoDB, regex("5575062at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("5575at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("5800476at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("5839at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("5857104at2759",ignore_case=F)) ~ "DNA",
                                                 str_detect(OrthoDB, regex("590761at2759",ignore_case=F)) ~ "Translation",
                                                 str_detect(OrthoDB, regex("6108017at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("613763at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("63267at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("6353017at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("64353at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("64767at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("64767at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("6500128at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("6513042at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("6600758at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("660555at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("6780543at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("68056at2759",ignore_case=F)) ~ "Translation",
                                                 str_detect(OrthoDB, regex("7851174at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("7875889at2759",ignore_case=F)) ~ "Translation",
                                                 str_detect(OrthoDB, regex("8068875at2759",ignore_case=F)) ~ "RNA",
                                                 str_detect(OrthoDB, regex("844at2759",ignore_case=F)) ~ "Transport/Structure",
                                                 str_detect(OrthoDB, regex("8962942at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("9332038at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("9984419at2759",ignore_case=F)) ~ "PTM",
                                                 str_detect(OrthoDB, regex("9991317at2759",ignore_case=F)) ~ "RNA",)) %>%
                                               arrange(result,Complex))[c(1,2,ncol(query_orthodb_filter)+27)],orthodb_heatmappositions[2:3])
  colnames(orthodb_heatmappositionsarranged[1])<-"label"
  
  # split resulting dataframe
  orthodb_splits <- split(orthodb_heatmappositionsarranged, orthodb_heatmappositionsarranged$result)
  
  labelspace = orthodb_heatmappositions$x[2]-orthodb_heatmappositions$x[1]
  
  #firstoffset = 8
  secondoffset = firstoffset + labelspace*nrow(orthodb_splits$DNA)
  thirdoffset = secondoffset + labelspace*nrow(orthodb_splits$PTM)
  fourthoffset = thirdoffset + labelspace*nrow(orthodb_splits$RNA)
  fifthoffset = fourthoffset + labelspace*nrow(orthodb_splits$Translation)
  sixthoffset = fifthoffset + labelspace*nrow(orthodb_splits$`Transport/Structure`)
  
  queryodbplot1 <- gheatmap(gheatmap(querytreeplot, query_orthodb_pivot2 %>% dplyr::select(where(~!any(. <= numhomologs-1))) %>% dplyr::select(all_of(pull(orthodb_splits$DNA, 1))), offset = firstoffset, colnames=F, width=widthno/nrow(orthodb_heatmappositionsarranged)*nrow(orthodb_splits$DNA)) + ylim(-10, NA) + geom_text(data=orthodb_splits$DNA, aes(x, y, label=Complex, colour = result), nudge_y=-0.5, angle=45, hjust = 1, size=2) + scale_fill_gradient(name="genes/OG", low="white", high="green4", limits = c(-20,300), trans="log") +theme(legend.position = "below") + geom_text(aes(label="DNA", y=45, x=orthodb_splits$DNA$x[12]), size=4), 
                            query_orthodb_pivot2 %>% dplyr::select(where(~!any(. <= numhomologs-1))) %>% dplyr::select(all_of(pull(orthodb_splits$DNA, 1))), offset = firstoffset, colnames_position = "top", colnames_offset_y=8, colnames_angle = "90", font.size = 2,width=widthno/nrow(orthodb_heatmappositionsarranged)*nrow(orthodb_splits$DNA)) + ylim(-10, NA) + geom_text(data=orthodb_splits$DNA, aes(x, y, label=Complex, colour = result), nudge_y=-0.5, angle=45, hjust = 1, size=2) + scale_fill_gradient(name="genes/OG", low="white", high="green4", limits = c(-20,300), trans="log") + theme(legend.position = "below") + geom_text(aes(label="DNA", y=45, x=orthodb_splits$DNA$x[12]), size=4)
  queryodbplot2 <- gheatmap(gheatmap(queryodbplot1 + new_scale_fill(), query_orthodb_pivot2 %>% dplyr::select(where(~!any(. <= numhomologs-1))) %>% dplyr::select(all_of(pull(orthodb_splits$PTM, 1))), offset = secondoffset, colnames=F, width=widthno/nrow(orthodb_heatmappositionsarranged)*nrow(orthodb_splits$PTM)) +  ylim(-10, NA) + geom_text(data=orthodb_splits$PTM, aes(x, y, label=Complex, colour = result), nudge_y=-0.5, angle=45, hjust = 1, size=2) + scale_fill_gradient(name="genes/OG", low="white", high="blue", limits = c(-20,300), trans="log") + guides(fill=FALSE) + theme(legend.position = "right", legend.key.height = unit(0.3, 'cm'), legend.box.margin=margin(0,0,0,0)) + geom_text(aes(label="PTM", y=45, x=orthodb_splits$PTM$x[10]), size=4) + new_scale_fill(), query_orthodb_pivot2 %>% dplyr::select(where(~!any(. <= numhomologs-1))) %>% dplyr::select(all_of(pull(orthodb_splits$PTM, 1))), offset = secondoffset, colnames_position = "top", colnames_offset_y=8, colnames_angle = "90", font.size = 2, width=widthno/nrow(orthodb_heatmappositionsarranged)*nrow(orthodb_splits$PTM)) +  ylim(-10, NA) + geom_text(data=orthodb_splits$PTM, aes(x, y, label=Complex, colour = result), nudge_y=-0.5, angle=45, hjust = 1, size=2) + scale_fill_gradient(name="genes/OG", low="white", high="blue", limits = c(-20,300), trans="log") + guides(fill=FALSE) + theme(legend.position = "right", legend.key.height = unit(0.3, 'cm'), legend.box.margin=margin(0,0,0,0)) + geom_text(aes(label="PTM", y=45, x=orthodb_splits$PTM$x[10]), size=4)
  queryodbplot3 <- gheatmap(gheatmap(queryodbplot2 + new_scale_fill(), query_orthodb_pivot2 %>% dplyr::select(where(~!any(. <= numhomologs-1))) %>% dplyr::select(all_of(pull(orthodb_splits$RNA, 1))), offset = thirdoffset, colnames=F, width=widthno/nrow(orthodb_heatmappositionsarranged)*nrow(orthodb_splits$RNA)) +  ylim(-10, NA) + geom_text(data=orthodb_splits$RNA, aes(x, y, label=Complex, colour = result), nudge_y=-0.5, angle=45, hjust = 1, size=2) + scale_fill_gradient(name="genes/OG", low="white", high="maroon", limits = c(-20,300), trans="log") + guides(fill=FALSE) + theme(legend.position = "right", legend.key.height = unit(0.3, 'cm'), legend.box.margin=margin(0,0,0,0)) + geom_text(aes(label="RNA", y=45, x=orthodb_splits$RNA$x[15]), size=4) + new_scale_fill(), query_orthodb_pivot2 %>% dplyr::select(where(~!any(. <= numhomologs-1))) %>% dplyr::select(all_of(pull(orthodb_splits$RNA, 1))), offset = thirdoffset, colnames_position = "top", colnames_offset_y=8, colnames_angle = "90", font.size = 2, width=widthno/nrow(orthodb_heatmappositionsarranged)*nrow(orthodb_splits$RNA)) +  ylim(-10, NA) + geom_text(data=orthodb_splits$RNA, aes(x, y, label=Complex, colour = result), nudge_y=-0.5, angle=45, hjust = 1, size=2) + scale_fill_gradient(name="genes/OG", low="white", high="maroon", limits = c(-20,300), trans="log") + guides(fill=FALSE) + theme(legend.position = "right", legend.key.height = unit(0.3, 'cm'), legend.box.margin=margin(0,0,0,0)) + geom_text(aes(label="RNA", y=45, x=orthodb_splits$RNA$x[15]), size=4)
  queryodbplot4 <- gheatmap(gheatmap(queryodbplot3 + new_scale_fill(), query_orthodb_pivot2 %>% dplyr::select(where(~!any(. <= numhomologs-1))) %>% dplyr::select(all_of(pull(orthodb_splits$Translation, 1))), offset = fourthoffset, colnames=F, width=widthno/nrow(orthodb_heatmappositionsarranged)*nrow(orthodb_splits$Translation)) +  ylim(-10, NA) + geom_text(data=orthodb_splits$Translation, aes(x, y, label=Complex, colour = result), nudge_y=-0.5, angle=45, hjust = 1, size=2) + scale_fill_gradient(name="genes/OG", low="white", high="khaki4", limits = c(-20,300), trans="log") + guides(fill=FALSE) + theme(legend.position = "right", legend.key.height = unit(0.3, 'cm'), legend.box.margin=margin(0,0,0,0)) + geom_text(aes(label="Translation", y=45, x=orthodb_splits$Translation$x[9]), size=4) + new_scale_fill(), query_orthodb_pivot2 %>% dplyr::select(where(~!any(. <= numhomologs-1))) %>% dplyr::select(all_of(pull(orthodb_splits$Translation, 1))), offset = fourthoffset, colnames_position = "top", colnames_offset_y=8, colnames_angle = "90", font.size = 2, width=widthno/nrow(orthodb_heatmappositionsarranged)*nrow(orthodb_splits$Translation)) +  ylim(-10, NA) + geom_text(data=orthodb_splits$Translation, aes(x, y, label=Complex, colour = result), nudge_y=-0.5, angle=45, hjust = 1, size=2) + scale_fill_gradient(name="genes/OG", low="white", high="khaki4", limits = c(-20,300), trans="log") + guides(fill=FALSE) + theme(legend.position = "right", legend.key.height = unit(0.3, 'cm'), legend.box.margin=margin(0,0,0,0)) + geom_text(aes(label="Translation", y=45, x=orthodb_splits$Translation$x[9]), size=4)
  queryodbplot5 <- gheatmap(gheatmap(queryodbplot4 + new_scale_fill(), query_orthodb_pivot2 %>% dplyr::select(where(~!any(. <= numhomologs-1))) %>% dplyr::select(all_of(pull(orthodb_splits$`Transport/Structure`, 1))), offset = fifthoffset, colnames=F, width=widthno/nrow(orthodb_heatmappositionsarranged)*nrow(orthodb_splits$`Transport/Structure`)) +  ylim(-10, NA) + geom_text(data=orthodb_splits$`Transport/Structure`, aes(x, y, label=Complex, colour = result), nudge_y=-0.5, angle=45, hjust = 1, size=2) + scale_fill_gradient(name="genes/OG", low="white", high="magenta", limits = c(-20,300), trans="log") + guides(fill=FALSE) + theme(legend.position = "right", legend.key.height = unit(0.3, 'cm'), legend.box.margin=margin(0,0,0,0)) + geom_text(aes(label="Transport/Structure", y=45, x=orthodb_splits$`Transport/Structure`$x[20]), size=4) + new_scale_fill(), query_orthodb_pivot2 %>% dplyr::select(where(~!any(. <= numhomologs-1))) %>% dplyr::select(all_of(pull(orthodb_splits$`Transport/Structure`, 1))), offset = fifthoffset, colnames_position = "top", colnames_offset_y=8, colnames_angle = "90", font.size = 2, width=widthno/nrow(orthodb_heatmappositionsarranged)*nrow(orthodb_splits$`Transport/Structure`)) +  ylim(-10, NA) + geom_text(data=orthodb_splits$`Transport/Structure`, aes(x, y, label=Complex, colour = result), nudge_y=-0.5, angle=45, hjust = 1, size=2) + scale_fill_gradient(name="genes/OG", low="white", high="magenta", limits = c(-20,300), trans="log") + guides(fill=FALSE) + theme(legend.position = "right", legend.key.height = unit(0.3, 'cm'), legend.box.margin=margin(0,0,0,0)) + geom_text(aes(label="Transport/Structure", y=45, x=orthodb_splits$`Transport/Structure`$x[20]), size=4)
  queryodbplot6 <- gheatmap(gheatmap(queryodbplot5 + new_scale_fill(), query_orthodb_pivot2 %>% dplyr::select(where(~!any(. <= numhomologs-1))) %>% dplyr::select(all_of(pull(orthodb_splits$xOther, 1))), offset = sixthoffset, colnames=F, width=widthno/nrow(orthodb_heatmappositionsarranged)*nrow(orthodb_splits$xOther)) +  ylim(-10, NA) + geom_text(data=orthodb_splits$xOther, aes(x, y, label=Complex, colour = result), nudge_y=-0.5, angle=45, hjust = 1, size=2) + scale_fill_gradient(name="genes/OG", low="white", high="orange", limits = c(-20,300), trans="log") + guides(fill=FALSE) + theme(legend.position = "right", legend.key.height = unit(0.3, 'cm'), legend.box.margin=margin(0,0,0,0)) + geom_text(aes(label="xOther", y=45, x=orthodb_splits$xOther$x[7]), size=4) + new_scale_fill(), query_orthodb_pivot2 %>% dplyr::select(where(~!any(. <= numhomologs-1))) %>% dplyr::select(all_of(pull(orthodb_splits$xOther, 1))), offset = sixthoffset, colnames_position = "top", colnames_offset_y=8, colnames_angle = "90", font.size = 2, width=widthno/nrow(orthodb_heatmappositionsarranged)*nrow(orthodb_splits$xOther)) +  ylim(-10, NA) + geom_text(data=orthodb_splits$xOther, aes(x, y, label=Complex, colour = result), nudge_y=-0.5, angle=45, hjust = 1, size=2) + scale_fill_gradient(name="genes/OG", low="white", high="orange", limits = c(-20,300), trans="log") + guides(fill=FALSE) + theme(legend.position = "right", legend.key.height = unit(0.3, 'cm'), legend.box.margin=margin(0,0,0,0)) + geom_text(aes(label="xOther", y=45, x=orthodb_splits$xOther$x[7]), size=4)
  
  query_repertoire_discordance <- (merge((merge(distinct(merge(
    filter_if(query_orthodb_stripped %>%
                distinct(OrthoDB_genes, .keep_all=T) %>%
                group_by(OrthoDB, TaxonID) %>%
                summarise(Count = n(), .groups = "drop") %>%
                pivot_wider(names_from = TaxonID, values_from = Count, values_fill = 0), is.numeric, all_vars((.) >= numhomologs)),
    Uniprot_annotated_complexes %>% mutate(CPXs = strsplit(as.character(ComplexPortal), ";")) %>% unnest(CPXs))[c(1:(ncol(query_orthodb_filter)),(ncol(query_orthodb_filter)+2):(ncol(query_orthodb_filter)+6))]),orthodb_heatmappositionsarranged) %>% distinct(OrthoDB, .keep_all=T) %>% mutate(Mean = rowMeans(across(2:(ncol(query_orthodb_filter))), na.rm = TRUE)) %>% rowwise() %>% mutate(sd = sd(c_across(2:(ncol(query_orthodb_filter))))) %>% pivot_longer(cols=seq(2,(ncol(query_orthodb_filter))), names_to = "Taxa")) 
    %>% mutate(IsElevated=pick(14) > Mean + 2 * sd & !Taxa %in% c("Allomyces macrogynus ATCC 38327","Paramecium tetraurelia","Ricinus communis","Physcomitrium patens","Arabidopsis thaliana","Emiliania huxleyi CCMP1516"))
    %>% mutate(IsRestricted=pick(14) < Mean - 1.2 * sd & !Taxa %in% c("Encephalitozoon cuniculi","Theileria parva","Babesia bovis","Toxoplasma gondii","Plasmodiophora brassicae"))
    %>% filter(IsRestricted == "TRUE" | IsElevated == "TRUE" ),queryodbplot6$data[c('label','y')],by.x ='Taxa', by.y='label', all.x=T) %>% mutate(xmin=x-(labelspace/2)) %>% mutate(xmax=x+(labelspace/2)) %>% mutate(ymin=y.y-0.5) %>% mutate(ymax=y.y+0.5) %>% mutate(alpha=0) %>% mutate(color=case_when(
      IsRestricted == "TRUE" ~ "blue",
      IsElevated == "TRUE" ~ "red")))
  
  query_repertoire_annotations <- lapply(seq_len(nrow(query_repertoire_discordance)), function(i) {
    annotate("rect",
             xmin  = query_repertoire_discordance$xmin[i],
             xmax  = query_repertoire_discordance$xmax[i],
             ymin  = query_repertoire_discordance$ymin[i],
             ymax  = query_repertoire_discordance$ymax[i],
             alpha = query_repertoire_discordance$alpha[i],
             color = query_repertoire_discordance$color[i])
  })
  
ggsave(file= paste0("Results/",filename, "_plot.svg"), plot=queryodbplot6+query_repertoire_annotations, width=18, height=12)