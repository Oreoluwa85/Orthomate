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

# i <- read.table("/Users/of2/Documents/Complex_Portal/ComPred/microsporidia.tab", sep="\t", header=F)
# j <- read.nhx("/Users/of2/Documents/Complex_Portal/ComPred/microsporidia.nwk")
# setwd("/Users/of2/Documents/Complex_Portal/ComPred/")

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

# 510 complexes are conserved across eukaryotes according to OrthoDB
# Orthodb_pivot2 %>% select(where(~!any(. == 0)))

numhomologs = 1

# Filter and build dataframe - negates need for "filter_if(Orthodb_merge[c(1,7:36)], is.numeric, all_vars((.) >= 2))"
query_orthodb_filter <- (filter_if(merge(query_orthodb_stripped,Uniprot_annotated_complexes) %>%
                                     distinct(OrthoDB_genes, .keep_all=T) %>%
                                     group_by(OrthoDB, TaxonID) %>%
                                     summarise(Count = n(), .groups = "drop") %>%
                                     pivot_wider(names_from = TaxonID, values_from = Count, values_fill = 0), is.numeric, all_vars((.) >= numhomologs)))

querytaxtotals <- as.data.frame(colSums(Filter(is.numeric, query_orthodb_pivot)))
querytaxtotals <- cbind(rownames(querytaxtotals), data.frame(querytaxtotals, row.names=NULL))

querytreeplot <- ggtree(j)
querytreeplot <- ggtree(j) + 
  geom_tiplab(size = 3, vjust = 0, aes(label = querytreeplot$data$sci_name))

colnames(querytaxtotals) <- c("label3","OG_genes")
querytaxcompOGtotals <- as.data.frame(colSums(filter_if(query_orthodb_pivot[2:ncol(query_orthodb_pivot)], is.numeric, all_vars((.) >= 1))))
querytaxcompOGtotals <- cbind(rownames(querytaxcompOGtotals), data.frame(querytaxcompOGtotals, row.names=NULL))
colnames(querytaxcompOGtotals) <- c("label3", "compOG_genes")
querytaxtotalsum <-merge(querytaxtotals, querytreeplot$data[!is.na(querytreeplot$data$label), ], by.x='label3', by.y='label')
querytaxtotalsum <- merge(querytaxtotalsum,querytaxcompOGtotals)
querytaxtotalsum <- querytaxtotalsum[c('label3','OG_genes','compOG_genes')]

queryOGperspecies <- querytreeplot + 
  ggtreeExtra::geom_fruit(data=as.data.frame(querytaxtotalsum),
                          geom=geom_bar,
                          mapping=aes(y=label3, x=as.numeric(compOG_genes)),#, color=tree31setcolors
                          axis.params=list(axis="xy", text.size = 4, title="Total genes mapped \nto compOGs"),
                          offset=12,
                          stat="identity",
                          orientation="y",
                          pwidth=4,
                          show.legend = F) +
  ggtreeExtra::geom_fruit(data=as.data.frame(querytaxtotalsum),
                          geom=geom_bar,
                          mapping=aes(y=label3, x=as.numeric(OG_genes)),#, color=tree31setcolors
                          axis.params=list(axis="xy", text.size = 4, title="Total genes mapped \nto OGs"),
                          offset=-8,
                          stat="identity",
                          orientation="y",
                          pwidth=4,
                          show.legend = F)

ggsave(file= paste0("Results/", filename, "totals_plot.svg"), plot=queryOGperspecies, width=8, height=8)