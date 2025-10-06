library(ggplot2)
library(scales)
library(tidyverse)
library(dplyr)
library(ggtext)
library(ggnewscale)
library(ggExtra)
library(gtExtras)
library(egg)
library(servr)
library(DiagrammeR)
library(webshot2)
library(data.table)
args <- commandArgs(TRUE)
i <- read.table(args[1], sep="\t", header=F)
j <- args[2]
filename <- gsub("TabFiles/","", args[1])
# setwd("/Users/of2/Documents/Complex_Portal/ComPred/")
# i <- read.table("/Users/of2/Documents/Complex_Portal/ComPred/microsporidia.tab", sep="\t", header=F)
# filename <- "microsporidia"

Orthodb_results <- read.table("/Users/of2/Downloads/20250317_31species_OGsnew.tab", sep="\t", header=F)
colnames(Orthodb_results) <- c("OrthoDB","OrthoDB_genes")
Orthodb_stripped <- Orthodb_results %>%
  mutate(TaxonID = sub("_.*", "", OrthoDB_genes))
Orthodb_stripped <- Orthodb_stripped[Orthodb_stripped$OrthoDB %like% "at2759", ]

Uniprot_annotated_complexes<- as.data.frame(read_csv("Resources/Uniprot_annotated_complexes.csv"))
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

query_cpxs_df <- Uniprot_annotated_complexes %>%
  mutate(CPXs = strsplit(as.character(ComplexPortal), ";")) %>%
  unnest(CPXs) %>%
  dplyr::select(OrthoDB, CPXs) %>%
  distinct() %>%
  mutate(across(everything(),~ gsub(';$', '', .)))

# Step 2: Generate the OrthoDB x TaxonID count matrix
unloadNamespace('plyr')
query_orthodb_counts <- filter_if(query_orthodb_stripped %>%
                                    distinct(OrthoDB_genes, .keep_all = TRUE) %>%
                                    group_by(OrthoDB, TaxonID) %>%
                                    summarise(Count = n(), .groups = "drop") %>%
                                    pivot_wider(names_from = TaxonID, values_from = Count, values_fill = 0), is.numeric, all_vars((.) >= 1))[c("OrthoDB",j)]

# Step 3: Join count data to complex info
query_cpxs_counts <- query_cpxs_df %>%
  inner_join(query_orthodb_counts, by = "OrthoDB")

# Step 4: For each CPX, get all OrthoDB pairs
query_pairwise_scores <- query_cpxs_counts %>%
  group_by(CPXs) %>%
  filter(n() > 1) %>%
  summarise(pairs = list(combn(OrthoDB, 2, simplify = FALSE)), .groups = "drop") %>%
  unnest(pairs) %>%
  mutate(OrthoDB1 = map_chr(pairs, 1),
         OrthoDB2 = map_chr(pairs, 2)) %>%
  dplyr::select(CPXs, OrthoDB1, OrthoDB2)

# Identify taxon columns
query_taxon_cols <- setdiff(colnames(query_orthodb_counts), "OrthoDB")

query_taxon_cols <- setdiff(colnames(query_orthodb_counts), "OrthoDB")[1]

# Step 5: Join back counts for both OrthoDBs and rename only taxon columns
query_joined_counts <- query_pairwise_scores %>%
  left_join(query_orthodb_counts, by = c("OrthoDB1" = "OrthoDB")) %>%
  rename_with(~ paste0(., "_1"), all_of(query_taxon_cols)) %>%
  left_join(query_orthodb_counts, by = c("OrthoDB2" = "OrthoDB")) %>%
  rename_with(~ paste0(., "_2"), all_of(query_taxon_cols))

# Step 6: Calculate pairwise difference of squares for each taxon
query_score_df <- query_joined_counts %>%
  rowwise() %>%
  #mutate(score = sum((c_across(ends_with("_1")) + c_across(ends_with("_2")))^2)) %>%
  mutate(score = sum((c_across(ends_with("_1")) + c_across(ends_with("_2"))))) %>%
  ungroup()

# Step 7: Heatmap-ready dataframe
query_heatmap_df <- query_score_df %>%
  dplyr::select(OrthoDB1, OrthoDB2, score)

# Optional: Create symmetric entries (OrthoDB2 vs OrthoDB1)
query_heatmap_df_sym <- bind_rows(
  query_heatmap_df,
  query_heatmap_df %>% rename(OrthoDB1 = OrthoDB2, OrthoDB2 = OrthoDB1)
)


## Load reference grid
cpxs_df <- Uniprot_annotated_complexes %>%
  mutate(CPXs = strsplit(as.character(ComplexPortal), ";")) %>%
  unnest(CPXs) %>%
  dplyr::select(OrthoDB, CPXs) %>%
  distinct() %>%
  mutate(across(everything(),~ gsub(';$', '', .)))
orthodb_counts <- filter_if(Orthodb_stripped %>%
                      distinct(OrthoDB_genes, .keep_all = TRUE) %>%
                      group_by(OrthoDB, TaxonID) %>%
                      summarise(Count = n(), .groups = "drop") %>%
                      pivot_wider(names_from = TaxonID, values_from = Count, values_fill = 0), is.numeric, all_vars((.) >= 1))
cpxs_counts <- cpxs_df %>%
  inner_join(orthodb_counts, by = "OrthoDB")
pairwise_scores <- cpxs_counts %>%
  group_by(CPXs) %>%
  filter(n() > 1) %>%
  summarise(pairs = list(combn(OrthoDB, 2, simplify = FALSE)), .groups = "drop") %>%
  unnest(pairs) %>%
  mutate(OrthoDB1 = map_chr(pairs, 1),
         OrthoDB2 = map_chr(pairs, 2)) %>%
  dplyr::select(CPXs, OrthoDB1, OrthoDB2)
taxon_cols <- setdiff(colnames(orthodb_counts), "OrthoDB")
joined_counts <- pairwise_scores %>%
  left_join(orthodb_counts, by = c("OrthoDB1" = "OrthoDB")) %>%
  rename_with(~ paste0(., "_1"), all_of(taxon_cols)) %>%
  left_join(orthodb_counts, by = c("OrthoDB2" = "OrthoDB")) %>%
  rename_with(~ paste0(., "_2"), all_of(taxon_cols))
score_df <- joined_counts %>%
  rowwise() %>%
  mutate(score = sum((c_across(ends_with("_1")) + c_across(ends_with("_2"))))) %>%
  ungroup()
heatmap_df <- score_df %>%
  dplyr::select(OrthoDB1, OrthoDB2, score)
heatmap_df_sym <- bind_rows(
  heatmap_df,
  heatmap_df %>% rename(OrthoDB1 = OrthoDB2, OrthoDB2 = OrthoDB1)
)

# Get all unique x-axis OrthoDB1 labels (in the correct order)
pairwise_x_labels <- unique(arrange(heatmap_df_sym, score)$OrthoDB1)

# Create a dataframe with labels and assign groups
# Example: manually set some groups, others as "Other"
label_group_df <- tibble(
  OrthoDB = pairwise_x_labels,
  Group = case_when(
    OrthoDB %in% c("10250970at2759","10260794at2759","10256122at2759","5575at2759","364892at2759","196131at2759") ~ "spliceosome",
    OrthoDB %in% c("431557at2759", "10255768at2759", "268428at2759","268763at2759", "429533at2759", "10252509at2759", "10256771at2759","1093at2759") ~ "proteasome",
    OrthoDB %in% c("10251154at2759","441223at2759","273340at2759","269804at2759","6780543at2759") ~ "sr_processome",
    OrthoDB %in% c("1744952at2759","844at2759","10251574at2759","1882346at2759","10036721at2759") ~ "MCM_complex_2_6",
    OrthoDB %in% c("64224at2759","256303at2759","49016at2759") ~ "Nuclear_pore_complex",
    OrthoDB %in% c("613763at2759","248779at2759","270392at2759","10248617at2759","270173at2759","282152at2759") ~ "RNA_polymerase_III_complex",
    OrthoDB %in% c("10256289at2759","542917at2759") ~ "COPII_vesicle_coat_complex",
    OrthoDB %in% c("1882346at2759") ~ "MCM3",
    OrthoDB %in% c("10249535at2759","64353at2759") ~ "Histone_pre-mRNA_cleavage_complex",
    OrthoDB %in% c("10264038at2759","1650at2759","372421at2759","2250022at2759") ~ "Exosome",
    OrthoDB %in% c("1698572at2759","238316at2759") ~ "Phenylalanyl-tRNA synthetase complex",
    # Add more groupings here if needed
    TRUE ~ "Other"
  )
)

group_colors <- c(
  "spliceosome" = "#990099", #purple
  "proteasome" = "#FF0000",  # red
  "sr_processome" = "#FF9999",  # pink
  "MCM_complex_2_6" = "#00FF00",  # green
  "Nuclear_pore_complex" = "#FF9900",  # orange
  "RNA_polymerase_III_complex" = "#FF00FF",  # magenta
  "COPII_vesicle_coat_complex" = "#BBBB00",  # khaki
  "MCM3" = "#009900",  # darkgreen
  "Histone_pre-mRNA_cleavage_complex" = "#0000FF",  # blue
  "Exosome" = "#00FFFF", # cyan
  "Phenylalanyl-tRNA synthetase complex" = "#993300", # brown
  "Other"  = "#000000"   # black
  # Add more group colors here if needed
)

# Create a named vector of colored labels
label_group_df <- label_group_df %>%
  mutate(label_colored = if_else(
    Group == "Other",
    OrthoDB,
    paste0("<span style='color:", group_colors[Group], ";'>", OrthoDB, "</span>")
  ))

heatmap_df_sym_labeled <- heatmap_df_sym %>%
  left_join(label_group_df, by = c("OrthoDB1" = "OrthoDB"))




x_label_map <- deframe(dplyr::select(label_group_df, OrthoDB, label_colored))

query_heatmap_df_sym_labeled <- query_heatmap_df_sym %>%
  left_join(label_group_df, by = c("OrthoDB1" = "OrthoDB"))

# Step 8: Plot as heatmap

query_orthodbpairs <-ggplot(arrange(query_heatmap_df_sym_labeled, score), 
                            aes(x = factor(query_heatmap_df_sym_labeled$OrthoDB1, ordered = TRUE, levels = pairwise_x_labels), 
                                y = factor(query_heatmap_df_sym_labeled$OrthoDB2, ordered = TRUE, levels = unique(arrange(heatmap_df_sym_labeled, score)$OrthoDB2)), 
                                fill = log10(query_heatmap_df_sym_labeled$score))) +
  geom_tile() +
  
  # Add invisible points just for the legend
  geom_point(aes(color = Group), x = Inf, y = Inf, size = 3, show.legend = TRUE) +
  
  scale_fill_viridis_c(name = "log10(SoS gene number)", limits=c(0,max(log10(query_heatmap_df_sym_labeled$score)))) +
  scale_color_manual(name = "Notable Complexes", values = group_colors) +
  scale_x_discrete(labels = x_label_map) +
  scale_y_discrete(labels = x_label_map) +# colored axis labels
  theme_minimal() +
  theme(
    axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1, size = 4.5),
    axis.text.y = ggtext::element_markdown(size = 4.5),
    legend.key.size = unit(5,"pt"),
    legend.position = "left"
  ) + labs(subtitle= j, x = "CompOG partner 1", y = "CompOG partner 2")

ggsave(file= paste0("Results/", filename,"singlepair_plot.svg"), plot=query_orthodbpairs, width=20, height=20)
