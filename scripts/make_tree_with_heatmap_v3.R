### MAKE GENUS LEVEL CIRCULAR TREE ###

pacman::p_load(ggtreeExtra, 
               ggtree,
               phyloseq, 
               dplyr, 
               ggplot2, 
               ggstar,
               treeio, 
               tidytree, 
               ggnewscale, 
               TDbook)


#############################################################
# MAKE FILE FOR TREE
#############################################################

# Load the midas database to select the relevant species
midas_db <- paste0(WD,"/data/MIDAS_53_taxa.txt")

# Select ASVs
set.seed(123)
top_20_per_WWTPID <- data_long %>% 
  filter(identy_cutoff %in% c(0.945)) %>%
  select(-total_reads) %>% 
  group_by(WWTP_ID2) %>% 
  sample_n(2) %>%  
  unnest(samples) %>%
  group_by(WWTP_ID2, Genus) %>%
  summarise(mean = mean(rel_abun_genus)) %>% 
  group_by(WWTP_ID2) %>% 
  slice_max(mean, n = 15) %>% 
  ungroup()

top_20_per_WWTPID_species <- data_long %>% 
  filter(identy_cutoff %in% c(0.945)) %>%
  select(-total_reads) %>% 
  group_by(WWTP_ID2) %>% 
  sample_n(2) %>%  
  unnest(samples) %>%
  filter(rel_abun_ASV > 0.001) %>% 
  inner_join(top_20_per_WWTPID, by = join_by(Genus, WWTP_ID2)) %>% 
  ungroup() %>% 
  distinct(Genus, OTU)

### Subset to one ASV based on genes level
top_ASV_per_genus <- data_long %>% 
  filter(identy_cutoff %in% c(0.945)) %>%
  select(-total_reads) %>% 
  group_by(WWTP_ID2) %>% 
  sample_n(2) %>%  
  unnest(samples) %>%
  filter(rel_abun_ASV > 0.001) %>% 
  inner_join(top_20_per_WWTPID, by = join_by(Genus, WWTP_ID2)) %>% 
  ungroup() %>% 
  group_by(Genus) %>% 
  slice_max(rel_abun_ASV) %>% 
  distinct(OTU, rel_abun_ASV, WWTP_ID2)


#############################################################
# READ TREE
#############################################################


# Load tree (generated on servers)
tree <- ape::read.tree(paste0(WD,"data/tree_data/20241126_MIDAS53_genera_min_0.1_mafft_output_trimmed.msa.treefile"))

# Rename the tips to coorespond with the species 
tree_names <- tree$tip.label %>% 
  str_remove("_tax_d") 
tree$tip.label <- tree_names
#plot(tree)
# Root tree by midpoint
new.tree<-phangorn::midpoint(tree)
#new.tree<-ape::root(tree, outgroup = "FLASV65.1450")



#############################################################
# LOAD MIDAS CORE  &  GUILDS
#############################################################

Nitrifiers = c("g__Nitrosomonas","g__Nitrosospira","g__Nitrospira","g__Nitrotoga")
Denitrifiers = c("g__Zoogloea", "g__Rhodoferax", "g__Thauera", "g__Rhodobacter", 
                 "g__Sulfuritalea", "g__Paracoccus", "g__Azoarcus")
PAOs <- c("g__Tetrasphaera","g__Dechloromonas","g__Ca_Accumulibacter", "g__Ca_Azonexus", "g__Azonexus",
          #"g__Ca_Propionivibrio", 
          #"g__Ca_Proximibacter", 
          "g__Ca_Lutibacillus", 
          "g__Ca_Phosphoribacter")
GAOs <- c("g__Ca_Competibacter","g__Defluviicoccus",
          "g__Propionivibrio","g__Micropruina", 
          "g__Ca_Contendobacter", "g__Kineosphaera")
Filaments <- c("g__Ca_Microthrix", "g__Leptothrix","g__Acinetobacter",
               "g__Sphaerotilus","g__Ca_Villigracilis","g__Trichococcus",
               "g__Thiothrix","g__Ca_Promineofilum","Haliscomenobacter",
               "s__Defluviicoccus_seviorii","g__Bellilinea", "g__Ca_Brachythrix",
               "g__Gordonia","g__Sarcinithrix","g__Ca_Amarolinea", "g__Ca_Amarobacter",
               "g__Kouleothrix","g__Ca_Alysiosphaera", "g__Ca_Amarofilum",
               "g__Nocardioides","g__midas_g_1668", "g__Ca_Flexicrinis", "g__Ca_Brevefilum",
               "g__Ca_Trichobacter", "g__Ca_Alysiosphaera", "g__Ca_Amarobacillus",
               "g__Anaerolinea","s__Tetrashaera_midas_s_328", "g__Ca_Defluviilinea",
               "g__Ca_Fredericiella", "g__Ca_Hadersleviella", "g__Ca_Leptovillus", 
               "g__Ca_Moduliflexus", "g__Ca_Nostocoida", "g__Ca_Pachofilum",
               "g__Ca_Promineofilum", "g__Ca_Ribeiella", "g__Ca_Sarcinithrix", 
               "g__Ca_Trichofilum", "g__Ca_Tricholinea","g__Haliscomenobacter",
               #"g__Leptolinea", "g__Nocardioides",
               "g__Leucothrix", "s__midas_s_328", "g__Mycobacterium", 
               "g__Longilinea",
               "g__midas_g_105","g__midas_g_2111","g__midas_g_344", "g__Ca_Catenibacter",
               "g__Skermania","g__Ca_Nostocoida", "g__Ca_Caldilinea", "g__Ca_Defluviifilum",
               "g__Neomegalonema","g__Beggiatoa", "g__Ca_Epilinea", "g__Ca_Flexifilum")
func <- tibble(
  Genus = c(Nitrifiers, Denitrifiers, PAOs, GAOs, Filaments),
  guild = c(rep("Nitrification", 4),
            rep("Denitrification", 7), 
            rep("PAO", 7), rep("GAO", 6), 
            rep("Filamentious", 57) 
                       )
)


#############################################################
# FORMAT DATA FOR GGTREE WITH PHYLOSEQ
#############################################################

# Filter midas database to the FLASV used for the tree
MIDAS_db_sub <- 
  read_csv2(file = midas_db, 
            col_names = F) %>% 
  select(X1,X6, X7, X5) %>% 
  rename("FLASV" = 1,"Species" = 3,
         "Genus" = 2, "Family" = 4) %>% 
  mutate(FLASV = str_remove(FLASV, "\tk__Bacteria")) %>% 
  filter(FLASV %in% tree_names) %>% 
  mutate(lab = paste0(Family, ";", Genus, ";",Species)) %>% 
  arrange(FLASV)


## Subset to the ASV/species in the tree
 ## Mean for each WWTP ID 
abun_merged_genus <- 
  data_long %>% 
  filter(!WWTP_ID2  %in% c("US", "PT", "IL", "PT", "NL(A)", "DK(B)")) %>% 
  mutate(WWTP_ID2  = if_else(WWTP_ID2  == "NL(B)", "NL", WWTP_ID2 ), 
         WWTP_ID2  = if_else(WWTP_ID2  == "DK(A)", "DK", WWTP_ID2 )) %>% 
  filter(ProjectName != "MIDAS4") %>% 
  filter(!is.na(WWTP_ID)) %>% 
  mutate(x = " ") %>% 
  mutate(samples = 
           map(.x = samples, ~
                 ungroup(.) %>% 
                 filter(Genus %in% unlist(unique(top_20_per_WWTPID_species$Genus))) %>% 
                 group_by(Genus) %>%
                 summarise(sum = sum(rel_abun_ASV)) 
           )) %>% 
  unnest(samples)  %>% 
  left_join(top_ASV_per_genus %>% select(Genus, OTU)) %>% 
  group_by(WWTP_ID2, Genus, OTU) %>%
  summarise(mean_abun_WWTP_ID = mean(sum), 
            .groups = "drop") %>% 
  filter(mean_abun_WWTP_ID > 0.01) %>% 
  pivot_wider(names_from = WWTP_ID2, values_from = mean_abun_WWTP_ID, values_fill = 0.0000001) %>% 
  column_to_rownames("OTU") %>% select(-Genus)

tree_subset <- d_norm %>% 
  #amp_subset_taxa(tax_vector = tree_names)
  amp_subset_taxa(tax_vector = unlist(top_ASV_per_genus$OTU))


OTU <- otu_table(abun_merged_genus, taxa_are_rows = TRUE)
TAX = tax_table(tree_subset$tax %>% select(-OTU) %>%
                  mutate(
                    # type = case_when(
                    #   Genus %in% unlist(Genus_core_midas5.3 %>% filter(Genus_type == "SC") %>% select(Genus)) ~ "Strict/General core",
                    #   Genus %in% unlist(Genus_core_midas5.3 %>% filter(Genus_type == "GC") %>% select(Genus)) ~ "Strict/General core",
                    #   Genus %in% unlist(Genus_core_midas5.3 %>% filter(Genus_type == "LC") %>% select(Genus)) ~ "Loose core/CRAT",
                    #   Genus %in% unlist(Genus_core_midas5.3 %>% filter(Genus_type == "CAT") %>% select(Genus)) ~ "Loose core/CRAT",
                    #   .default = "No core group"),
                    functional =  case_when(
                      Genus %in% unlist(func %>% filter(guild == "Nitrification") %>% select(Genus)) ~ "Nitrification",
                      Genus %in% unlist(func %>% filter(guild == "Denitrification") %>% select(Genus)) ~ "Denitrification",
                      Genus %in% unlist(func %>% filter(guild == "GAO") %>% select(Genus)) ~ "GAO",
                      Genus %in% unlist(func %>% filter(guild == "PAO") %>% select(Genus)) ~ "PAO",
                      Genus %in% unlist(func %>% filter(guild == "Filamentious") %>% select(Genus)) ~ "Filamentious",
                      .default = "Unknown/No function"),
                    ) %>% as.matrix())
SAMPLE_DATA <- sample_data(
  data_long %>% 
    filter(!WWTP_ID2  %in% c("US", "PT", "IL", "PT", "NL(A)", "DK(B)")) %>% 
    mutate(WWTP_ID2  = if_else(WWTP_ID2  == "NL(B)", "NL", WWTP_ID2 ), 
           WWTP_ID2  = if_else(WWTP_ID2  == "DK(A)", "DK", WWTP_ID2 )) %>% 
    filter(!is.na(WWTP_ID)) %>%
    filter(ProjectName != "MIDAS4") %>% 
    distinct(WWTP_ID2, country, continent, country_code) %>% 
    mutate(WWTP_ID_1 = WWTP_ID2, 
           continent = factor(continent), 
           WWTP_ID2 = fct_reorder(WWTP_ID2, as.numeric(continent)),
           ekx = "x",
           ) %>%
    column_to_rownames("WWTP_ID_1"))

phylo_object <- 
  phyloseq(OTU, TAX, new.tree, SAMPLE_DATA
  ) 


colors_vector <- c("#D65B44", "#5366A9", "#88B661", "#8F5C8C", "#DABD5F", "#78A39C", "#B852B8", "#B88454",  "#D18EAF", "#5BB6B6")


p1 <- ggtree(phylo_object, layout="fan",
            open.angle=7, size = 0.5)



#p1
p2 <- p1 +
  scale_fill_manual(values = c("#66C2A5", "#8DA0CB",  "#D18EAF","#E5C494", "#B3B3B3", "white"), 
                    guide = "none",
                      # guide_legend(title = "Funtionality", 
                      #                    override.aes = list(
                      #                      alpha = 1
                      #                      #shape = 22, 
                      #                                        #color = c("#66C2A5", "#8DA0CB",  "#E5C494", "#B3B3B3", "white")
                      #                                        ))
                    ) +
  new_scale_fill() +
  new_scale_color() + 
  geom_tippoint(mapping=aes(fill=str_sub(Phylum,4)
                            #,shape = type,
                            ), 
                size=4, shape = 21,
                show.legend=T) +
  scale_fill_manual(values = colors_vector, 
                    guide = "none"
                      # guide_legend(
                      # override.aes = list(size = 4), 
                      #                                            title = "Phylum"
                      # )
                    ) +
  # scale_shape_manual(values = c(21,22,25), 
  #                    guide = guide_legend(title = "Type", 
  #                                         override.aes = list(fill = c("grey40", "grey40", "grey40"), size = 4)
  #                    )) + 
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) 
#p2
p3 <- p2 +
  new_scale_fill() +
  new_scale_color() +
  geom_fruit(
    #data = subset_samples(phylo_object, continent!="Africa"),
    geom=geom_tile,
             mapping=aes(y=label, 
                         x=fct_reorder(paste0((WWTP_ID2)), as.numeric(continent)),
                         #color = functional,
                         #alpha=Abundance, 
                         fill= cut(
                           Abundance,
                           breaks = c(0, 0.01, 0.1, 1, 5, Inf),  # Define bins (upper bound inclusive)
                           labels = c("0", "0.01", "0.1", "1", ">5"),
                           include.lowest = TRUE
                         ), 
                         #group = Abundance
             ),
             color = "white",#"grey30", 
             offset = 0.07, 
             width = 0.10,  
    size = 1, 
    pwidth=1, 
    linetype = 1, 
             axis.params=list(
               axis = "x",
               text.size = 5.5,  # Adjusted for readability
               #hjust = "outward",    # Centered labels
               #vjust = "outward",    # Moved slightly away from the heatmap
               hjust = 0.8, 
               vjust = 0.1,
               text.angle = 45, # Changed to a diagonal angle for better spacing
               line.size = 10,
               #nudge_x = 5, 
               #nudge_y = 0,
               line.color = "white",
               colour = "black"
             )
  ) +
  scale_fill_manual(
    values = c(
      "0" = "white", 
      "0.01" = "#FDEDEC", 
      "0.1" = "#F5B7B1", 
      "1" = "#D98880", 
      ">5" = "#922B21"
    ), 
    breaks = c(
      "0-0.01", 
      "0.01-0.1", 
      "0.1-1", 
      "1-5", 
      ">5" 
    ), 
    guide = "none",
    na.value = "grey"
  )
p3 
p4 <- p3 +
  geom_tiplab(aes(label = case_when(str_detect(Genus, "midas") ~ str_replace(str_sub(Genus, 4), "_", "_"),
                                    !str_detect(Genus, "midas") & !str_detect(Genus, "Ca_") ~ paste0("italic('",str_replace_all(str_sub(Genus, 4), "_", " "), "')"),
                                    !str_detect(Genus, "midas") & str_detect(Genus, "Ca_") ~ paste0("italic('Ca.')~",str_remove(str_sub(Genus, 4), "Ca_")),
                                    .default = str_sub(Genus, 4)),
                  #fill = functional #str_sub(Phylum,4), color = functional
                  ),
              color = "black", #"black", # color for label font
              geom = "label",  # labels not text
              label.padding = unit(0.05, "lines"), # amount of padding around the labels
              label.size = 0, parse = TRUE,
              #fill = "white",
              alpha = 1, 
              check.overlap = T,
              size = 6, show.legend = T, 
              align = T, 
              hjust = 0, 
              vjust = 0.5, 
              linesize = 0.1,linetype = "79", #linecolor = "grey",
              offset = 1.65   ## HOW FAR ARE THE LABELS FROM THE CENTER
              ) 
p5 <- p4 +
  new_scale_fill() +
  new_scale_color() +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(y=label, 
                x=ekx, 
                fill=str_sub(Phylum,4), 
    ),
    color = "white",#alpha = 0.7,
    offset = 0.015, #0.03, 
    width = 0.04, show.legend = FALSE,
    size = 0.1,  pwidth=0.001, linetype = 0,
  ) + 
  scale_fill_manual(values = colors_vector, 
                    guide = "none" 
                    # guide_legend(title = "Phylum"#, override.aes = list(shape = 22, alpha = 1)
                    #                    )
  )


ggsave(plot = p5, filename = paste0(OutputPath, "20250425_tree_top15.png"), width = 21, height = 21, dpi = 800)  
ggsave(plot = p5, paste0(OutputPath, "20250425_tree_top15.pdf"), 
       limitsize = F,  useDingbats=FALSE, 
       width = 21, height = 21)



####################
# Make legends 

## Tip labels functionallity 
ggtree(phylo_object, layout="fan", 
       open.angle=8, size = 0.1) + 
  geom_tiplab(aes(label = case_when(str_detect(Genus, "midas") ~ str_replace(str_sub(Genus, 4), "_", "_"),
                                    !str_detect(Genus, "midas") & !str_detect(Genus, "Ca_") ~ paste0("italic('",str_replace_all(str_sub(Genus, 4), "_", " "), "')"),
                                    !str_detect(Genus, "midas") & str_detect(Genus, "Ca_") ~ paste0("italic('Ca.')~",str_remove(str_sub(Genus, 4), "Ca_")),
                                    .default = str_sub(Genus, 4)),
                  fill = functional #str_sub(Phylum,4), color = functional
  ),
  color = "black", #"black", # color for label font
  geom = "label",  # labels not text
  label.padding = unit(0.15, "lines"), # amount of padding around the labels
  label.size = 0, parse = TRUE,
  #fill = "white",
  alpha = 1, 
  check.overlap = T,
  size = 4.5, show.legend = T, 
  align = T, hjust = 0, vjust = 0.5, linesize = 0.1,linetype = "79", #linecolor = "grey",
  offset = 0.88#1.13
  ) +
  scale_fill_manual(values = c("#66C2A5", "#8DA0CB",  "#D18EAF","#E5C494", "#B3B3B3", "white"), 
                    guide = #"none",
                    guide_legend(title = "Function", 
                                       override.aes = list(
                                         alpha = 1, 
                                         shape = 22, #color = "black"
                                                           color = c("#66C2A5", "#8DA0CB",  "#D18EAF","#E5C494", "#B3B3B3", "white")
                                                           ))
  ) 
ggsave(paste0(OutputPath, "20241205_function_label.png"), width = 16.5, height = 14)  


phylum_lab <- ggtree(phylo_object, layout="fan", 
       open.angle=8, size = 0.1) + 
  geom_tippoint(mapping=aes(color =str_sub(Phylum,4),
                            #shape = type,
  ), 
  size=3, shape = 15,
  show.legend=T) +
  scale_color_manual(values = colors_vector, guide = guide_legend(override.aes = list(size = 4), 
                                                                 title = "Phylum")) +
  scale_shape_manual(values = c(21,22,25), 
                     guide = guide_legend(title = "Type", 
                                          override.aes = list(fill = c("grey40", "grey40", "grey40"), size = 4)
                     )) + 
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))


# Extract the legend. Returns a gtable
leg_phylum_lab <- ggpubr::get_legend(phylum_lab)
# Convert to a ggplot and print
ggpubr::as_ggplot(leg_phylum_lab)
ggsave(plot = leg_phylum_lab, filename = paste0(OutputPath, "20250425_phylum_label.png"), width = 4, height = 6)  
ggsave(plot = leg_phylum_lab, paste0(OutputPath, "20250425_phylum_label.pdf"), 
       limitsize = F,  useDingbats=FALSE, 
       width = 4, height = 6)


###### ABUNDANCE
abundace_leg <- ggtree(phylo_object, layout="fan", 
       open.angle=30, size = 0.1) + 
    geom_fruit(
    #data = subset_samples(phylo_object, continent!="Africa"),
    geom=geom_tile,
    mapping=aes(y=label, 
                x=fct_reorder(paste0((WWTP_ID2)), as.numeric(continent)),
                #color = functional,
                #alpha=Abundance, 
                fill= cut(
                  Abundance,
                  breaks = c(0, 0.01, 0.1, 1, 5, Inf),  # Define bins (upper bound inclusive)
                  labels = c("0", "0.01", "0.1", "1", ">5"),
                  include.lowest = TRUE
                ), 
                #group = Abundance
    ),
    color = "white",#"grey30", 
    offset = 0.03, #0.03, 
    #width = 0.8,
    size = 0.8,  pwidth=1, linetype = 1, 
    axis.params=list(
      #geom = "label",
      axis       = "x",
      text.size  = 7,
      hjust      = 0,
      vjust      = 0.5, #nbreak = 10,
      text.angle = -90,
      line.size = 2, line.color = "white",
      colour = "black",
      #c(rep("black",1), rep("#FFB2B2", 1),rep("#FFFFCC", 5), rep("#BFD8FF", 7), rep("#B2FFB2", 2))
    )
  ) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  scale_fill_manual(
    values = c(
      #"white", "#FDEDEC",  "#F5B7B1", "#D98880", "#922B21")
      "0" = "white", 
      "0.01" = "#FDEDEC", 
      "0.1" = "#F5B7B1", 
      "1" = "#D98880", 
      ">5" = "#922B21"
    ), 
    labels = c(
      "0" = "0-0.01", 
      "0.01"="0.01-0.1", 
      "0.1" = "0.1-1", 
      "1" = "1-5", 
      ">5" 
    ), 
    guide = guide_legend(#color = "grey", fill = "grey", 
                         title = "Abundance [%]",
                              #                            theme = theme(legend.key = element_rect(fill = "transparent", color = "black"),
                              #                                          legend.background =  element_rect(fill = "transparent"),
                              #                                          legend.text = element_text(color = "black")),
    na.value = "grey"
  ))


# Extract the legend. Returns a gtable
abundace_leg_lab <- ggpubr::get_legend(abundace_leg)
# Convert to a ggplot and print
ggpubr::as_ggplot(abundace_leg_lab)
ggsave(plot = abundace_leg_lab, filename = paste0(OutputPath, "20250425_abundance_label.png"), width = 4, height = 6)  
ggsave(plot = abundace_leg_lab, filename = paste0(OutputPath, "20250425_abundance_label.pdf"), 
       limitsize = F,  useDingbats=FALSE, 
       width = 4, height = 6)

