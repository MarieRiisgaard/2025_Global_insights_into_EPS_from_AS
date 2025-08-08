# ADONIS from long dataframe

ADONIS_from_long_df <- function(
    data = data_merged_all2,
    tax_level = "Species", 
    abun_limit = 0.1,
    test_varibles = c("Investigator") 
){
  
  if (tax_level == "Species") {
    df_tax_level <- {{data}} %>% 
      mutate(samples = 
               map(.x = samples, ~
                     ungroup(.) %>% 
                     distinct(Species, rel_abun_species) %>% 
                     filter(!str_detect(Species, "unclassified")) %>% 
                     mutate(over_limit = if_else(rel_abun_species >= abun_limit, T, F))
               )) %>% 
      unnest(samples) %>% 
      rename(
        "tax_level" = "Species",
        "tax_abun" = "rel_abun_species"
      )
  } else if (tax_level == "Genus") {
    df_tax_level <- {{data}} %>% 
      mutate(samples = 
               map(.x = samples, ~
                     ungroup(.) %>% 
                     distinct(Genus, rel_abun_genus) %>% 
                     filter(!str_detect(Genus, "unclassified")) %>% 
                     mutate(over_limit = if_else(rel_abun_genus >= abun_limit, T, F))
               )) %>% 
      unnest(samples) %>% 
      rename(
        "tax_level" = "Genus",
        "tax_abun" = "rel_abun_genus"
      )
  } else if (tax_level == "OTU") {
    df_tax_level <- {{data}} %>% 
      mutate(samples = 
               map(.x = samples, ~
                     ungroup(.) %>% 
                     distinct(OTU, rel_abun_ASV) %>% 
                     #filter(!str_detect(Genus, "unclassified")) %>% 
                     mutate(over_limit = if_else(rel_abun_ASV >= abun_limit, T, F))
               )) %>% 
      unnest(samples) %>% 
      rename(
        "tax_level" = "OTU",
        "tax_abun" = "rel_abun_ASV"
      )
  } else if (tax_level == "Family") {
    
    df_tax_level <- {{data}} %>% 
      ungroup() %>% 
      mutate(samples = map(.x = samples, ~group_by(.x, Family) %>% mutate(rel_abun_family = sum(rel_abun_ASV)) %>% ungroup())) %>% 
      mutate(samples = 
               map(.x = samples, ~
                     ungroup(.) %>% 
                     distinct(Family, rel_abun_family) %>% 
                     filter(!str_detect(Family, "unclassified")) %>% 
                     mutate(over_limit = if_else(rel_abun_family >= abun_limit, T, F))
               )) %>% 
      unnest(samples) %>% 
      rename(
        "tax_level" = "Family",
        "tax_abun" = "rel_abun_family"
      )
  } else if (tax_level == "Order") {
    
    df_tax_level <- {{data}} %>% 
      ungroup() %>% 
      mutate(samples = map(.x = samples, ~group_by(.x, Order) %>% mutate(rel_abun_order = sum(rel_abun_ASV)) %>% ungroup())) %>% 
      mutate(samples = 
               map(.x = samples, ~
                     ungroup(.) %>% 
                     distinct(Order, rel_abun_order) %>% 
                     filter(!str_detect(Order, "unclassified")) %>% 
                     mutate(over_limit = if_else(rel_abun_order >= abun_limit, T, F))
               )) %>% 
      unnest(samples) %>% 
      rename(
        "tax_level" = "Order",
        "tax_abun" = "rel_abun_order"
      )
  } else if (tax_level == "Class") {
    
    df_tax_level <- {{data}} %>% 
      ungroup() %>% 
      mutate(samples = map(.x = samples, ~group_by(.x, Class) %>% mutate(rel_abun_class = sum(rel_abun_ASV)) %>% ungroup())) %>% 
      mutate(samples = 
               map(.x = samples, ~
                     ungroup(.) %>% 
                     distinct(Class, rel_abun_class) %>% 
                     filter(!str_detect(Class, "unclassified")) %>% 
                     mutate(over_limit = if_else(rel_abun_class >= abun_limit, T, F))
               )) %>% 
      unnest(samples) %>% 
      rename(
        "tax_level" = "Class",
        "tax_abun" = "rel_abun_class"
      )
  } else if (tax_level == "Phylum") {
    
    df_tax_level <- {{data}} %>% 
      ungroup() %>% 
      mutate(samples = map(.x = samples, ~group_by(.x, Phylum) %>% mutate(rel_abun_phylum = sum(rel_abun_ASV)) %>% ungroup())) %>% 
      mutate(samples = 
               map(.x = samples, ~
                     ungroup(.) %>% 
                     distinct(Phylum, rel_abun_phylum) %>% 
                     filter(!str_detect(Phylum, "unclassified")) %>% 
                     mutate(over_limit = if_else(rel_abun_phylum >= abun_limit, T, F))
               )) %>% 
      unnest(samples) %>% 
      rename(
        "tax_level" = "Phylum",
        "tax_abun" = "rel_abun_phylum"
      )
  } else {
    print("Check spelling of tax_level ")
  }
  
  ## Filter low abundant Species
  tax_over_limit <- 
    df_tax_level %>% 
    group_by(tax_level) %>% 
    reframe(sum = sum(over_limit)) %>% 
    filter(sum > 0) %>% 
    distinct(tax_level) %>% unlist()
  
  df_tax_level_filtered <- 
    df_tax_level %>% 
    filter(tax_level %in% tax_over_limit)
  
  ## Calculate distance matrix 
  abun_matrix <- df_tax_level_filtered %>% 
    distinct(SampleID, tax_level, tax_abun) %>% 
    pivot_wider(names_from = tax_level, values_from = tax_abun, values_fill = 0) %>% 
    column_to_rownames("SampleID") %>% 
    as.matrix()
  
  bray_curtis_dist <- vegan::vegdist((abun_matrix), method = "bray")
  
  
  # 
  # bc.dist.matrix <- 
  #   vegan::vegdist(t(ampvis_object_fil$abund), method = distance_matrix)
  # 
  
  adonis_of_differnet_varibles <- function(varibles_character){
    
    metadata <- 
      data %>% 
      rename("test_var" = varibles_character) %>% 
      select(SampleID, test_var)
    
    var <- if_else(is.numeric(metadata$test_var),
                   paste0("in the range: ", min(as.numeric(unique(metadata$test_var)), na.rm = T), "-", 
                          max(as.numeric(metadata$test_var))),
                   paste0(unique(as.character(metadata$test_var)),collapse = ", "))
    
    x <- vegan::adonis2(bray_curtis_dist ~ test_var, data = metadata)
    
    #a = cat("Test of **" varibles_character,"** R^2=", round(x$R2[1],2), ", p<", x$`Pr(>F)`[1])
    b <- (paste0("<br>",tax_level," ADONIS analysis for **", varibles_character, "**. <br>", 
                 "**R^2^=", round(x$R2[1],2), ", p<", x$`Pr(>F)`[1], "**. <br> ", 
                 "The varibles tested are: *", 
                 var, "*.<br>"))
    
    d <- tibble(
      test = varibles_character, 
      result = b
    )
    d
  }
  
  list <- lapply(test_varibles, 
                 FUN = adonis_of_differnet_varibles)
  
  do.call(rbind, list)
  
}








