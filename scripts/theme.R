
# ----- colors & factor levels -----

spec_category_levels <- 
  c('tissue enriched',
    'group enriched', 
    'tissue enhanced', 
    'low tissue specificity',  
    'not detected',
    
    'Tissue enriched',
    'Group enriched', 
    'Tissue enhanced', 
    'Low tissue specificity',  
    'Not detected')

dist_category_levels <- 
  c('detected in all', 
    'detected in many', 
    'detected in some', 
    'detected in single', 
    'not detected',
    
    'Detected in all', 
    'Detected in many', 
    'Detected in some', 
    'Detected in single', 
    'Not detected')

enrichment_overlap_levels <- 
  c("full overlap", 
    "partial overlap", 
    "no overlap",
    
    "Full overlap", 
    "Partial overlap", 
    "No overlap")

shared_category_levels <- 
  c("shared", 
    "minor difference", 
    "medium difference",
    "major difference")

enrichment_overlap_pal <-
  set_names(c(viridis(3), inferno(4),
              viridis(3), inferno(4)),
            c("full overlap", 
              "partial overlap", 
              "no overlap", 
              "shared", 
              "minor difference", 
              "medium difference",
              "major difference",
              
              "Full overlap", 
              "Partial overlap", 
              "No overlap", 
              "Shared", 
              "Minor difference", 
              "Medium difference",
              "Major difference"))


spec_category_overlap_levels <-
  c('tissue enriched full overlap', 
    'tissue enriched partial overlap', 
    'tissue enriched no overlap', 
    'group enriched full overlap',
    'group enriched partial overlap', 
    'group enriched no overlap', 
    'tissue enhanced full overlap',
    'tissue enhanced partial overlap', 
    'tissue enhanced no overlap', 
    'low tissue specificity no overlap', 
    'not detected no overlap')

spec_category_overlap_levels_short <-
  c('tissue enriched FO', 
    'tissue enriched PO', 
    'tissue enriched NO', 
    'group enriched FO',
    'group enriched PO', 
    'group enriched NO', 
    'tissue enhanced FO',
    'tissue enhanced PO', 
    'tissue enhanced NO', 
    'low tissue specificity NO', 
    'not detected NO')

spec_category_overlap_pal <- 
  c('tissue enriched full overlap' = "#E41A1C", 
    'tissue enriched partial overlap' = "#ED6667" , 
    'tissue enriched no overlap' = "#F6B2B3" , 
    'group enriched full overlap' = "#FF9D00",
    'group enriched partial overlap' = "#FFBD55", 
    'group enriched no overlap' = "#FFDEAA", 
    'tissue enhanced full overlap' = "#984EA3",
    'tissue enhanced partial overlap' = "#BA89C1", 
    'tissue enhanced no overlap' = "#DCC3E0", 
    'low tissue specificity no overlap' = "grey40", 
    'not detected no overlap' = "grey")

spec_category_overlap_short_pal <- 
  c('tissue enriched FO' = "#E41A1C", 
    'tissue enriched PO' = "#ED6667" , 
    'tissue enriched NO' = "#F6B2B3" , 
    'group enriched FO' = "#FF9D00",
    'group enriched PO' = "#FFBD55", 
    'group enriched NO' = "#FFDEAA", 
    'tissue enhanced FO' = "#984EA3",
    'tissue enhanced PO' = "#BA89C1", 
    'tissue enhanced NO' = "#DCC3E0", 
    'low tissue specificity NO' = "grey40", 
    'not detected NO' = "grey")


gene_category_pal <- 
  c("tissue enriched" = "#e41a1c",
    "group enriched" = "#FF9D00",
    "tissue enhanced" = "#984ea3",
    "low tissue specificity" = "grey40",
    
    "detected in all" = "#253494",
    "detected in many" = "#2c7fb8",
    "detected in some" = "#41b6c4",
    "detected in single" = "#a1dab4",
    
    "not detected" = "grey", 
    "not detected " = "grey", 
    
    "Tissue enriched" = "#e41a1c",
    "Group enriched" = "#FF9D00",
    "Tissue enhanced" = "#984ea3",
    "Low tissue specificity" = "grey40",
    
    "Detected in all" = "#253494",
    "Detected in many" = "#2c7fb8",
    "Detected in some" = "#41b6c4",
    "Detected in single" = "#a1dab4",
    
    "Not detected" = "grey", 
    "Not detected " = "grey")


gene_category_pal_human_pig <- 
  gene_category_pal %>% 
  enframe() %>% 
  do(bind_rows(mutate(., name = paste(name, "human")),
               mutate(., name = paste(name, "pig")))) %$% 
  set_names(value, name) 

gene_category_pal_comparison <- 
  gene_category_pal %>% 
  enframe() %>% 
  do(bind_rows(mutate(., name = paste(name, "canon")),
               mutate(., name = paste(name, "comparison")))) %$% 
  set_names(value, name) 

anova_pal <- 
  c("Residuals and species" = "gray",
    "Residuals" = "gray", 
    "tissue_ID" = "#D95E37", 
    "comparison_tissue" = "#D95E37", 
    "Tissue" = "#D95E37", 
    "individual" = "#149684", 
    "species" = "#F9C770",
    "Species" = "#F9C770",
    "sex" = "#E8C9C7")

anova_levels <- 
  c("Residuals and species",
    "Residuals", 
    "tissue_ID", 
    "comparison_tissue", 
    "Tissue",
    "individual", 
    "species",
    "Species",
    "sex")

sex_palette <- 
  c("male" = "#D96767", "female" = "#214785")

elevation_identity_levels <- 
  c("identity", "overlapping tissues", "different tissues", "no enrichment")

elevation_identity_pal <- 
  set_names(rev(pals::ocean.haline(9)[c(3, 4, 6, 7)]), elevation_identity_levels)

elevation_overlap_levels <- 
  c("Pig", "Both", "Human")

elevation_overlap_pal <- 
  c(set_names(rev(pals::ocean.haline(9)[c(3, 4, 6)]), elevation_overlap_levels),
    set_names(rev(pals::ocean.haline(9)[c(3, 4, 6)]), c("Elevated in Pig", "Elevated in both", "Elevated in Human")))

elevation_overlap_pal2 <-
  c("Human" = "#DB1F48", 
    "Overlap" = "#004369", 
    "Pig" = "#01949A",
    "Elevated in Human" = "#DB1F48", 
    "Elevated in both" = "#004369", 
    "Elevated in Pig" = "#01949A",
    "Enriched in Human" = "#DB1F48", 
    "Enriched in both" = "#004369", 
    "Enriched in Pig" = "#01949A")

# protein_type_pal <- 
#   c("secreted" = '#911D51',
#     "membrane" = '#6D4BAA', 
#     "other" = '#008490', 
#     "cd_marker" = '#318F1E', 
#     "transcription_factors" = '#B8801B', 
#     "mitochondrial" = '#E371B4', 
#     "ribosomal" = '#89A0F3', 
#     "none" = "black",
#     '#00C9BC', '#97C542', '#FFA05E')
# 
# protein.localization.palette2 <- c("membrane" = "#CE70A4",
#                                    "secreted" = 	"#FCAC3B",
#                                    "membrane and secreted isoforms" = "#755A85")

# protein.localization.palette <- c("intracellular and membrane isoforms" = "#858141",
#                                   "membrane" = "#6DB9C6",
#                                   "intracellular" = "#FCAC3B",
#                                   "secreted" = 	"#CE70A4",
#                                   "intracellular and secreted isoforms" = "#CF5734",
#                                   "membrane and secreted isoforms" = "#755A85",
#                                   "intracellular, membrane, secreted isoforms" = "#794A39")

# ----- themes -----

simple_theme <- 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

heatmap_palette = viridis::inferno(20, direction = -1)



theme_stripped <-
  theme(panel.background = element_rect(fill = NA, colour = NA),
        plot.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.key = element_rect(colour = NA),
        #legend.position = "bottom",
        #legend.direction = "horizontal",
        legend.key.size= unit(0.3, "cm"),
        legend.title = element_text(face="italic"),
        axis.line = element_line(colour="black",size=0.5))

theme_stripped_frame <-
  theme(panel.background = element_rect(fill = NA, colour = "gray"),
        plot.background = element_rect(fill = NA, color = "gray"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.key = element_rect(colour = NA),
        #legend.position = "bottom",
        #legend.direction = "horizontal",
        legend.key.size= unit(0.3, "cm"),
        legend.title = element_text(face="italic"),
        axis.line = element_line(colour="black",size=0.5))

theme_angletext <- theme(axis.text.x = element_text(angle = 60, hjust = 1))

# Make plot theme
stripped_theme <-
  theme(panel.background = element_rect(fill = NA, colour = NA),
        plot.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.key = element_rect(colour = NA),
        #legend.position = "bottom",
        #legend.direction = "horizontal",
        legend.key.size= unit(0.3, "cm"),
        legend.title = element_text(face="italic"),
        axis.line = element_line(colour="black",size=0.5))



# stripped theme facet
stripped_theme_facet <-
  stripped_theme+
  theme(legend.position = "right",
        panel.border = element_rect(color = "gray", fill = NA),
        strip.background = element_rect(fill = "#003253",
                                        color = "#003253"),
        strip.text = element_text(color = "white"))

# stripped theme facet
stripped_theme_HPA <-
  stripped_theme+
  theme(legend.position = "right",
        panel.border = element_rect(color = "gray", fill = NA),
        strip.background = element_rect(fill = "#313131",
                                        color = "#313131"),
        strip.text = element_text(color = "white"), 
        plot.background = element_rect(fill = "#D9D9D9",
                                       color = "#D9D9D9"), 
        panel.background = element_rect(fill = "white",
                                        color = "black"))


