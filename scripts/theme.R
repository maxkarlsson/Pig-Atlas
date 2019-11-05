
# ----- colors & factor levels -----

spec_category_levels <- 
  c('tissue enriched',
    'group enriched', 
    'tissue enhanced', 
    'low tissue specificity',  
    'not detected')

dist_category_levels <- 
  c('detected in all', 
    'detected in many', 
    'detected in some', 
    'detected in single', 
    'not detected')

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
    "not detected " = "grey")


gene_category_pal_human_pig <- 
  gene_category_pal %>% 
  enframe() %>% 
  do(bind_rows(mutate(., name = paste(name, "human")),
               mutate(., name = paste(name, "pig")))) %$% 
  set_names(value, name) 

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