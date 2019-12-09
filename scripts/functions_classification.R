
# expression_col = "nx"
# tissue_col = "consensus_tissue_name" 
# gene_col = "enssscg_id"  
# enr_fold = 4
# max_group_n = 5
# det_lim = 1
# data <- pig_atlas_consensus

hpa_gene_classification <- 
  function(data, expression_col, tissue_col, gene_col, enr_fold, max_group_n, det_lim = 1) {
    data_ <- 
      data %>% 
      select(gene = gene_col,
             expression = expression_col,
             tissue = tissue_col) %>% 
      mutate(expression = round(expression, 4)) 
    
    if(any(is.na(data_$expression))) stop("NAs in expression column")
    if(any(is.na(data_$gene))) stop("NAs in gene column")
    if(any(is.na(data_$tissue))) stop("NAs in tissue column")
    
    n_groups <- length(unique(data_$tissue))
  
    gene_class_info <- 
      data_ %>%
      group_by(gene) %>%
      summarise(
        
        # Gene expression distribution metrics
        mean_exp = mean(expression, na.rm = T),
        min_exp = min(expression, na.rm = T),
        max_exp = max(expression, na.rm = T), 
        max_2nd = sort(expression)[length(expression)-1],
        
        # Expression frequency metrics
        n_exp = length(which(expression >= det_lim)),
        frac_exp = n_exp/length(expression[!is.na(expression)])*100,
        
        # Limit of enhancement metrics
        lim = max_exp/enr_fold, 
        
        exps_over_lim = list(expression[which(expression >= lim & expression >= det_lim)]),
        n_over = length(exps_over_lim[[1]]), 
        mean_over = mean(exps_over_lim[[1]]),
        min_over = ifelse(n_over == 0, NA,
                          min(exps_over_lim[[1]])),
        
        max_under_lim = max(expression[which(expression < min_over)], det_lim*0.1),
        
        
        exps_enhanced = list(which(expression/mean_exp >= enr_fold & expression >= det_lim)),
        
        
        
        
        # Expression patterns
        enrichment_group = paste(sort(tissue[which(expression >= lim & expression >= det_lim)]), collapse=";"),
        
        n_enriched = length(tissue[which(expression >= lim & expression >= det_lim)]),
        n_enhanced = length(exps_enhanced[[1]]), 
        enhanced_in = paste(sort(tissue[exps_enhanced[[1]]]), collapse=";"),
        n_na = n_groups - length(expression),
        max_2nd_or_lim = max(max_2nd, det_lim*0.1),
        tissues_not_detected = paste(sort(tissue[which(expression < det_lim)]), collapse=";"),
        tissues_detected = paste(sort(tissue[which(expression >= det_lim)]), collapse=";")) 
      
    
    gene_categories <- 
      gene_class_info %>%
      
      mutate(
        spec_category = case_when(n_exp == 0 ~ "not detected", 
                                  
                                  # Genes with expression fold times more than anything else are tissue enriched
                                  max_exp/max_2nd_or_lim >= enr_fold ~ "tissue enriched", 
                                  
                                  # Genes with expression fold times more than other tissues in groups of max group_n - 1 are group enriched
                                  max_exp >= lim &
                                    n_over <= max_group_n & n_over > 1 &
                                    mean_over/max_under_lim >= enr_fold ~ "group enriched", 
                                  
                                  # Genes with expression in tissues fold times more than the mean are tissue enhance
                                  n_enhanced > 0 ~ "tissue enhanced", 
                                  
                                  # Genes expressed with low tissue specificity
                                  T ~ "low tissue specificity"), 
        
        
        dist_category = case_when(frac_exp == 100 ~ "detected in all",
                                  frac_exp >= 31 ~ "detected in many",
                                  n_exp > 1 ~ "detected in some",
                                  n_exp == 1 ~ "detected in single",
                                  n_exp == 0 ~ "not detected"),
        
        spec_score = case_when(spec_category == "tissue enriched" ~ max_exp/max_2nd_or_lim,
                               spec_category == "group enriched" ~ mean_over/max_under_lim, 
                               spec_category == "tissue enhanced" ~ max_exp/mean_exp)) 
    
      
    
    
    ##### Rename and format
    gene_categories %>%
      mutate(enriched_tissues = case_when(spec_category %in% c("tissue enriched", "group enriched") ~ enrichment_group,
                                          spec_category == "tissue enhanced" ~ enhanced_in),
             n_enriched = case_when(spec_category %in% c("tissue enriched", "group enriched") ~ n_enriched,
                                    spec_category == "tissue enhanced" ~ n_enhanced)) %>%
      select(gene, 
             spec_category, 
             dist_category, 
             spec_score,
             n_expressed = n_exp, 
             fraction_expressed = frac_exp,
             max_exp = max_exp,
             enriched_tissues,
             n_enriched,
             n_na = n_na,
             tissues_not_detected,
             tissues_detected) 
      
    
    
  }	



calc_gene_correlations <- 
  function(data, var1, var2, val1, val2, cor_method = "spearman", p_adjust_method = "BH", alternative = "two.sided") {
    
    data_ <- 
      data %>%
      rename(var1 = var1, 
             var2 = var2, 
             val1 = val1, 
             val2 = val2)
      
    
    data_ %>%
      group_by(var1, var2) %>% 
      do(if(cor_method == "pearson") {
        cor.test(.$val1, .$val2, method = cor_method, alternative = alternative) %$%
          tibble(pval = p.value, 
                 cor = estimate, 
                 lo_confint = conf.int[1], 
                 hi_confint = conf.int[2])
      } else if(cor_method == "spearman") {
        cor.test(.$val1, .$val2, method = cor_method, alternative = alternative) %$%
          tibble(pval = p.value, 
                 cor = estimate)
        }) %>% 
      ungroup() %>%
      mutate(padj = p.adjust(pval, method = p_adjust_method), 
             significant = padj <= 0.05, 
             log10P = -log10(padj)) %>% 
      arrange(padj) %>% 
      set_colnames(c(var1, var2, colnames(.)[-c(1, 2)]))
  }

calc_gene_distance <- 
  function(data, var1, var2, val1, val2) {
    
    data_ <- 
      data %>%
      rename(var1 = var1, 
             var2 = var2, 
             val1 = val1, 
             val2 = val2)
    
    
    data_ %>%
      group_by(var1, var2) %>% 
      summarise(dist = rbind(val1, val2) %>% 
                  dist() %>% 
                  as.numeric(), 
                mean_var1 = mean(val1), 
                mean_var2 = mean(val2), 
                common_mean = mean(c(val1, val2))) %>% 
      ungroup() %>%
      arrange(dist) %>% 
      set_colnames(c(var1, var2, "dist", paste0("mean_", c(var1, var2)), "common_mean"))
  }












