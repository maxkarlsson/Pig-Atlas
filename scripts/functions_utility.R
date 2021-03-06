
savepath <- 
  function(savename, version_ = F) { 
    wd <- getwd()
    result_folder <- paste0(str_extract(wd, ".*Pig Atlas"), "/results/", Sys.Date())
    dir.create(result_folder, showWarnings = FALSE)
    
    suffix <- 
      savename %>% 
      strsplit("\\.") %>%
      {.[[1]] %>%
          .[length(.)]}
    
    savename <-
      savename %>% 
      strsplit("\\.") %>%
      {.[[1]] %>%
          .[-length(.)] %>%
          paste(collapse = "\\.")}
    
    
    
    
    if(version_) {
      vers <- 
        list.files(result_folder) %>% 
        enframe() %>%
        filter(grepl(paste0(savename, "v\\d*.", suffix), value)) %>%
        mutate(str = str_extract(value, paste0("v\\d*\\.", suffix, "$")) %>%
                 gsub(paste0("\\.", suffix, "$"), "", .),
               n = gsub("v", "", str) %>% 
                 as.numeric()) %>%
        arrange(-n) %>% 
        slice(1) %>% 
        pull(n) 
      
      if(length(vers) == 0) {
        vers <- "v1"
      } else {
        vers <- paste0("v", vers + 1)
      }
      
      savename <- 
        paste0(savename, vers)
      
    } 
    
    outp <- paste0(result_folder, "/", savename, ".", suffix)
    
    return(outp)
    
  }

omega_sq <- function(aov_in, neg2zero=T){
  aovtab <- summary(aov_in)[[1]]
  n_terms <- length(aovtab[["Sum Sq"]]) - 1
  output <- rep(-1, n_terms)
  SSr <- aovtab[["Sum Sq"]][n_terms + 1]
  MSr <- aovtab[["Mean Sq"]][n_terms + 1]
  SSt <- sum(aovtab[["Sum Sq"]])
  for(i in 1:n_terms){
    SSm <- aovtab[["Sum Sq"]][i]
    DFm <- aovtab[["Df"]][i]
    output[i] <- (SSm-DFm*MSr)/(SSt+MSr)
    if(neg2zero & output[i] < 0){output[i] <- 0}
  }
  output <- c(output, 1 - sum(output))
  names(output) <- c(rownames(aovtab)[1:n_terms], "Residuals")
  
  return(output)
}

multispread <- function(df, key, value) {
  
  # quote key
  keyq <- rlang::enquo(key)
  # break value vector into quotes
  valueq <- rlang::enquo(value)
  s <- rlang::quos(!!valueq)
  df %>% gather(variable, value, !!!s) %>%
    unite(temp, !!keyq, variable) %>%
    spread(temp, value)
}

list_all_object_sizes <- 
  function() {
    
    
    
    objects(envir = .GlobalEnv) %>%
      enframe() %>%
      mutate(byte = sapply(value,
                           function(x) {
                             eval(parse(text = paste0("pryr::object_size(",
                                                      x,
                                                      ")")))
                           }),
             Mb = byte / 1e6,
             Gb = byte / 1e9)
    
  }


bigcorPar <- function(x, nblocks = 10, verbose = TRUE, ncore="all", ...){
  library(ff, quietly = TRUE)
  require(doMC)
  if(ncore=="all"){
    ncore = multicore:::detectCores()
    registerDoMC(cores = ncore)
  } else{
    registerDoMC(cores = ncore)
  }
  
  NCOL <- ncol(x)
  
  ## test if ncol(x) %% nblocks gives remainder 0
  if (NCOL %% nblocks != 0){stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")}
  
  ## preallocate square matrix of dimension
  ## ncol(x) in 'ff' single format
  corMAT <- ff(vmode = "single", dim = c(NCOL, NCOL))
  
  ## split column numbers into 'nblocks' groups
  SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))
  
  ## create all unique combinations of blocks
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)
  
  ## iterate through each block combination, calculate correlation matrix
  ## between blocks and store them in the preallocated matrix on both
  ## symmetric sides of the diagonal
  results <- foreach(i = 1:nrow(COMBS)) %dopar% {
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")
    flush.console()
    COR <- cor(x[, G1], x[, G2], ...)
    corMAT[G1, G2] <- COR
    corMAT[G2, G1] <- t(COR)
    COR <- NULL
  }
  
  gc()
  return(corMAT)
}


pdf_multiple <- 
  function(plots, filepath, widths, heights) {
    
    filename <- gsub("^.*/", "", filepath)
    filepath_ <- gsub(paste0(filename, "$"), "", filepath)
    
    temp_folder <- paste0(filepath_, "TEMP_FOLDER")
    
    dir.create(temp_folder, showWarnings = FALSE)
    
    for(i in 1:length(plots)) {
      # pdf(paste0(temp_folder, "/TEMP", i, ".pdf"), width = widths[i], height = heights[i])
      ggsave(paste0(temp_folder, "/TEMP", i, ".pdf"), width = widths[i], height = heights[i], plot = plots[[i]], 
             limitsize = FALSE)
      # dev.off()
    }
    
    wd <- getwd()
    
    setwd(temp_folder)
    
    system(paste0('pdftk *pdf cat output "../', filename, '"'), 
           intern = F)
    
    setwd(wd)
    
    unlink(temp_folder, recursive = T)
    
    
  }

