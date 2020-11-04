SL_CCLE <- function()
  
{
  
  ###########################################################################
  # This function implements the SL genome scan for CCLE_20Q2 dataset
  # It scans through the genome looking for which genes are co-mutated or WT
  # based on tumour samples that contain information for AGO2 and MYC CNA
  ###########################################################################
  
  # source the data-set ; the .zip folder contais all data in //inputs subfolder
  # change as appropriate ...
  input_dir <- paste0(getwd(),"//inputs")
  
  
  print("Loading CCLE data-sets...")
  load(file=paste0(input_dir,"/CCLE_20Q2.RData"))
  print("Done...")
  disease_name <- readline(prompt="Enter type of disease (e.g. breast, lung etc) - 'all' for full CCLE run: ")
  disease_name_bkp <- disease_name
  gene <- readline(prompt="Enter the gene to analyze (e.g. AGO2): ")
  
  if (disease_name == "all") { disease_name = "CCLE.csv"}
  else{disease_name = paste0(disease_name,".csv")}
  
  # Sink is a command to transfer all console output to a .txt file here
  # This prevents any output while running - can deactivate with command : CloseAllConnections()
  sink(paste0("tp_",disease_name_bkp,".txt"))
  
  lib_dir <- paste0(getwd(),"/libs")
  .libPaths( c( lib_dir , .libPaths() ) )
  
  print("Loading libraries required...")
  list.of.packages <- c("dplyr","ggplot2","ggrepel","ggpubr","tibble","stringr",
                        "tidyverse", "data.table", 
                        "stringi","openxlsx", "data.table" )
  invisible(lapply(list.of.packages, library, character.only = TRUE))
  colnames(expr_matrix_csv)[1] <- "CELLLINE"
  
  
  if (disease_name == "CCLE.csv") {
    print("Running through all cancer cell lines in CCLE...")
  }
  else{
    print(paste0("Running for ",disease_name))
    cell_line_dir =  paste0(input_dir,"/CELL_LINES")
    command <- paste0(cell_line_dir,"//",disease_name)
    command <- shQuote(command)
    command <- paste0("read.csv(",command)
    command <- paste0(command,",header = TRUE)")
    command <-paste0("as.data.frame(",command)
    command <- paste0(command,")")
    disease_csv <- eval(parse(text = command))
    disease_cell_lines = disease_csv[, 1]
    expr_matrix_csv <- expr_matrix_csv %>%  dplyr::filter(CELLLINE %in% disease_cell_lines)
    mut_matrix_csv <- mut_matrix_csv %>% dplyr::filter(Tumor_Sample_Barcode %in% disease_cell_lines)
    
  }
  # initialize lists to be populated later on
  master_pvalues <- c()
  ratio1 <- c()
  ratio2 <- c()
  mt <- c() # mutational frequency for each gene
  tp_flag <- c() # twenty percent mutated flag, true or false
  all_genes <- colnames(expr_matrix_csv)
  #all_genes <- unique(mut_matrix_csv$Hugo_Symbol)
  
  gene_list <-   sapply(strsplit(all_genes[2:length(all_genes)], split='..', fixed=TRUE),function(x) (x[1]))
  gene_list_IDs <- sapply(strsplit(all_genes[2:length(all_genes)], split='..', fixed=TRUE),function(x) (x[2]))
  gene_list_IDs <- sapply(str_sub(gene_list_IDs,1,-2),function(x) (x[1]))
  user <- length(gene_list) # can change to scan first x genes (1:x)
  limit <- 100 # how many pvalues and genes to plot (top in ascending order)
  
  title_cancer = paste0("% done - Analyzing genes in  ",disease_name_bkp," ...")
  
  #####################################################################################################
  
  cell_line_names <- as.data.frame(read.csv(paste0(input_dir,"/cell_line_names.csv"), header = TRUE))
  names_mat <-  cell_line_names  %>% dplyr::select(c("BROAD_ID","CCLE_ID"))
  colnames(names_mat) <- c("Tumor_Sample_Barcode","NAME")
  
  calls <- read.table(file = paste0(input_dir,"/calls.txt"), sep = '\t', header = TRUE)
  colnames(calls)[3] <- "TARGET"
  colnames(calls)[4] <- "MYC"
  calls <- calls[,c(1,2,3,4)]
  colnames(expr_matrix_csv)[1] <- "Tumor_Sample_Barcode"
  colnames(calls)[2] <- "NAME"
  pb <- winProgressBar(title = "progress bar", min = 0,max = user,width = 300)
  
  df_t <- expr_matrix_csv %>%  dplyr::select(c("Tumor_Sample_Barcode"))
  df_t <- merge(df_t, names_mat, by="Tumor_Sample_Barcode")
  df_t <- merge(df_t,calls, by  = "NAME")
  
  # filter by MYC diploidy:
  
  df_t <- df_t %>% dplyr::filter(MYC == 0)
  colnames(df_t)[2] <- "Tumor_Sample_Barcode"
  
  for (i in 1:user) {
    setWinProgressBar(pb,i,title = paste(round(i / user * 100, 0),  title_cancer))
    c_gene <- gene_list[i]
    target_mutations <- mut_matrix_csv %>%  dplyr::filter(Hugo_Symbol %in% c_gene &
                                                            !(Variant_Classification %in% "Silent") )
    check_1 <- nrow(target_mutations)
    df <- merge(df_t,target_mutations, by = "Tumor_Sample_Barcode", all = TRUE)
    check_2 <- nrow(df)
    check <- check_1/check_2 # mutational frequency for the gene
    
    df <- df[order(df[,'Variant_Classification']), ]
    df[,'Variant_Classification'][is.na(df[,'Variant_Classification'])] <- "WT"
    df[,'isDeleterious'][is.na(df[,'isDeleterious'])] <- "FALSE"
    df <- tibble::add_column(df, status = "MT")
    df$status  <- ifelse(df$Variant_Classification == "WT", "WT", "MT")
    df_small <- df %>%  dplyr::select(c("Tumor_Sample_Barcode",
                                        "Variant_Classification", "TARGET","status","isDeleterious"))
    
    
    # -2 or Deep Deletion indicates a deep loss, possibly a homozygous deletion
    # -1 or Shallow Deletion indicates a shallow loss, possibley a heterozygous deletion
    # 0 is diploid
    # 1 or Gain indicates a low-level gain (a few additional copies, often broad)
    # 2 or Amplification indicate a high-level amplification (more copies, often focal)
    
    df_small <-  df_small  %>%  dplyr::mutate (TARGET=replace(TARGET, TARGET==-2, "deep_deletion"))   
    df_small <-  df_small  %>%  dplyr::mutate (TARGET=replace(TARGET, TARGET==-1, "shallow_deletion"))   
    df_small <-  df_small  %>%  dplyr::mutate (TARGET=replace(TARGET, TARGET==0, "diploid")) 
    df_small <-  df_small  %>%  dplyr::mutate (TARGET=replace(TARGET, TARGET==1, "gain")) 
    df_small <-  df_small  %>%  dplyr::mutate (TARGET=replace(TARGET, TARGET==2, "amplification")) 
    
    df_small <- df_small[order(df_small[,"TARGET"]), ]
    
    ##################################
    
    df_small <- tibble::add_column(df_small, AMPL_MUT = 0)
    df_small <- tibble::add_column(df_small, NON_AMPL_MUT = 0)
    df_small <- tibble::add_column(df_small, AMPL_WT = 0)
    df_small <- tibble::add_column(df_small, NON_AMPL_WT = 0)
    
    df_small$AMPL_MUT  <- ifelse(df_small$status == "MT" & df_small$TARGET == "amplification", 1, 0)
    df_small$NON_AMPL_MUT  <- ifelse(df_small$status == "MT" & df_small$TARGET != "amplification", 1, 0)
    df_small$AMPL_WT  <- ifelse(df_small$status == "WT" & df_small$TARGET == "amplification", 1, 0)
    df_small$NON_AMPL_WT  <- ifelse(df_small$status == "WT" & df_small$TARGET != "amplification", 1, 0)
    
    df_small <- na.omit(df_small)
    df_small <- df_small[order(df_small[,'status']), ]
    
    # prepare matrix for Fisher's exact test:
    
    c1 <- table(df_small$AMPL_MUT)
    c2 <- table(df_small$AMPL_WT)
    c3 <- table(df_small$NON_AMPL_MUT)
    c4 <- table(df_small$NON_AMPL_WT)
    
    c <- matrix(c(c1[2],c2[2],c3[2],c4[2]),2,2)
    
    colnames(c) <- c(paste0(gene," amplified"),paste0(gene," not-amplified (deep_deletion/diploid)"))
    rownames(c) <- c(paste0(c_gene," mutated"),paste0(c_gene," WT"))
    c[is.na(c)] <- 0
    print("'######################################################################")
    print(c)
    res <- fisher.test(c)
    
    
    # now filter out anything that is mutant in less than 20% of samples:
    ratio1 <- c(ratio1, c[1,1]/(c[1,1]+c[2,1]))
    ratio2 <- c(ratio2,c[1,2]/(c[1,2]+c[2,2]))
    
    
    if ( check >= 0.2 ) {
      
      tp_flag <- c(tp_flag,"TRUE")
    }
    else{
      tp_flag <- c(tp_flag,"FALSE")
    }
    
    
    master_pvalues <- c(master_pvalues, res$p.value)
    mt <- c(mt, check)
    
    print(res)
    
    
  }
  close(pb)
  names(master_pvalues) <- gene_list[1:user]
  
  df <- matrix(data = , nrow = user, ncol = 1) 
  df[, 1] = master_pvalues
  row.names(df) <- names(master_pvalues)
  colnames(df) <- "pvalues"
  
  df <- data.frame(df)
  
  # add ratio columns
  df <- add_column(df, ratios1 = 0)
  df <- add_column(df, ratios2 = 0)
  df <- add_column(df,tp_flags = "FALSE")
  df <- add_column(df,Mutational_Freq = 0)
  df$ratios1 <- ratio1
  df$ratios2 <- ratio2
  df$tp_flags <- tp_flag
  df$Mutational_Freq <- mt
  
  df <- df %>% arrange(desc(-pvalues))
  # adjust p values
  dfshort <- df %>% dplyr::filter(tp_flags=="TRUE") # & pvalues < 0.05)
  dfshort <- add_column(dfshort, adjusted_pvalue = 0)
  dfshort$adjusted_pvalue <- p.adjust(dfshort$pvalues, method = "BH")
  write.csv(dfshort,paste0("CCLE_MT_VS_WT_adjpvalues_",disease_name))
  
  colnames(df)[2] <- "AMPLIFIED_MUTATED/(AMPLIFIED_MUTATED+AMPLIFIED_WT)"
  colnames(df)[3] <- "NON_AMPLIFIED_MUTATED/(NON_AMPLIFIED_MUTATED+NON_AMPLIFIED_WT)"
  write.csv(as.data.frame(df), paste0("genome_results_",disease_name))
  
  if (user>limit) {df2 <- head(df,limit)}
  png(paste0("sp_",disease_name_bkp,".png"), width=1920, height=1080, res=100)
  sp <- ggplot(data = df2,aes(x=factor(row.names(df2), levels=row.names(df2)), y = pvalues ))+
    geom_point(size = 4,color = dplyr::case_when(df2$pvalues > 0.05 ~ "#FF0000",
                                                 df2$pvalues < 0.05 ~ "#00CC00",
                                                 TRUE ~ "#00CC00"), alpha = 0.8) +
    geom_hline(yintercept = 0.05, color = "#00CC00") +
    
    geom_label_repel(aes(label=row.names(df2)),
                     box.padding   = 0.5,
                     point.padding = 0.005, size = 2) +
    
    
    ylab(  c("P-value of Fisher exact test")  )  +
    xlab(  c("Hugo Symbol") ) +
    
    font("xylab",size=10)+
    font("xy",size=10)+
    font("xy.text", size = 10) +
    font("legend.text",size = 10) +
    theme(axis.text.x=element_text(size=10, angle=90,hjust=0.95,vjust=0.02)) +
    ggtitle(paste0("Genome-wide comparison of ", gene, " amplification VS MT/WT"))
  
  
  print(sp)
  dev.off()
  return(df)
  
  
}