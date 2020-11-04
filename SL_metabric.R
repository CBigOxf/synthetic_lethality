SL_metabric <- function()
  
{
  
  ###########################################################################
  # This function implements the SL genome scan for TCGA_METABRIC dataset
  # It scans through the genome looking for which genes are co-mutated or WT
  # based on tumour samples that contain information for AGO2 and MYC CNA
  ###########################################################################
  
  
  # source the data-set ; the .zip folder contais all data in //inputs subfolder
  # change as appropriate ...
  input_dir <- paste0(getwd(),"//inputs")
  
  # Sink is a command to transfer all console output to a .txt file here
  # This prevents any output while running - can deactivate with command : CloseAllConnections()
  sink("tp_metabric.txt")
  lib_dir <- paste0(getwd(),"/libs")
  .libPaths( c( lib_dir , .libPaths() ) )
  
  print("Loading libraries required...")
  list.of.packages <- c("dplyr","ggplot2","ggrepel","ggpubr","tibble","stringr",
                        "tidyverse", "data.table", 
                        "stringi","openxlsx", "data.table" )
  
  invisible(lapply(list.of.packages, library, character.only = TRUE))
  
  # initialize lists to be populated later on
  master_pvalues <- c()
  ratio1 <- c()
  ratio2 <- c()
  mt <- c() # mutational frequency for each gene
  tp_flag <- c() # twenty percent mutated flag, true or false
  print("Reading METABRIC data-set mutations...")
  mut_matrix_csv <- data.table::fread(file = paste0(input_dir,"/metabric_mutations.txt"))
  
  print("Loading expression data...")
  load(file=paste0(input_dir,"/METABRIC_expression.RData"))
  print("Done...")
  
  print("Reading METABRIC CNA calls...")
  metabric_calls <- as.data.frame(read.table(file = paste0(input_dir,"/metabric_CNA.txt"), sep = '\t', header = TRUE))
  metabric_calls <- metabric_calls[,-2]
  metabric_calls <- metabric_calls %>%  dplyr::filter(Hugo_Symbol %in% c("AGO2","MYC"))
  
  temp <- data.table:: transpose(metabric_calls)
  rownames(temp) <- colnames(metabric_calls)
  colnames(temp) <- c("AGO2","MYC")
  calls <- temp
  colnames(calls)[1] <- "TARGET"
  calls <- add_column(calls, Tumor_Sample_Barcode = colnames(metabric_calls))
  calls <- calls[2:nrow(calls),]
  calls <- calls %>% dplyr::mutate(Tumor_Sample_Barcode = str_replace_all(Tumor_Sample_Barcode, "\\.", "-"))
  calls <-  calls %>%  dplyr::filter(MYC == 0)
  all_genes <- unique(colnames(expr_matrix_csv)[2:ncol(expr_matrix_csv)])
  
  gene_list <- all_genes
  gene_list_mutations <- unique(mut_matrix_csv$Hugo_Symbol)
  
  user <- length(gene_list) # can change to scan first x genes (1:x)
  limit <- 100 # how many pvalues and genes to plot (top in ascending order)
  
  title_cancer = "% done - Analyzing genes... "
  
  #####################################################################################################
  gene <- "AGO2"
  pb <- winProgressBar(title = "progress bar", min = 0,max = user,width = 300)
  
  for (i in 1:user) {
    
    setWinProgressBar(pb,i,title = paste(round(i / user * 100, 0),  title_cancer))
    c_gene <- gene_list[i]
    df <- expr_matrix_csv %>%  dplyr::select(c("Tumor_Sample_Barcode"))
    df <- merge(df,calls,by = "Tumor_Sample_Barcode", all = TRUE)
    target_mutations <- mut_matrix_csv %>%  dplyr::filter(Hugo_Symbol %in% c_gene  &
                                                            !(Variant_Classification %in% "Silent") )
    
    check_1 <- nrow(target_mutations)
    df <- merge(df,target_mutations, by = "Tumor_Sample_Barcode", all = TRUE)
    
    check_2 <- nrow(df)
    check <- check_1/check_2 # mutational frequency for the gene
    
    df <- df[order(df[,'Variant_Classification']), ]
    df[,'Variant_Classification'][is.na(df[,'Variant_Classification'])] <- "WT"
    df <- tibble::add_column(df, status = "MT")
    
    df$status  <- ifelse(df$Variant_Classification == "WT", "WT", "MT")
    df_small <- df
    
    df_small <- tibble::add_column(df_small, AMPL_MUT = 0)
    df_small <- tibble::add_column(df_small, NON_AMPL_MUT = 0)
    df_small <- tibble::add_column(df_small, AMPL_WT = 0)
    df_small <- tibble::add_column(df_small, NON_AMPL_WT = 0)
    
    
    df_small <-  df_small  %>%  dplyr::mutate (TARGET=replace(TARGET, TARGET==-2, "deep_deletion"))   
    df_small <-  df_small  %>%  dplyr::mutate (TARGET=replace(TARGET, TARGET==-1, "shallow_deletion"))   
    df_small <-  df_small  %>%  dplyr::mutate (TARGET=replace(TARGET, TARGET==0, "diploid")) 
    df_small <-  df_small  %>%  dplyr::mutate (TARGET=replace(TARGET, TARGET==1, "gain")) 
    df_small <-  df_small  %>%  dplyr::mutate (TARGET=replace(TARGET, TARGET==2, "amplification")) 
    
    df_small$AMPL_MUT  <- ifelse(df_small$status == "MT" & df_small$TARGET == "amplification", 1, 0)
    df_small$NON_AMPL_MUT  <- ifelse(df_small$status == "MT" & df_small$TARGET != "amplification", 1, 0)
    df_small$AMPL_WT  <- ifelse(df_small$status == "WT" & df_small$TARGET == "amplification", 1, 0)
    df_small$NON_AMPL_WT  <- ifelse(df_small$status == "WT" &  df_small$TARGET != "amplification", 1, 0)
    
    df_small <- df_small[order(df_small[,'status']), ]
    
    
    
    # prepare matrix for Fisher's exact test:
    
    c1 <- table(df_small$AMPL_MUT)  # A
    c2 <- table(df_small$AMPL_WT)  # C
    c3 <- table(df_small$NON_AMPL_MUT) # B
    c4 <- table(df_small$NON_AMPL_WT)  # D
    
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
  write.csv(dfshort,"METABRIC_MT_VS_WT_adjpvalues.csv")
  
  colnames(df)[2] <- "AMPLIFIED_MUTATED/(AMPLIFIED_MUTATED+AMPLIFIED_WT)"
  colnames(df)[3] <- "NON_AMPLIFIED_MUTATED/(NON_AMPLIFIED_MUTATED+NON_AMPLIFIED_WT)"
  write.csv(as.data.frame(df), "genome_results_metabric.csv")
  
  if (user>limit) {df2 <- head(df,limit)}
  png("sp_metabric.png", width=1920, height=1080, res=100)
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
    ggtitle(paste0("Genome-wide comparison of ", gene, " amplification versus MT/WT"))
  
  
  print(sp)
  dev.off()
  return(df)
  
  
}