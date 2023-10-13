#' Analyze pooled sequencing data
#'
#' This function generate a table including mutations density, deleterious ratio, and KS p-value of every coding gene.
#' @param 1) Input the original VEP .txt file generated from Ensembl, 2) Input the SIFT .csv file generayed from SIFT 4G, 3) Remove mutations in duplicated regions (T/F), and 4) Name the output file.
#' @keywords mutations density, deleterious ratio, KS p-value
#' @export
#' @examples
#' pooled_analysis()

pooled_analysis <- function(raw_VEP, SIFT_file, Dup, output){
  library(dplyr)
  library(tidyr)
  data(hardmask_list)
  # Remove mutations in low complicity region (duplicated region)
  if(Dup == TRUE){
  list <- hardmask_list %>% unite("merge", c(CHR, Seq), remove = F, sep = ":") %>%
    unite("Location", c(merge, Seq), remove = T, sep = "-") %>%
    select(Location)
  df_rm <- anti_join(raw_VEP, list, by = "Location")
  df2_list <- SIFT_file %>% unite("merge", c(CHROM, POS), remove = F, sep = ":") %>%
    unite("Location", c(merge, POS), remove = F, sep = "-")
  df2_rm <- anti_join(df2_list, list, by = "Location") %>% select(-Location, -merge)
  a <- df_rm %>% select(2,4,5,7) %>% group_by(Location) %>% unique() %>% ungroup()
  working_df <- data.frame(WBGENEID=a$Gene, IMPACT=a$IMPACT, Consequence=a$Consequence) %>% filter(!IMPACT=="MODIFIER")
  rm(a,list,df2_list)
  # Count total mutations of each gene (includes up/down stream and science mutation)
  a <- data.frame(table(working_df$WBGENEID))
  totalFreq_df <- data.frame(WBGENEID = a$Var1, TotalFreq = a$Freq)
  rm(a)
  # Count Non-Synonymous mutations in each gene
  non_syn_df <- filter(working_df, IMPACT == "HIGH"| IMPACT == "MODERATE") %>% select(1,2,3)
  a <- data.frame(table(non_syn_df$WBGENEID))
  NonSynMut_df <- data.frame(WBGENEID = a$Var1, Freq = a$Freq)
  rm(a)
  a <- subset(totalFreq_df, WBGENEID %in% NonSynMut_df$WBGENEID)
  b <- a$TotalFreq - NonSynMut_df$Freq
  c <- data.frame(WBGENEID = a$WBGENEID, SynMut = b, NonSynMut = NonSynMut_df$Freq)
  d <- subset(totalFreq_df, !(WBGENEID %in% NonSynMut_df$WBGENEID))
  e <- data.frame(WBGENEID = d$WBGENEID, SynMut = d$TotalFreq, NonSynMut = 0)
  temp1 <- rbind(c,e)
  rm(a,b,c,d,e)
  # High list
  a <- filter(working_df, IMPACT == "HIGH") %>% select(1) %>% unique()
  b <- data.frame(WBGENEID = a$WBGENEID, Impact = "High")
  # Moderate
  c <- filter(working_df, IMPACT == "MODERATE") %>% select(1) %>% unique()
  d <- subset(totalFreq_df, WBGENEID %in% c$WBGENEID)
  e <- subset(d, !(WBGENEID %in% a$WBGENEID))
  f <- data.frame(WBGENEID = e$WBGENEID, Impact = "Moderate")
  # Low
  g <- subset(totalFreq_df, !(WBGENEID %in% NonSynMut_df$WBGENEID))
  h <- data.frame(WBGENEID = g$WBGENEID, Impact = "Low")
  i <- rbind(b,f)
  j <- rbind(i,h)
  Semi_df <- merge(temp1, j, by ="WBGENEID")
  rm(a,b,c,d,e,f,g,h,i,j,temp1)
  # Calculate Nonsense mutation
  a <- filter(working_df, Consequence == "stop_gained") %>%
    select(1,2,3)
  b <- data.frame(table(a$WBGENEID))
  Nonsense_df <- data.frame(WBGENEID = b$Var1, NonsenseMut = b$Freq)
  c <- setdiff(totalFreq_df$WBGENEID, Nonsense_df$WBGENEID)
  d <- data.frame(WBGENEID = c,  NonsenseMut= 0)
  e <- rbind(Nonsense_df, d)
  semisemi_df <- merge(Semi_df, e, by = "WBGENEID")
  rm(a,b,c,d,e,Nonsense_df,Semi_df)
  # Add gene information such as position and expression pattern
  data(geneinfo)
  final_df <- merge(semisemi_df, geneinfo, by = "WBGENEID") %>%
    mutate(NonSynMut_perKb=NonSynMut/(CDS_Length/1000)) %>%
    mutate_if(is.factor, as.character)
  rm(geneinfo,non_syn_df, NonSynMut_df, totalFreq_df, working_df)
  # KS test
  a <- raw_VEP %>% filter(Consequence == "missense_variant" |
                            Consequence == "missense_variant,splice_region_variant"|
                            Consequence == "splice_acceptor_variant,missense_variant") %>%
    select(c("Location","Gene", "CDS_position", "Feature")) %>%
    group_by(Location) %>% unique() %>% dplyr::slice(which.max(CDS_position)) %>% ungroup() %>%
    group_by(Gene) %>% mutate(duplicate = n()) %>% filter(duplicate > 1) %>% select(-duplicate) %>% ungroup()
  c <- a %>% select("Gene") %>% unique() %>% as.data.frame()
  x <- final_df %>% filter(NonSynMut > 1) %>% select(WBGENEID, CDS_Length) %>% unique()
  temp_df <- data.frame(WBGENEID = character(),
                        p_kstest = numeric(),
                        stringsAsFactors = FALSE)
  for(i in 1:nrow(c)){
    step1 <- x %>% filter(WBGENEID == c[[1]][[i]]) %>% select("CDS_Length") %>% unique() %>% as.numeric()
    step2 <- a %>% filter (Gene == c[[1]][[i]]) %>% select("CDS_position") %>% mutate_all(as.numeric)
    step3 <- ks.test(step2, "punif", 1, step1)$p.value
    step4 <- data.frame(WBGENEID = c[[1]][[i]], p_kstest = step3)
    temp_df <- rbind(temp_df, step4)
  }
  rm(a,c,x,step1,step2,step3,step4)
  ks_df <- merge(final_df,temp_df,by="WBGENEID",all = TRUE)
  rm(temp_df,final_df)
  # Deleterious assay
  a <- df2_rm %>% na.omit() %>% select(POS,AMINO_POS,GENE_ID, SIFT_PREDICTION) %>%
    group_by(POS) %>% unique() %>% dplyr::slice(which.max(AMINO_POS)) %>% ungroup() %>%
    select(GENE_ID, SIFT_PREDICTION)
  c <- a$GENE_ID[duplicated(a$GENE_ID)] %>% unique() %>% as.data.frame()
  temp_df <- data.frame(WBGENEID = character(),
                        DeleteriousRatio = numeric(),
                        Total = numeric(),
                        stringsAsFactors = FALSE)
  for(i in 1:nrow(c)){
    step1 <- a %>% filter (GENE_ID == c[[1]][[i]]) %>% nrow() %>% as.numeric()
    step2 <- a %>% filter (GENE_ID == c[[1]][[i]]) %>%
      filter(SIFT_PREDICTION == "DELETERIOUS (*WARNING! Low confidence)" | SIFT_PREDICTION == "DELETERIOUS") %>% nrow() %>% as.numeric()
    step3 <- step2/step1
    step4 <- data.frame(WBGENEID = c[[1]][[i]],
                        DeleteriousRatio = step3,
                        Total = step1)
    temp_df <- rbind(temp_df, step4)
  }
  SIFT_df <- merge(ks_df,temp_df,by="WBGENEID",all = TRUE) %>% select(
    WBGENEID,
    Gene_Name,
    Gene_Length,
    CDS_ID,
    CDS_Length,
    CHR,
    POS,
    SeqStart,
    Impact,
    SynMut,
    NonSynMut,
    NonsenseMut,
    NonSynMut_perKb,
    DeleteriousRatio,
    p_kstest,
    PharynxEx,
    IntesEx,
    GermlineEx)
  write.table(SIFT_df, output, quote = FALSE, sep="\t", row.names = FALSE)
  rm(a,c,step1,step2,step3,step4,temp_df,SIFT_df,ks_df,df_rm,df2_rm)
  } else
    {
    a <- raw_VEP %>% select(2,4,5,7) %>% group_by(Location) %>% unique() %>% ungroup()
    working_df <- data.frame(WBGENEID=a$Gene, IMPACT=a$IMPACT, Consequence=a$Consequence) %>% filter(!IMPACT=="MODIFIER")
    rm(a)
    # Count total mutations of each gene (includes up/down stream and science mutation)
    a <- data.frame(table(working_df$WBGENEID))
    totalFreq_df <- data.frame(WBGENEID = a$Var1, TotalFreq = a$Freq)
    rm(a)
    # Count Non-Synonymous mutations in each gene
    non_syn_df <- filter(working_df, IMPACT == "HIGH"| IMPACT == "MODERATE") %>% select(1,2,3)
    a <- data.frame(table(non_syn_df$WBGENEID))
    NonSynMut_df <- data.frame(WBGENEID = a$Var1, Freq = a$Freq)
    rm(a)
    a <- subset(totalFreq_df, WBGENEID %in% NonSynMut_df$WBGENEID)
    b <- a$TotalFreq - NonSynMut_df$Freq
    c <- data.frame(WBGENEID = a$WBGENEID, SynMut = b, NonSynMut = NonSynMut_df$Freq)
    d <- subset(totalFreq_df, !(WBGENEID %in% NonSynMut_df$WBGENEID))
    e <- data.frame(WBGENEID = d$WBGENEID, SynMut = d$TotalFreq, NonSynMut = 0)
    temp1 <- rbind(c,e)
    rm(a,b,c,d,e)
    # High list
    a <- filter(working_df, IMPACT == "HIGH") %>% select(1) %>% unique()
    b <- data.frame(WBGENEID = a$WBGENEID, Impact = "High")
    # Moderate
    c <- filter(working_df, IMPACT == "MODERATE") %>% select(1) %>% unique()
    d <- subset(totalFreq_df, WBGENEID %in% c$WBGENEID)
    e <- subset(d, !(WBGENEID %in% a$WBGENEID))
    f <- data.frame(WBGENEID = e$WBGENEID, Impact = "Moderate")
    # Low
    g <- subset(totalFreq_df, !(WBGENEID %in% NonSynMut_df$WBGENEID))
    h <- data.frame(WBGENEID = g$WBGENEID, Impact = "Low")
    i <- rbind(b,f)
    j <- rbind(i,h)
    Semi_df <- merge(temp1, j, by ="WBGENEID")
    rm(a,b,c,d,e,f,g,h,i,j,temp1)
    # Calculate Nonsense mutation
    a <- filter(working_df, Consequence == "stop_gained") %>%
      select(1,2,3)
    b <- data.frame(table(a$WBGENEID))
    Nonsense_df <- data.frame(WBGENEID = b$Var1, NonsenseMut = b$Freq)
    c <- setdiff(totalFreq_df$WBGENEID, Nonsense_df$WBGENEID)
    d <- data.frame(WBGENEID = c,  NonsenseMut= 0)
    e <- rbind(Nonsense_df, d)
    semisemi_df <- merge(Semi_df, e, by = "WBGENEID")
    rm(a,b,c,d,e,Nonsense_df,Semi_df)
    # Add gene information such as position and expression pattern
    data(geneinfo)
    final_df <- merge(semisemi_df, geneinfo, by = "WBGENEID") %>%
      mutate(NonSynMut_perKb=NonSynMut/(CDS_Length/1000)) %>%
      mutate_if(is.factor, as.character)
    rm(non_syn_df, NonSynMut_df, totalFreq_df, working_df)
    # KS test
    a <- raw_VEP %>% filter(Consequence == "missense_variant" |
                              Consequence == "missense_variant,splice_region_variant"|
                              Consequence == "splice_acceptor_variant,missense_variant") %>%
      select(c("Location","Gene", "CDS_position", "Feature")) %>%
      group_by(Location) %>% unique() %>% dplyr::slice(which.max(CDS_position)) %>% ungroup() %>%
      group_by(Gene) %>% mutate(duplicate = n()) %>% filter(duplicate > 1) %>% select(-duplicate) %>% ungroup()
    c <- a %>% select("Gene") %>% unique() %>% as.data.frame()
    x <- final_df %>% filter(NonSynMut > 1) %>% select(WBGENEID, CDS_Length) %>% unique()
    temp_df <- data.frame(WBGENEID = character(),
                          p_kstest = numeric(),
                          stringsAsFactors = FALSE)
    for(i in 1:nrow(c)){
      step1 <- x %>% filter(WBGENEID == c[[1]][[i]]) %>% select("CDS_Length") %>% unique() %>% as.numeric()
      step2 <- a %>% filter (Gene == c[[1]][[i]]) %>% select("CDS_position") %>% mutate_all(as.numeric)
      step3 <- ks.test(step2, "punif", 1, step1)$p.value
      step4 <- data.frame(WBGENEID = c[[1]][[i]], p_kstest = step3)
      temp_df <- rbind(temp_df, step4)
    }
    rm(a,c,x,step1,step2,step3,step4)
    ks_df <- merge(final_df,temp_df,by="WBGENEID",all = TRUE)
    rm(temp_df,final_df)
    # Deleterious assay
    a <- SIFT_file %>% na.omit() %>% select(POS,AMINO_POS,GENE_ID, SIFT_PREDICTION) %>%
      group_by(POS) %>% unique() %>% dplyr::slice(which.max(AMINO_POS)) %>% ungroup() %>%
      select(GENE_ID, SIFT_PREDICTION)
    c <- a$GENE_ID[duplicated(a$GENE_ID)] %>% unique() %>% as.data.frame()
    temp_df <- data.frame(WBGENEID = character(),
                          DeleteriousRatio = numeric(),
                          Total = numeric(),
                          stringsAsFactors = FALSE)
    for(i in 1:nrow(c)){
      step1 <- a %>% filter (GENE_ID == c[[1]][[i]]) %>% nrow() %>% as.numeric()
      step2 <- a %>% filter (GENE_ID == c[[1]][[i]]) %>%
        filter(SIFT_PREDICTION == "DELETERIOUS (*WARNING! Low confidence)" | SIFT_PREDICTION == "DELETERIOUS") %>% nrow() %>% as.numeric()
      step3 <- step2/step1
      step4 <- data.frame(WBGENEID = c[[1]][[i]],
                          DeleteriousRatio = step3,
                          Total = step1)
      temp_df <- rbind(temp_df, step4)
    }
    SIFT_df <- merge(ks_df,temp_df,by="WBGENEID",all = TRUE) %>% select(
      WBGENEID,
      Gene_Name,
      Gene_Length,
      CDS_ID,
      CDS_Length,
      CHR,
      POS,
      SeqStart,
      Impact,
      SynMut,
      NonSynMut,
      NonsenseMut,
      NonSynMut_perKb,
      DeleteriousRatio,
      p_kstest,
      PharynxEx,
      IntesEx,
      GermlineEx)
    write.table(SIFT_df, output, quote = FALSE, sep="\t", row.names = FALSE)
    rm(a,c,step1,step2,step3,step4,temp_df,SIFT_df,ks_df)
  }
  # Reorganize data as CHR, POS, ID, Syn, NonSyn, Nonsense, HighImpact, TotalFreq
}
