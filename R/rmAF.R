#' Remove High Allele Frequency
#'
#' This function does two things: 1) Separate a single position (nucleotide) with a maximum three different allele frequency reported from CRISP .vcf file into three different rows with the corresponding altered result and frequency, and 2) remove SNVs that have equal or greater than the designated AF.
#' @param 1) Input the raw .vcf files in a .txt format, 2) Name the output file, and 3) Designate the AF that you wanted to remove (from 0.0 to 1.0).
#' @keywords Allele Frequency
#' @export
#' @examples
#' rmAF()

rmAF <- function(input, output, af){
  library(dplyr)
  a <- read.table(input, header=F, stringsAsFactor = FALSE)
  res <- str_match(a$V8, "AF=\\s*(.*?)\\s*;EM") %>% as.data.frame() %>% select(V6 =V2)
  b<- a %>% select(V1,V2,V3,V4,V5) %>% cbind(.,res)
  temp_df <- data.frame(CHROM = as.character(),
                        POS = as.character(),
                        ID = as.character(),
                        REF = as.character(),
                        ALT = as.character(),
                        QUAL = as.character(),
                        FILTER = as.character(),
                        INFO = as.character(),
                        FREQ = as.numeric())
  for(i in 1:nrow(b)) {
    if(str_count(b[i,6],",") == 0){
      temp <- data.frame(CHROM = b[i,1],
                         POS = b[i,2],
                         ID = b[i,3],
                         REF = b[i,4],
                         ALT = b[i,5],
                         QUAL = ".",
                         FILTER = ".",
                         INFO = ".",
                         FREQ = b[i,6])
    } else if (str_count(b[i,6],",") == 1){
      temp <- data.frame(CHROM = c(b[i,1], b[i,1]),
                         POS = c(b[i,2], b[i,2]),
                         ID = c(b[i,3], b[i,3]),
                         REF = c(b[i,4], b[i,4]),
                         ALT = c(str_match(b[i,5],"\\s*(.*?)\\s*,")[,2],
                                 str_match(b[i,5],"^.+,(.*)")[,2]),
                         QUAL = c(".","."),
                         FILTER = c(".","."),
                         INFO = c(".","."),
                         FREQ = c(str_match(b[62,6],"\\s*(.*?)\\s*,")[,2],
                                  str_match(b[62,6],"^.+,(.*)")[,2]))
    } else if (str_count(b[i,6],",") == 2){
      temp <- data.frame(CHROM = c(b[i,1], b[i,1], b[i,1]),
                         POS = c(b[i,2], b[i,2], b[i,2]),
                         ID = c(b[i,3], b[i,3], b[i,3]),
                         REF = c(b[i,4], b[i,4], b[i,4]),
                         ALT = c(str_match(b[i,5],"\\s*(.*?)\\s*,")[,2],
                                 str_match(b[i,5],"^.+,(.*),")[,2],
                                 str_match(b[i,5],"^.+,(.*)")[,2]),
                         QUAL = c(".",".","."),
                         FILTER = c(".",".","."),
                         INFO = c(".",".","."),
                         FREQ = c(str_match(b[i,6],"\\s*(.*?)\\s*,")[,2],
                                  str_match(b[i,6],"^.+,(.*),")[,2],
                                  str_match(b[i,6],"^.+,(.*)")[,2]))
    } else if (str_count(b[i,6],",") == 3){
      temp <- data.frame(CHROM = c(b[i,1], b[i,1], b[i,1], b[i,1]),
                         POS = c(b[i,2], b[i,2], b[i,2], b[i,2]),
                         ID = c(b[i,3], b[i,3], b[i,3], b[i,3]),
                         REF = c(b[i,4], b[i,4], b[i,4], b[i,4]),
                         ALT = c(str_match(b[i,5],"\\s*(.*?)\\s*,")[,2],
                                 str_match(b[i,5],"^.+,(.*),.+,")[,2],
                                 str_match(b[i,5],"^.+,(.*),.+")[,2],
                                 str_match(b[i,5],"^.+,(.*)")[,2]),
                         QUAL = c(".",".",".","."),
                         FILTER = c(".",".",".","."),
                         INFO = c(".",".",".","."),
                         FREQ = c(str_match(b[i,6],"\\s*(.*?)\\s*,")[,2],
                                  str_match(b[i,6],"^.+,(.*),.+,")[,2],
                                  str_match(b[i,6],"^.+,(.*),.+")[,2],
                                  str_match(b[i,6],"^.+,(.*)")[,2]))
    } else {
      print("more than 4 different allele frequencies exist")
    }
    temp_df <- rbind(temp_df, temp)
  }
  temp_df2<- temp_df %>% filter(FREQ < af) %>% select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)
  write.table(temp_df2, output, quote = FALSE, sep="\t", row.names = FALSE)
  rm(a,b,res,temp,temp_df)
}
