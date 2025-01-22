cat("\n ###### DNPcall ######")
cat("\n")
cat("\n")
cat("\n")
cat("\n")
# package needed
suppressPackageStartupMessages(library('tidyverse'))
suppressPackageStartupMessages(library('stringr'))

# read the arguments given
args <- commandArgs(trailingOnly = TRUE)

error0 = c()
if(length(grep("bamlist",args)) != 0) {
  bamlist = strsplit(args[grep("--bamlist*",args)],"=")[[1]][[2]]
} else {
  error0 = c(error0,"bamlist")
}
if(length(grep("DNPs",args)) != 0) {
  DNPlist = strsplit(args[grep("--DNPs*",args)],"=")[[1]][[2]]
} else {
  error0 = c(error0,"DNPs")
}
if(length(grep("reference",args)) != 0) {
  reference = strsplit(args[grep("--reference*",args)],"=")[[1]][[2]]
} else {
  error0 = c(error0,"reference")
}
if(length(error0)>0) {
  cat("\n\t # ERROR: These required arguments could not be found:\n")
  cat(paste0("\t\t"),paste0(error0,collapse = ", "))
  cat("\n")
  quit()
}

# read the df
r_bamlist<-read.table(bamlist,header=F,col.names = "samples")
r_DNPlist<-read.table(DNPlist,header=F,col.names = c("chr","pos1","pos2","REF","ALT"))

cat("Samples analyzed:\n")
print(r_bamlist)
cat("\n")

# make .bed file for calling
bedfile <- as.data.frame(cbind(r_DNPlist$chr, r_DNPlist$pos1-2, r_DNPlist$pos2+1))
write.table(bedfile,file = "bedfile_for_mpileup.bed",row.names = F,col.names = F,quote=F,sep="\t")

# make a df to assign a unique label for each DNP, necessary if there are very close DNPs falling in the same read
ForUniqReads1 <- as.data.frame(cbind(paste0(r_DNPlist$chr,":",r_DNPlist$pos1-1),paste0(r_DNPlist$chr,":",r_DNPlist$pos1),
                                     paste0(r_DNPlist$chr,":",r_DNPlist$pos2), paste0(r_DNPlist$chr,":",r_DNPlist$pos2+1),
                                     seq.int(nrow(r_DNPlist))))
df_int <- gather(data=ForUniqReads1,key="col", value="pos", 1:4)
ForUniqReads <- df_int[,c(3,1)]
colnames(ForUniqReads) <- c("df.chr","uniq")

# make a DNP-based df
RefAltDNP <- as.data.frame(cbind(paste0(r_DNPlist$chr,":",r_DNPlist$pos1,"|",r_DNPlist$chr,":",r_DNPlist$pos2),
                                 gsub('^(.{1})(.*)$', '\\1|\\2',r_DNPlist$REF),
                                 gsub('^(.{1})(.*)$', '\\1|\\2',r_DNPlist$ALT)))
colnames(RefAltDNP)<-c("POS","ref","alt")

# pileup calling with samtools
cat("\nMaking mpileup...\n")
commA = paste0("samtools mpileup -b ",bamlist,
               " -l bedfile_for_mpileup.bed -f ",reference,
               " --output-QNAME --no-output-ins --no-output-del --no-output-ends > mpileup.txt")

system(commA, ignore.stderr=T)

# read sample list, remove path and ".bam" extension in the name
list0<-scan(bamlist,character(),quiet=T)
list <- gsub("(.*/\\s*(.*$))", "\\2", list0)
samplelist<-gsub(".bam", "", list)

# reading pileup file in r
cat("\nReading mpileup...\n")
pileup_bamlist <- read.csv("mpileup.txt", quote="", row.names=NULL, stringsAsFactors=F, header =F, sep="\t")
common_columns <- 3
final_merged_data <- NULL

# calling processing for each sample
for (i in samplelist) {

  cat(paste0("\nProcessing sample: ", i, "\n"))

  # extract sample-specific columns
  depth_col <- common_columns + (which(samplelist == i) - 1) * 4 + 1
  reads_col <- depth_col + 1
  qual_col <- reads_col + 1
  id_col <- qual_col + 1

  indiv_df <- pileup_bamlist %>%
    select(V1 = V1, V2 = V2, V3 = V3,
           V4 = !!sym(paste0("V", depth_col)),
           V5 = !!sym(paste0("V", reads_col)),
           V6 = !!sym(paste0("V", qual_col)),
           V7 = !!sym(paste0("V", id_col)))

  # substitute . and , with reference base
  indiv_df$V5 <- mapply(function(v3, v5) gsub("[.,]", v3, v5), indiv_df$V3, indiv_df$V5)

  # create "chr:position" column
  indiv_df$chr <- paste(indiv_df$V1, indiv_df$V2, sep = ":")

  # select necessary columns
  df2 <- data.frame(df.chr = indiv_df$chr, df.V5 = indiv_df$V5, df.V7 = indiv_df$V7)

  # substitution of characters to modify the df
  df2$df.V7 <- paste0(",", df2$df.V7, ",")
  df2$df.V5 <- gsub('\\*', '@', df2$df.V5)
  df2$df.V5 <- gsub('[actgACTG]([+-]\\d+)+', '@', df2$df.V5)
  df2$df.V5 <- gsub('', ' ', df2$df.V5)
  df2$df.V7 <- gsub(',', ' ', df2$df.V7)

  # split bases and reads columns in single cells
  df3 <- df2 %>% separate_rows(df.V5, df.V7, sep = " ")
  df4 <- df3[,c(3,1,2)]
  df4.1 <-  df4[!(df4$df.V7==""), ]

  # get unique label for each DNP
  df5 <- left_join(df4.1,ForUniqReads,by="df.chr")
  df5$readUniq <- paste(df5$df.V7,df5$uniq,sep="_")
  df6 <- df5[,c(5,2,3)]

  # group by new read label
  DF_chr <- df6 %>%
    group_by(readUniq) %>%
    summarise(BASE = paste(df.chr, collapse = "|"))
  DF_pos <- df6 %>%
    group_by(readUniq) %>%
    summarise(BASE = paste(df.V5, collapse = "|"))

  DF <- merge(DF_chr,DF_pos, by='readUniq')
  names(DF) <- c("ID_DNP","chr","pos")

  # substitute characters strand specific characters (unecessary) with only uppercase letters
  DF$pos = gsub('g', 'G', DF$pos)
  DF$pos = gsub('a', 'A', DF$pos)
  DF$pos = gsub('c', 'C', DF$pos)
  DF$pos = gsub('t', 'T', DF$pos)
  DF$pos = gsub(',', '.', DF$pos)

  # eliminate reads with "@", indicating reads with INDELs
  DFnoINDEL <- DF %>%
    filter_all(all_vars(!grepl("@", .)))

  # keep only micro haplotype
  righe_modificate <- DFnoINDEL[nchar(DFnoINDEL$pos)>6,]

  # extract DNP info
  df_3coppie <- righe_modificate %>%
    mutate(pos1 = substr(pos, 1, 2),
           pos2 = substr(pos, 3, 5),
           pos3 = substr(pos, 6, 7))
  df_3coppie <- subset(df_3coppie, select = -pos)

  df_count <- df_3coppie %>%
    group_by(chr, pos2) %>%
    summarise(n = n())

  df_count$chr <- str_extract(df_count$chr, "\\|.*\\|")
  df_count$chr <- substr(df_count$chr, 2, nchar(df_count$chr) - 1)

  data_info <- merge(df_count, RefAltDNP, by.x = "chr", by.y = "POS")

  df_grouped <- df_count %>%
    group_by(chr) %>%
    summarise(
      combined_bases = paste(unique(pos2), collapse = "-"),
      combined_counts = paste(n, collapse = "-"),
      total_count = sum(n)
    )

  # create df for reference alleles
  ref_dt <- subset(data_info, pos2 == ref)
  colnames(ref_dt) <- c("DNPs", "geno", "n_ref", "ref", "alt")

  # create df for alternative alleles
  alt_dt <- subset(data_info, pos2 == alt)
  colnames(alt_dt) <- c("DNPs", "geno", "n_alt", "ref", "alt")

  # create df for other (error) alleles
  altro_dt <- subset(data_info, pos2 != ref & pos2 != alt)
  altro_dt2 <- data.frame(DNPs = altro_dt$chr, n = altro_dt$n)
  altro_dt3 <- altro_dt2 %>% group_by(DNPs) %>% summarise(n_altro = sum(n))

  # merge the 3 dfs
  merged1 <- merge(ref_dt, alt_dt, by = "DNPs", all = TRUE)
  merged2 <- merge(merged1, altro_dt3, by = "DNPs", all = TRUE)
  merged3 <- merge(merged2, RefAltDNP, by.x = "DNPs", by.y = "POS", all = TRUE)

  merged3[is.na(merged3)] <- 0
  dataframe_final <- merged3[, c(1, 11, 12, 3, 7, 10)]
  colnames(dataframe_final) <- c("DNPs", "Reference", "Alternative", "N_Reference", "N_Alternative", "N_Altro")

  # compute ratio reference allele/total for genotyping
  dataframe_final <- dataframe_final %>%
    mutate(ratio_Ref=ifelse(N_Altro/(N_Altro+N_Reference+N_Alternative) >= 0.1, NA, 
                             N_Reference / (N_Reference + N_Alternative))) %>%
    select(everything(), ratio_Ref)
  
  
  


  FINAL_df <- merge(dataframe_final, df_grouped, by.x = "DNPs", by.y = "chr", all = TRUE)

  FINAL_df2 <- separate(FINAL_df, col = "DNPs", into = c("pos1", "pos2"), sep = "\\|")
  FINAL_df3 <- separate(FINAL_df2, col = "pos1", into = c("chr", "pos1"), sep = "\\:")
  FINAL_df4 <- separate(FINAL_df3, col = "pos2", into = c("chr", "pos2"), sep = "\\:")
  final <- FINAL_df4[,c(2,1,3:12)]

  final$Reference <- gsub('\\|', '', final$Reference)
  final$Alternative <- gsub('\\|', '', final$Alternative)
  final$combined_bases <- gsub('\\|', '', final$combined_bases)

  colnames(final) <- c("chr","pos1","pos2",
                       "REF", "ALT", "N_REF",
                       "N_ALT","N_Other",
                       "allele_ratio","combined_bases",
                       "combined_counts","total_count")

  indiv_final <- final[, c("chr", "pos1", "pos2", "REF", "ALT", "allele_ratio")]
  colnames(indiv_final)[6] <- paste0("GT_", i)

  # start making simil-VCF file
  if (is.null(final_merged_data)) {
    final_merged_data <- indiv_final
  } else {
    final_merged_data <- merge(final_merged_data, indiv_final,
                               by = c("chr", "pos1", "pos2", "REF", "ALT"),
                               all = TRUE)}
  
  # save individual output files
  write.table(final, paste0("full_output_", i, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)

}

# finalise simil-VCF file substituting ratio ref with genotype
final_merged_data <- final_merged_data %>%
  mutate_at(vars(starts_with("GT_")),
            ~ ifelse(. <= 0.1, "1/1",
                     ifelse(. >= 0.4 & . <= 0.6, "0/1",  # change here for more strict/relaxed heterozygous ratio
                            ifelse(. >= 0.9, "0/0", "NA"
                                   ))))

# save simil-VCF file
write.table(final_merged_data, "AllGenotypes.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# remove now useless files
commB = paste0("rm mpileup.txt")
system(commB, ignore.stderr=T)

commC = paste0("rm bedfile_for_mpileup.bed")
system(commC, ignore.stderr=T)

cat("\nAll done!\n")
