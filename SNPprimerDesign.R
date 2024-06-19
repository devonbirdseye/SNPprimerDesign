# Find sequences for primers with degenerate bases according to selected panel of varieties. Input affymetrix of alleles for each variety within a selected region, position of SNP for which you want surrounding sequence (+/- 100bp), indication of whether it's plus or minus strand, and the reference sequence for this region.
# When you get sequences from GT they use the IUPAC alleles taking ALL of the lines into consideration. This script will get the IUPAC sequences taking only selected lines into consideration.
require(spgs)

############################################################################ INPUTS ##########################################################################################
affy <- "LT_DM_EZ_V10/9f1291de-5011-438a-bbe4-45d1a2167d8a-V10_119640122.matrix.affymetrix" #GT's affymetrix formated export for the region of the SNP plus/minus 100bp
snp_pos <- 119640122 #position of SNP (bp)
strand <- "minus" #indicate whether the sequence is on the plus or minus strand
#copy+paste the reference sequence of snp plus/minus 100bp (should be 201bp total). Be sure to copy the minus strand sequence if you indicated "minus" above
snp_seq <- "AAACAGGGGAAGATCTTAAACCTGAATTAGAGAGTGCAATTAAGGCATCTCGGGCTTCTATTATTGTGTTGTCTAAGAACTACGCTTCTTCAACATGGTGTCTGGATGAACTGGTGTTGATCCTTGAGCAGCGTATCACATCCAACCACATTGTCATCCCCATATTTTATCATGTGGAGCCCACCCATGTCAGGAAGCAAC"

#select the lines you want to pull alleles from
line_table <- read.csv("SKT_Lettuce_reseq100_v8_GTprojecLineInfo.csv") #import info table about lines in the resequencing panel on GT
#lines <- line_table$NAME
#lines <- line_table$NAME[line_table$Species %in% c("L. serriola")]
lines <- line_table$NAME[line_table$Crop.type %in% c("Cos", "Crisp")]
#lines <- line_table$NAME[line_table$Crop.type == "Crisp"]
#lines <- line_table$NAME[line_table$Crop.type %in% c("Cos", "Crisp")|line_table$Accession == "CGN04683"]
#lines <- line_table$NAME[line_table$Crop.type %in% c("Cos", "Crisp")|line_table$Species == "L. virosa"]
#lines <- line_table$NAME[line_table$Accession %in% c("Salinas", "CGN15705")]
##############################################################################################################################################################################
#After adjusting above inputs, highlight all code and run!


if(strand=="minus"){snp_seq <- toupper(paste0(spgs::reverseComplement(snp_seq), collapse = ""))}

#get affymetrix table from region of interest
snp <- read.table(affy, sep = "\t", header = T, tryLogical = F) #import GT's affymetrix
snp <- snp[,-c(1,3)] #remove contig and ranking columns
snp <- snp[,-grep("FREQUENCY|COVERAGE", colnames(snp))] #subset to just the allele info
rownames(snp) <- snp$Position
snp[snp==""] <- NA

#filter snp table to just the chosen lines
lines <- paste0(sub("-", ".", lines), "_ALLELE") #make line names match column names of snp table
snp <- snp[,colnames(snp) %in% lines]

#get the full sequence from the reference genome
snp_seq <- unlist(strsplit(snp_seq, "")) #convert ref sequence to vector
length(snp_seq)==201 #check that it's the right length

#make dataframe with each nucleoride as a row in the "seq" column and the bp position in the genome as the "bp" column
snp_df <- data.frame("bp"=c((snp_pos-100):(snp_pos+100)), "ref"=snp_seq)
#merge full sequence with SNPs
snp$bp <- as.numeric(rownames(snp)) #create bp column to merge by
snp <- merge.data.frame(snp_df, snp, by = "bp", all.x = T)
snp <- snp[,-1] #remove bp column

#function to convert vector of nucleotides to iupac letter
allele2iupac <- function(allele){
  allele = as.character(na.omit(as.character(allele)))
  allele = unlist(strsplit(allele, "/")) #split heterozygous alleles so that both are considered
  if(all(unique(allele)=="A")){"A"}
  else if(all(unique(allele)=="T")){"T"}
  else if(all(unique(allele)=="C")){"C"}
  else if(all(unique(allele)=="G")){"G"}
  else if("A" %in% allele & "T" %in% allele & !("C" %in% allele) & !("G" %in% allele)){"W"}
  else if("C" %in% allele & "G" %in% allele & !("A" %in% allele) & !("T" %in% allele)){"S"}
  else if("A" %in% allele & "G" %in% allele & !("C" %in% allele) & !("T" %in% allele)){"R"}
  else if("C" %in% allele & "T" %in% allele & !("A" %in% allele) & !("G" %in% allele)){"Y"}
  else if("G" %in% allele & "T" %in% allele & !("C" %in% allele) & !("A" %in% allele)){"K"}
  else if("A" %in% allele & "C" %in% allele & !("T" %in% allele) & !("G" %in% allele)){"M"}
  else if("C" %in% allele & "G" %in% allele & "T" %in% allele & !("A" %in% allele)){"B"}
  else if("A" %in% allele & "G" %in% allele & "T" %in% allele & !("C" %in% allele)){"D"}
  else if("A" %in% allele & "C" %in% allele & "T" %in% allele & !("G" %in% allele)){"H"}
  else if("A" %in% allele & "C" %in% allele & "G" %in% allele & !("T" %in% allele)){"V"}
  else if("A" %in% allele & "C" %in% allele & "G" %in% allele & "T" %in% allele){"N"}
}

#convert each row (locus) to IUPAC allele
snp$iupac <- apply(snp, 1, allele2iupac)

#return sequence with iupac alleles
if(strand=="plus"){paste0(snp$iupac, collapse = "")}
if(strand=="minus"){toupper(spgs::reverseComplement(paste0(snp$iupac, collapse = ""), case="upper"))}