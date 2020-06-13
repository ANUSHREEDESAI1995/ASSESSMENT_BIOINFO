#Start of assignment- PART 1/ QUESTIONS 1 TO 5

# Downloaded the data of gene expression file - ANSWER 1
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part1_files/gene_expression.tsv",
              destfile = "assessment3.tsv")
# Reading the table into R
originaldata <- read.table("assessment3.tsv")

#ANSWER 1 
# Making sure that gene accession names are the row names
originaldata <- read.table("assessment3.tsv", header= TRUE, stringsAsFactors= FALSE, row.names=1) 
# Data imported properly and printing first six rows and structure as well.
head(originaldata)
str(originaldata)

#ANSWER 2

# A new column is made with the name of Meansofothercolumn to store these means.
Meansofothercolumn <- rowMeans(originaldata)
# Adding the new column of means.
# Permanent change to data. 
originaldata <- cbind(originaldata, Meansofothercolumn)
# Observe the first six genes to check if new column is added
head(originaldata)

#ANSWER 3

# Permanent change of the descending order to Meansofothercolumn column has been made after arranging in descending step using the command order
originaldata <- originaldata[order( - originaldata$Meansofothercolumn ) , ]
# Extracting just the first ten rows with highest means
row.names(originaldata[1:10,])

#ANSWER 4
# Using the subset command to only extract the genes with means less than 10
# And simultaneously added nrow to count how many genes have means less than 10
nrow(subset(originaldata, Meansofothercolumn < 10))

#ANSWER 5
# In order to plot a histogram command hist is used. 
# i tried plotting it just normally first
hist(originaldata$Meansofothercolumn, xlab="Mean of the genes", 
     main = "MEAN V/S FREQUENCY")
# But because the genes with mean say less than 25,000 are maxiumum the entire hishtogram was just plotted round that values.
# Hence by limiting the genes with means in x axis from 20,000 to 25,000 
# and frequency till 60 we could get a good histogram
hist(originaldata$Meansofothercolumn, xlim= c(20000,250000), ylim= c(0,60), xlab="Mean of the genes", 
     main = "MEAN V/S FREQUENCY")

#START OF PART1/ QUESTIONS 6 TO 10
#ANSWER 6
# Downloading the file and reading it into R
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part1_files/growth_data.csv", 
              destfile="assessment3part2.csv")
originaldatapart2 <- read.csv("assessment3part2.csv")
# Command colnames will give us the result to check the column names.
colnames(originaldatapart2)

#ANSWER 7 
# by the command subset i would separate the data for both the sites: northeast and southwest
# Making permanent changes to the sites
NORTHEAST <- subset(originaldatapart2, Site== 'northeast')
SOUTHWEST <- subset(originaldatapart2, Site== 'southwest')
# Calculating the mean of both the sites using command mean, for years 2004 and 2019
# Calculating the standard deviation for both the sites using command sd, for years 2004 and years 2019
sd(NORTHEAST$Circumf_2004_cm)
sd(NORTHEAST$Circumf_2019_cm)
sd(SOUTHWEST$Circumf_2004_cm)
sd(SOUTHWEST$Circumf_2019_cm)

#ANSWER 8
# Making a box plot of tree circumference at both sites
boxplot ( NORTHEAST$Circumf_2004_cm, NORTHEAST$Circumf_2019_cm, names = c("year2004", "year2019"),
          ylab = "CIRCUMFERENCE IN CMS", main= "CIRCUMFERENCE OF YEAR 2004 AND 2019 OF NORTHEAST SITE")
boxplot ( SOUTHWEST$Circumf_2004_cm, SOUTHWEST$Circumf_2019_cm, names = c("year2004", "year2019"),
          ylab = "CIRCUMFERENCE IN CMS", main= "CIRCUMFERENCE OF YEAR 2004 AND 2019 OF NORTHEAST SITE")

#ANSWER 9
# Calcuating the mean for past 10  years i.e. 2009 and 2019.
mean(NORTHEAST$Circumf_2009_cm)
mean(NORTHEAST$Circumf_2019_cm)
mean(SOUTHWEST$Circumf_2009_cm)
mean(SOUTHWEST$Circumf_2019_cm)

# Obtaining the difference by subtracting circumference of year 2009 from year 2019 at site Northeast
DiffNortheast <- NORTHEAST$Circumf_2019_cm - NORTHEAST$Circumf_2009_cm
# Making permanent change by adding a new column DiffNorthEast that is the difference in circumference of year 2009 and 2019
NORTHEAST <-cbind(NORTHEAST,DiffNortheast)
# Checking first six rows
head(NORTHEAST)

# Obtaining the difference by subtracting circumference of year 2009 from year 2019 at site SouthWest
DiffSouthwest <- SOUTHWEST$Circumf_2019_cm - SOUTHWEST$Circumf_2009_cm
# Making permanent change by adding a new column DiffSouthWest that is the difference in circumference of year 2009 and 2019
SOUTHWEST <- cbind(SOUTHWEST,DiffSouthwest)
# Checking first six rows
head(SOUTHWEST)

# Mean of the difference of circumference of year 2009 and 2019 at site NorthEast
meanNorth <- mean(DiffNortheast)
# Stored as new variable meanNorth
meanNorth

# Mean of the difference of circumference of year 2009 and 2019 at site SouthWest
meanSouth <- mean(DiffSouthwest)
# Stored as new variable meanSouth
meanSouth

# ANSWER 10
# In order to find the P- value I used the new columns made DiffNorthEast and DiffSouthWest
t.test(DiffSouthwest,DiffNortheast)
wilcox.test(DiffNortheast,DiffSouthwest)

# PART 2 OF THE ASSIGNMENT STARTS HERE

#ANSWER 1
#Some of the basic libraries have been run to apply the functions and run the commands 
library("seqinr")
library("R.utils")
library("rBLAST")
library("ape")
library("ORFik")
library("Biostrings")

# There are some premade functions by our supervisor that have been run as well
source("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/mutblast_functions.R")

# The whole e.coli file has been run
download.file("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-42/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cds/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz",
              destfile = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz")

#Shows the gunzip() function to decompress the file.
R.utils:: gunzip("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz", overwrite="TRUE")

# makeblastdb() function is used to make a blast database and calculate number of genes.
makeblastdb("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa",dbtype = "nucl","-parse_seqids")

#ANSWER 2

# Downloading new file and naming it as sample.fa
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/sample.fa", destfile = "sample.fa" )
makeblastdb("sample.fa",dbtype = "nucl","-parse_seqids")
# A new sequence variable is made that stores sample.fa data
sequence <- read.fasta("sample.fa")
# Shows structure of my sequence
str(sequence)
# My number is 44 so extractiing 44 number sequence and storing it in new variable SEQUENCE
SEQUENCE <- sequence[[44]]
str(SEQUENCE)
# Using command getlength()
seqinr::getLength(SEQUENCE)
# Using command GC()
seqinr::GC(SEQUENCE)

# ANSWER 3
# myblastn_tab() command is used. It tries to find closest match to my SEQUENCE
res <- myblastn_tab(myseq = SEQUENCE, db = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
res
head(res)
hits<- res[1:3, c(3, 11, 12)]
hits

SEQUENCEMUTATE <- mutator(myseq=SEQUENCE,100)
str(SEQUENCEMUTATE)
SEQUENCE
SEQUENCEMUTATE
SEQUENCE <- DNAString(c2s(SEQUENCE))
SEQUENCEMUTATE <- DNAString(c2s(SEQUENCEMUTATE))
CHANGE <- Biostrings::pairwiseAlignment(SEQUENCE,SEQUENCEMUTATE)
pid(CHANGE)
nmismatch(CHANGE)
