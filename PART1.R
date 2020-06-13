#Start of assignment- PART 1/ QUESTIONS 1 TO 5

# Downloaded the data of gene expression file - ANSWER 1
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part1_files/gene_expression.tsv",
              destfile = "assessment3.tsv")
# Reading the table and making sure that gene accession names are the row names
originaldata <- read.table("assessment3.tsv")

#ANSWER 1 
originaldata <- read.table("assessment3.tsv", header= TRUE, stringsAsFactors= FALSE, row.names=1) 
# Data imported properly and printing first six rows and structure as well.
head(originaldata)
str(originaldata)

#ANSWER 2
# Command rowmeans is used to find the means of other columns 
rowMeans(originaldata)
# A new column is made with the name of Meansofothercolumn to store these means
Meansofothercolumn <- rowMeans(originaldata)
# This column is bind to the main data now
cbind(originaldata, Meansofothercolumn)
# Permanent change to data. Adding the new column of means
originaldata <- cbind(originaldata, Meansofothercolumn)
# Observe the first six genes to check if new column is added
head(originaldata)

#ANSWER 3
# The order commands is used in order to adjust any of the columns in ascending or descending order
order( originaldata $ Meansofothercolumn)
# The order of column Meansofothercolumn is arranged in descending order with the highest values first
originaldata[order( - originaldata$Meansofothercolumn ) , ]
# Permanent change of the descending order to Meansofothercolumn column has been made after this step
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
hist(originaldata$Meansofothercolumn)
# But because the genes with mean say less than 100 are maxiumum the entire hishtogram was just plotted round that values.
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
# Observing the structure of the data and using command head to print the data of first 6 rows
head(originaldatapart2)
str(originaldatapart2)
# Command colnames will give us the result to check the column names.
colnames(originaldatapart2)

#ANSWER 7 
# by the command subset i would separate the data for both the sites: northeast and southwest
subset(originaldatapart2, Site== 'northeast')
# Making permanent changes to the sites
NORTHEAST <- subset(originaldatapart2, Site== 'northeast')
subset(originaldatapart2, Site== 'southwest')
SOUTHWEST <- subset(originaldatapart2, Site== 'southwest')
mean(NORTHEAST$Circumf_2010_cm)
mean(NORTHEAST$Circumf_2019_cm)
mean(SOUTHWEST$Circumf_2004_cm)
mean(SOUTHWEST$Circumf_2019_cm)
sd(NORTHEAST$Circumf_2004_cm)
sd(NORTHEAST$Circumf_2019_cm)
sd(SOUTHWEST$Circumf_2004_cm)
sd(SOUTHWEST$Circumf_2019_cm)
