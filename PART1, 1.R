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
subset( originaldata, Meansofothercolumn < 10)
nrow