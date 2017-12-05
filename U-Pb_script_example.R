

## Apatite U-Pb concordia and WM plot generator ##

# version 1.0
# 25/09/2017



#set library - packages are installed here
.libPaths('D:/r/lib')

#install packages (only necessary first time on a particular computer)
install.packages("readxl")
install.packages("readr")
install.packages("tibble")
install.packages("tidyr")
install.packages("dplyr")
install.packages("IsoplotR")
install.packages("qdap")
install.packages("ggplot2")
install.packages("reshape2")

#activate libraries (every time you start R)
library(readxl)
library(readr)
library(tibble)
library(tidyr)
library(dplyr)
library(IsoplotR)
library(qdap)
library(ggplot2)
library(reshape2)
library(xlsx)


# set working directory - access files, save files here by default
setwd("E:/")



path = "D:/r/stijn/ALL"



####set name of sample (identical to excel sheet) for data import, plot generation
unknown.name <- "UZ_39"
####

sheet <- read.delim("ZU-Pb_example.txt")


# filter rows so only unknowns remain
unknowns.sheet <- filter(sheet, !grepl("MAD|McL|McC|N610|610|DUR", Comments))

# remove empty rows
unknowns.sheet <- unknowns.sheet[!is.na(unknowns.sheet$Total.points),]

# select columns with values of interest for plotting concordia 
plotting.sheet <- select(unknowns.sheet, Final207_235, Final207_235_Prop2SE, Final206_238, Final206_238_Prop2SE, ErrorCorrelation_6_38vs7_35)

# divide errors by two to make 1se 
plotting.sheet$Final207_235_Int2SE <- plotting.sheet$Final207_235_Prop2SE/2
plotting.sheet$Final206_238_Int2SE <- plotting.sheet$Final206_238_Prop2SE/2

# divide errors by two to make 1se 
sheet$Final207_235_Prop2SE <- sheet$Final207_235_Prop2SE/2
sheet$Final206_238_Prop2SE <- sheet$Final206_238_Prop2SE/2

# converts to matrix so isoplotr can read it
plotting.matrix <- data.matrix(sheet)



## use this line to remove data points by row (note: grain number n on concordia is row n+1) ##

################
plotting.matrix <- plotting.matrix[-c(2,22),]
################
# view data
print(plotting.matrix)





# puts data into isoplot format
data.isoplot <- read.data(plotting.matrix, method = "U-Pb", format = 1)                          



# plot concordia diagram
concordia(data.isoplot, wetherill=TRUE, show.numbers = FALSE, tlim = c(250,2100))

# saves most recently generated plot to working directory
dev.copy(pdf, paste0(unknown.name, 'concordia.pdf'))
dev.off()

# plot KDE 
kde(data.isoplot, binwidth = 20, log = FALSE, type = 4, cutoff.76 = 1000)

# saves most recently generated plot to working directory
dev.copy(pdf, paste0(unknown.name, 'KDE'))
dev.off()


age(data.isoplot, type = 3)


# calculates weighted mean age plot
weightedmean(data.isoplot, common.Pb = 1)

# saves most recently generated plot to working directory
dev.copy(pdf, paste0(unknown.name, 'WMplot.pdf'))
dev.off()


