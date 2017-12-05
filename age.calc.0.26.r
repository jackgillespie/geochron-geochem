#### Apatite fission track age calculation ####


# 11/03/2017 #
# v.0.26 #

# manual input of unknowns
# simple input and extraction of durango and automatic zeta factor calculation
# implemented fasttracks-iolite row matching for both unknowns and durango 
# name input as variable to generate output .csv and dataframes for eventual merging for all session unknowns
# set path to reduce mistakes
# moved all inputs to one area to increase efficiency

# working version



## further development plans ##

#create loop to input all unknowns
#loop calculation
#loop writing output tables
#radial plots for each sample
#output final table for publication




## References ##

# Age calculation after Vermeesch (2017) by zeta-calibration method

# Session-specific zeta-calibration calculated  following Hasebe et al (2004), Vermeesch (2017) 
# using the Durango apatite age standard (McDowell et al 2005)

# zeta uncertainty calculation following Donelick (2005) and Hasebe et al (2013)







### Start ###


#ensure 'readr', 'tibble', 'tidyr', 'dplyr', 'IsoplotR', 'xlsx', 'qdap' packages are installed via installed.packages()

#library path                                                                                     &&&&&&
.libPaths('D:/r/lib')


#library(xlsx)
library(readr)
library(tibble)
library(tidyr)
library(dplyr)
library(IsoplotR)
library(qdap)
library(ggplot2)
library(reshape2)


### Constants ###



#238U total decay constant
ld <- 1.55125*10^-10

#238U spontaneous fission decay constant
lf <- 8.46*10^-17

#half of the mean etchable spontaneous fission-track length in unannealed apatite
Rsp <- 0.000725

#mass of 238U
M <- 238.0508

#apatite density (gleadow et al 2015 after carlson 1999)
d <- 3.21

#k (hasebe 2004, 2013)
k <- 1

#efficiency factor (gleadow et al 2015)
q <- 0.96

#avogadro's number
Na <- 6.0221409*10^23

#Durango reference age (McDowell et al 2005)
tstd <- 31.44

#Durango reference age SD (McDowell et al 2005)
tstd.se <- 0.18

# absolute calibration
abs.cal <- (M/(lf*Na*d*Rsp*q))




### Age Equations ###



## Vermeesch 2016 zeta approach ##

#zeta calculation - pooled durango 
ver16.zeta <- function(DUR.Ns, DUR.area, DUR.uconc){
  (((sum(DUR.area)*mean(DUR.uconc))/(sum(DUR.Ns)*ld))*(exp(ld*tstd)-1))
}
ver16.zetasd <- function(zeta, DUR.Ns, DUR.area, DUR.uconc, DUR.uconc.se){
  zeta*sqrt((1/(sum(DUR.Ns)))+(sum((mean(DUR.area)/sum(DUR.area))*mean(DUR.uconc.se))^2/mean(DUR.uconc)^2)+((tstd.se/tstd)^2))
}

#age calculation 
ver16z <- function(track.density, uconc, zeta){
  (1/ld)*log(1+((ld*zeta*(track.density/uconc))))
}
ver16zsd <- function(uncor.age, Ns, uconc, uconc.se, zeta, zeta.SD){
  uncor.age*sqrt((1/(Ns))+((uconc.se/uconc)^2)+((zeta.SD/zeta)^2))
}


## Hasebe 2004 ##

has04 <- function(track.density, uconc){
  (1/ld)*log(1+((track.density*ld*M)/(lf*Na*uconc*(10^-6)*d*Rsp*k)))/(10^6)
}
has04sd <- function(uncor.age, Ns, uconc, unconc.se){
  uncor.age*sqrt((1/(Ns))+((unconc.se/uconc)^2))
}



## pooled age equation ## - laicpms version after Hasebe 2013
pooled.age <- function(zeta, Ns, area, uconc){
  (1/ld)*log(1+(ld)*zeta*((sum(Ns))/(sum(area)*mean(uconc))))
}



## central age ## - doesnt work
#central.age <- function(z, z.se){
#  (sum(z/z.se^2)/(sum(1/z.se^2)))
#}







### INPUTS ###


#set your own working directory and library path                                                  &&&&&&
path = "D:/r/pb corr/"#                                                                       &&&&&&

setwd(path)






## import durango ##
# input counts file name                           

DUR.counts <- read.csv(paste0(path, "DURf count data.csv"), skip = 4, stringsAsFactors = FALSE)


# import all iolite files for zeta calculation  - extract durangos #
# input directory                                                                                

#searches for all iolite outputs in folder following naming convention below
all.iolite<-""

files <- dir(path, pattern = glob2rx("***** AFT_All_Integrations.txt"))


#edit this to output durango table for entire session once calculated
write.csv(DUR.comb, 'DURj.csv')





## input sample name ##
unknown.name <- '8'

#import raw iolite and count data tables

iolite <- read.delim(paste0(path, "YP-8 MAD_All_Integrations.txt"))
counts <- read.csv(paste0(path, "YP-8 count data.csv"), skip = 4, stringsAsFactors = FALSE)








### Session Durango ###


#create dataframe of all iolite outputs
for(i in 1:length(files)) {
  file <- read.delim(files[i],head=T,stringsAsFactors = FALSE)
  all.iolite <- rbind(all.iolite, file)
}

# extract durangos from combined file
DUR.iolite <- filter(all.iolite, grepl('DUR', Source.file))

# remove empty rows
DUR.counts <- filter(DUR.counts, Area.cm2. > 0)

# remove all characters and set source.file to numeric to allow correct row matching
DUR.iolite$Source.file <- as.numeric(as.character(gsub("\\D", "", DUR.iolite$Source.file)))

# remove all characters and set grain.mica to numeric to allow correct row matching
DUR.counts$Grain.Mica <- as.numeric(as.character(gsub("\\D", "", DUR.counts$Grain.Mica)))

# remove characters and set source.file to numeric to allow correct sorting
#DUR.iolite$Source.file <- as.numeric(as.character(genXtract(DUR.iolite$Source.file, "R", ".")))


# order dataframe by durango number
DUR.iolite <- DUR.iolite[order(DUR.iolite$Source.file),]

#combine iolite and count data - matches by row values - can handle non-matching number of rows and different row order
DUR.comb <- left_join(DUR.iolite, DUR.counts, by=c("Source.file"="Grain.Mica"))

#ensure number, not factor so that age equation works
DUR.comb$U_ppm_m238_Int2SE <- as.numeric(as.character(DUR.comb$U_ppm_m238_Int2SE))
DUR.comb$U_ppm_m238 <- as.numeric(as.character(DUR.comb$U_ppm_m238))

DUR.comb$Th_ppm_m232_Int2SE <- as.numeric(as.character(DUR.comb$Th_ppm_m232_Int2SE))
DUR.comb$Th_ppm_m232 <- as.numeric(as.character(DUR.comb$Th_ppm_m232))

DUR.comb$Cl_ppm_m35_Int2SE <- as.numeric(as.character(DUR.comb$Cl_ppm_m35_Int2SE))
DUR.comb$Cl_ppm_m35 <- as.numeric(as.character(DUR.comb$Cl_ppm_m35))

DUR.comb <- filter(DUR.comb, Area.cm2. > 0)

#convert to numeric preserves exponent 
DUR.comb$Density.tracks.cm2. <- as.numeric(as.character(DUR.comb$Density.tracks.cm2.))






### zeta calculation ###



#calculate pooled age 
DURf.pooled <- pooled.age(abs.cal, DUR.comb$Tracks, DUR.comb$Area.cm2., DUR.comb$U_ppm_m238)

#age calc using function
DUR.comb <- mutate(DUR.comb, has04.Uncor.Age = has04(Density.tracks.cm2., U_ppm_m238))

#calculate sd - equation after Hasebe et al. 2004
DUR.comb <- mutate(DUR.comb, has04.SD = has04sd(has04.Uncor.Age, Tracks, U_ppm_m238, U_ppm_m238_Int2SE))

# calculate zeta values

DUR.zeta <- ver16.zeta(DUR.comb$Tracks, DUR.comb$Area.cm2., DUR.comb$U_ppm_m238)

DUR.zeta.SD <- ver16.zetasd(DUR.zeta, DUR.comb$Tracks, DUR.comb$Area.cm2., DUR.comb$U_ppm_m238, DUR.comb$U_ppm_m238_Int2SE)






### Unknowns ###

## input sample name ##                                                                               &&&&&&&&&&&
unknown.name <- '8'

#import raw iolite and count data tables

iolite <- read.delim(paste0(path, "YP-8 MAD_All_Integrations.txt"))
counts <- read.csv(paste0(path, "YP-8 count data.csv"), skip = 4, stringsAsFactors = FALSE)


## data manipulation ##

#extract from iolite separate data files for unknowns and each standard
iolite_unknowns <- filter(iolite, grepl(unknown.name, Source.file), Ca43_CPS > 1)

#iolite_MAD <- filter(iolite, grepl('MAD', Source.file))

#iolite_DUR <- filter(iolite, grepl('DUR', Source.file))

#iolite_McC <- filter(iolite, grepl('McC|McL', Source.file))

#separate from iolite important variables - probably unnecessary in hindsight
#iolite_unknowns_var <- select(iolite_unknowns, Source.file, Cl_ppm_m35,Cl_ppm_m35_Int2SE, Th_ppm_m232, Th_ppm_m232_Int2SE, U_ppm_m238, U_ppm_m238_Int2SE)

#iolite_MAD_var <- select(iolite_MAD, Source.file, Cl_ppm_m35,Cl_ppm_m35_Int2SE, Th_ppm_m232, Th_ppm_m232_Int2SE, U_ppm_m238, U_ppm_m238_Int2SE)

#iolite_DUR_var <- select(iolite_DUR, Source.file, Cl_ppm_m35,Cl_ppm_m35_Int2SE, Th_ppm_m232, Th_ppm_m232_Int2SE, U_ppm_m238, U_ppm_m238_Int2SE)

#iolite_McC_var <- select(iolite_McC, Source.file, Cl_ppm_m35,Cl_ppm_m35_Int2SE, Th_ppm_m232, Th_ppm_m232_Int2SE, U_ppm_m238, U_ppm_m238_Int2SE)



## data cleaning + combine track count and LAICPMS data ##

# ensure only grains not lengths
counts <- filter(counts, grepl('Grain', Grain.Mica))

# remove characters and set grain.mica to numeric to allow correct sorting
counts$Grain.Mica <- as.numeric(as.character(gsub("\\D", "", counts$Grain.Mica)))

# remove characters and set source.file to numeric to allow correct sorting
iolite_unknowns$Source.file <- as.numeric(as.character(genXtract(iolite_unknowns$Source.file, "_", ".")))

#combine iolite and count data
comb <- left_join(iolite_unknowns, counts, by=c("Source.file"="Grain.Mica"))

#ensure number, not factor so that age equation works
comb$U_ppm_m238 <- as.numeric(as.character(gsub("no value", "0", comb$U_ppm_m238))) 
comb$U_ppm_m238_Int2SE <- as.numeric(as.character(gsub("NaN", "0", comb$U_ppm_m238_Int2SE)))

comb$Th_ppm_m232 <- as.numeric(as.character(gsub("no value", "0", comb$Th_ppm_m232)))
comb$Th_ppm_m232_Int2SE <- as.numeric(as.character(gsub("NaN", "0", comb$Th_ppm_m232_Int2SE)))


# remove uncounted grains - after combine so also removes LAICPMS values and ensure grain matching data
comb <- filter(comb, Area.cm2. > 0)
comb <- filter(comb, U_ppm_m238 > 0)

# remove this line once zero track solution implemented
comb <- filter(comb, Tracks > 0)


#convert to numeric preserves exponent 
comb$Density.tracks.cm2. <- as.numeric(as.character(comb$Density.tracks.cm2.))

#if using cl for radplots - do not use for final output
#levels(comb$Cl_ppm_m35)[levels(comb$Cl_ppm_m35)=='Below LOD'] <- '0'
#levels(comb$Cl_ppm_m35_Int2SE)[levels(comb$Cl_ppm_m35_Int2SE)=='Below LOD'] <- '0'




### AGE CALCULATION ###



#age calc using function
comb <- mutate(comb, has04.Uncor.Age = has04(Density.tracks.cm2., U_ppm_m238))

#calculate sd - equation after Hasebe et al. 2004
comb <- mutate(comb, has04.SD = has04sd(has04.Uncor.Age, Tracks, U_ppm_m238, U_ppm_m238_Int2SE))


#age calc using function
comb <- mutate(comb, ver16z.Age = ver16z(Density.tracks.cm2., U_ppm_m238, DUR.zeta))

#calculate sd 
comb <- mutate(comb, ver16z.SD = ver16zsd(ver16z.Age, Tracks, U_ppm_m238, U_ppm_m238_Int2SE, DUR.zeta, DUR.zeta.SD))




### Output table format ###





# results output and automatically named
write.csv(comb, paste0(unknown.name, '.csv')) 













##shit's fucked, yo##


# select columns for supplementary info output
#comb.output <- select(comb, Source.file, Tracks, Area.cm2., Density.tracks.cm2., U_ppm_m238, U_ppm_m238_Int2SE, Cl_ppm_m35, Cl_ppm_m35_Int2SE, Average.DPar.µmm., ver16z.Age, ver16z.SD)


# summarise columns of interest for publication table #

# format columns to make suitable for summary
comb$Average.DPar.µmm. <- as.numeric(as.character(gsub("No DPars", "", comb$Average.DPar.µmm.)))

# Cl conc. currently as factor - how to deal with Below LOD values??

# also needs column for: sample name, lat/long, altitude, lithology, analyst, chisquared, central age and uncertainty,
# pooled age and uncertainty, Nlength, MTL, Std.devlength

#calculate central age
#comb <- mutate(comb, z = log(ver16z.Age))
#comb <- mutate(comb, z0 = central.age(z, z.se))


summary.comb <- summarise(comb, n(), sum(Tracks), mean(Density.tracks.cm2.), mean(U_ppm_m238), mean(Cl_ppm_m35, na.rm=TRUE), mean(Average.DPar.µmm., na.rm = TRUE))



# results output and automatically named
write.csv(comb, paste0(unknown.name, '.csv')) 

# creates automatically named tables for supplementary info and for published table
assign(paste0(unknown.name, ".supp"),comb.output)

assign(paste0(unknown.name, ".final"), summary.comb)

# 




