# Trace element and REE plotter #
# version 0.3.1  
# 27/10/2017
# basic edition, removed apatite specific discrimination plot code
# group and age now different variables to allow non-unique ages


# library path
.libPaths('D:/r/lib')

#run these first time on a computer
install.packages("readr")
install.packages("tibble")
install.packages("tidyr")
install.packages("dplyr")
install.packages("qdap")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("scales")
install.packages("viridis")

# run these every time you open R
library(readr)
library(tibble)
library(tidyr)
library(dplyr)
library(qdap)
library(ggplot2)
library(reshape2)
library(scales)
library(viridis)


#values for CL chondrite from Sun and McDonough 1989

Cs <- 0.188
Tl <- 0.14
Rb <- 2.32
Ba <- 2.41
W <- 0.095
Th <- 0.029
U <- 0.008
Nb <- 0.246
Ta <- 0.014
K <- 545
La <- 0.237
Ce <- 0.612
Pb <- 2.47
Pr <- 0.095
Mo <- 0.92
Sr <- 7.26
P <- 1220
Nd <- 0.467
F <- 60.7
Sm <- 0.153
Zr <- 3.87
Hf <- 0.1066
Eu <- 0.058
Sn <- 1.72
Sb <- 0.16
Ti <- 445
Gd <- 0.2055
Tb <- 0.0374
Dy <- 0.2540
Li <- 1.57
Y <- 1.57
Ho <- 0.0566
Er <- 0.1655
Tm <- 0.0255
Yb <- 0.170
Lu <- 0.0254

# working directory
setwd('D:/r/Yen final plots')


# import trace element data
# specify sample name - identical to that found in 'source.file' iolite output
unknown.name <- "YP4"

iolitedata <- read.delim("allyentrace.txt")

trace.data <- filter(iolitedata, grepl(paste(unknown.name), Source.file))


# import table with age data - can use for 206/238 age just as for 207 corrected age. 
# format - column titled 'source.file', identical grain names for trace element and U-pb iolite output 
# and second column titled 'age'
age.data <- read.csv("rough207ages1.csv")

age.data <- filter(age.data, grepl(paste(unknown.name), Source.file))

# merges trace and age sheets
all.data <- inner_join(trace.data, age.data, by = "Source.file")



unknowns <- all.data


# concentration data
# quality control to ensure values are numeric
unknowns$Mn_ppm_m55 <- as.numeric(as.character(unknowns$Mn_ppm_m55))

unknowns$Sr_ppm_m88 <- as.numeric(as.character(unknowns$Sr_ppm_m88))

unknowns$La_ppm_m139 <- as.numeric(as.character(unknowns$La_ppm_m139))

unknowns$Ce_ppm_m140 <- as.numeric(as.character(unknowns$Ce_ppm_m140))

unknowns$Pr_ppm_m141 <- as.numeric(as.character(unknowns$Pr_ppm_m141))

unknowns$Nd_ppm_m146 <- as.numeric(as.character(unknowns$Nd_ppm_m146))

unknowns$Sm_ppm_m147 <- as.numeric(as.character(unknowns$Sm_ppm_m147))

unknowns$Eu_ppm_m153 <- as.numeric(as.character(unknowns$Eu_ppm_m153))

unknowns$Gd_ppm_m157 <- as.numeric(as.character(unknowns$Gd_ppm_m157))

unknowns$Tb_ppm_m159 <- as.numeric(as.character(unknowns$Tb_ppm_m159))

unknowns$Dy_ppm_m163 <- as.numeric(as.character(unknowns$Dy_ppm_m163))

unknowns$Y_ppm_m89 <- as.numeric(as.character(unknowns$Y_ppm_m89))

unknowns$Ho_ppm_m165 <- as.numeric(as.character(unknowns$Ho_ppm_m165))

unknowns$Er_ppm_m166 <- as.numeric(as.character(unknowns$Er_ppm_m166))

unknowns$Tm_ppm_m169 <- as.numeric(as.character(unknowns$Tm_ppm_m169))

unknowns$Yb_ppm_m172 <- as.numeric(as.character(unknowns$Yb_ppm_m172))

unknowns$Lu_ppm_m175 <- as.numeric(as.character(unknowns$Lu_ppm_m175))

unknowns$Th_ppm_m232 <- as.numeric(as.character(unknowns$Th_ppm_m232))

unknowns$U_ppm_m238 <- as.numeric(as.character(unknowns$U_ppm_m238))

unknowns$U_ppm_m238_Int2SE <- as.numeric(as.character(unknowns$U_ppm_m238_Int2SE))


# edit this to specify the isotopes you are interested in/have measured
unknowns.conc <- select(unknowns, Source.file, Mn_ppm_m55, Sr_ppm_m88, Y_ppm_m89, La_ppm_m139, Ce_ppm_m140, Pr_ppm_m141, 
                    Nd_ppm_m146, Sm_ppm_m147, Eu_ppm_m153, Gd_ppm_m157,
                    Tb_ppm_m159, Dy_ppm_m163, Y_ppm_m89, Ho_ppm_m165, 
                    Er_ppm_m166, Tm_ppm_m169, Yb_ppm_m172, Lu_ppm_m175, Th_ppm_m232, U_ppm_m238, age)



# normalise REE values
# quality control to ensure values are numeric

unknowns$La_ppm_m139 <- as.numeric(as.character(unknowns$La_ppm_m139))
unknowns$La_ppm_m139 <- unknowns$La_ppm_m139/La

unknowns$Ce_ppm_m140 <- as.numeric(as.character(unknowns$Ce_ppm_m140))
unknowns$Ce_ppm_m140 <- unknowns$Ce_ppm_m140/Ce

unknowns$Pr_ppm_m141 <- as.numeric(as.character(unknowns$Pr_ppm_m141))
unknowns$Pr_ppm_m141 <- unknowns$Pr_ppm_m141/Pr

unknowns$Nd_ppm_m146 <- as.numeric(as.character(unknowns$Nd_ppm_m146))
unknowns$Nd_ppm_m146 <- unknowns$Nd_ppm_m146/Nd

unknowns$Sm_ppm_m147 <- as.numeric(as.character(unknowns$Sm_ppm_m147))
unknowns$Sm_ppm_m147 <- unknowns$Sm_ppm_m147/Sm

unknowns$Eu_ppm_m153 <- as.numeric(as.character(unknowns$Eu_ppm_m153))
unknowns$Eu_ppm_m153 <- unknowns$Eu_ppm_m153/Eu

unknowns$Gd_ppm_m157 <- as.numeric(as.character(unknowns$Gd_ppm_m157))
unknowns$Gd_ppm_m157 <- unknowns$Gd_ppm_m157/Gd

unknowns$Tb_ppm_m159 <- as.numeric(as.character(unknowns$Tb_ppm_m159))
unknowns$Tb_ppm_m159 <- unknowns$Tb_ppm_m159/Tb

unknowns$Dy_ppm_m163 <- as.numeric(as.character(unknowns$Dy_ppm_m163))
unknowns$Dy_ppm_m163 <- unknowns$Dy_ppm_m163/Dy

unknowns$Y_ppm_m89 <- as.numeric(as.character(unknowns$Y_ppm_m89))
unknowns$Y_ppm_m89 <- unknowns$Y_ppm_m89/Y

unknowns$Ho_ppm_m165 <- as.numeric(as.character(unknowns$Ho_ppm_m165))
unknowns$Ho_ppm_m165 <- unknowns$Ho_ppm_m165/Ho

unknowns$Er_ppm_m166 <- as.numeric(as.character(unknowns$Er_ppm_m166))
unknowns$Er_ppm_m166 <- unknowns$Er_ppm_m166/Er

unknowns$Tm_ppm_m169 <- as.numeric(as.character(unknowns$Tm_ppm_m169))
unknowns$Tm_ppm_m169 <- unknowns$Tm_ppm_m169/Tm

unknowns$Yb_ppm_m172 <- as.numeric(as.character(unknowns$Yb_ppm_m172))
unknowns$Yb_ppm_m172 <- unknowns$Yb_ppm_m172/Yb

unknowns$Lu_ppm_m175 <- as.numeric(as.character(unknowns$Lu_ppm_m175))
unknowns$Lu_ppm_m175 <- unknowns$Lu_ppm_m175/Lu


# extracts only those columns necessary for creating REE profile plot
unknowns.norm <- select(unknowns, Source.file, age, La_ppm_m139, Ce_ppm_m140, Pr_ppm_m141, 
                    Nd_ppm_m146, Sm_ppm_m147, Eu_ppm_m153, Gd_ppm_m157,
                    Tb_ppm_m159, Dy_ppm_m163, Y_ppm_m89, Ho_ppm_m165, 
                    Er_ppm_m166, Tm_ppm_m169, Yb_ppm_m172, Lu_ppm_m175)


# alters columns names for nicer display on plot
colnames(unknowns.norm) <- c("Source.file", "age", "La", "Ce", "Pr", 
                             "Nd", "Sm", "Eu", "Gd",
                             "Tb", "Dy", "Y", "Ho", 
                             "Er", "Tm", "Yb", "Lu")




## REE profile plot ##

# gather data for spider plot, set category to display on legend as age of grain
gather.norm <- gather(unknowns.norm, element, value, -age, -Source.file)



# create spider plot

# edit size, stroke of text, line, and point to suit your needs

gather.norm %>%
  ggplot(aes(x=element, y=value, group = Source.file, colour = age)) +
  geom_line(size=1) +
  geom_point(stroke = 1, size=1, shape=21, fill="white", colour = "black") +
  scale_y_log10("Apatite / Chondrite", limits= c(0.1, 10000)) +
  theme(axis.title.x = element_blank(),text = element_text(size=10)) +
  scale_x_discrete(limits=c("La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd","Tb", "Dy", "Y", "Ho","Er", "Tm", "Yb", "Lu")) +
  labs(colour="Age (Ma)") +
  scale_colour_viridis(option = 'plasma', direction = -1, limits=c(500,2500))

#save spider plot - change plot size if necessary
ggsave(paste0(unknown.name, "spiderplot.pdf"), width = 12, height = 9)

