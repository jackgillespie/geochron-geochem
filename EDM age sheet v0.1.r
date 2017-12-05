# EDM zeta calibration

# v0.1


setwd('D:/r/Tunka EDM')
.libPaths('D:/r/lib')


#238U total decay constant
ld <- 1.55125*10^-10

#Mt Drom reference age (Williams et al 1982)
tstd <- 98.8

#Mt Drom reference age SD (Williams et al 1982)
tstd.se <- 0.6

theta <- read.csv("tunka theta interpol.csv")
cn513 <- read.csv('CN5-13 count data.csv', skip = 4, stringsAsFactors = FALSE)
mt.drom <- read.csv('Mt Drom count data.csv', skip = 4, stringsAsFactors = FALSE)
unknown <- read.csv('LS-4 count data.csv', skip = 4, stringsAsFactors = FALSE)#         &&&&&&


##  functions ##


edm.zeta <- function(rhos, rhoi, rhod){
  (exp(ld*tstd)-1)/(ld*(mean(rhos)/mean(rhoi))*0.5*rhod)
}

edm.zetasd <- function(zeta, Ns, Ni, Nd){
  zeta*sqrt((1/sum(Ns))+(1/sum(Ni))+(1/sum(Nd))+((tstd.se/tstd)^2))
}

edm.age <- function(rhos, rhoi, zeta, rhod){
  (1/ld)*log(1+(ld*zeta*0.5*(rhos/rhoi)*rhod))
}

edm.agesd <- function(age, zeta, zeta.se, Ns, Ni, Nd){
  age*sqrt((1/Ns)+(1/Ni)+(1/sum(Nd))+((zeta.se/zeta)^2))
}


## calculations ##

#zeta
jack.zeta <- edm.zeta(mt.drom$Density.tracks.cm2., mt.drom$Mica.Density.tracks.cm2., theta[12,2])


jack.zetase <- edm.zetasd(jack.zeta, mt.drom$Tracks, mt.drom$Mica.Count, cn513$Tracks)


#age
unknown <- mutate(unknown, age = edm.age(unknown$Density.tracks.cm2., unknown$Mica.Density.tracks.cm2., jack.zeta, theta[2,2]))

unknown <- mutate(unknown, age.se = edm.agesd(unknown$age, jack.zeta, jack.zetase, unknown$Tracks, unknown$Mica.Count, cn513$Tracks))


#write output
write.csv(unknown, "LS-4 output.csv")

