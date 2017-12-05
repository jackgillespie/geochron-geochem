setwd('D:/r/biryusa')
#import raw data table from fastracks output
lengths <- read.csv('D:/r/biryusa/14-45 lengths.csv', skip =  4)

#define mean, sd, n
Mean <- mean(lengths$True.Length)
SD <- sd(lengths$True.Length)
n <- nrow(lengths)

#create histogram of 'True Length' column
h= hist(lengths$True.Length, main = 'Track length distribution', 
        col='grey', breaks=seq(6,17, by = 1), xlim=c(6,18), ylab = 'Frequency', xlab = 'Track length (µm)')

#convert to relative frequency
h$density=h$counts/sum(h$counts)
plot(h, freq=FALSE, main = 'Track length distribution', 
     col='grey', xlim=c(6,18), ylim=c(0,0.4), ylab = 'Frequency', xlab = 'Track length (µm)', las=1)

#add text - use legend
legend("topleft", inset = 0.05, legend = c(paste("Mean =", round(Mean, 1)),
                                           paste("Std. Dev. =", round(SD, 1)),
                                           paste("n =", n)), bty = 'n')

#save as pdf
#remember to update for each sample
dev.copy(pdf,'D:/r/biryusa/14-45 length.pdf')
dev.off()
