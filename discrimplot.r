
### discrimination plots ###

# v0.1
# 5/12/17


# script for making discrimination plots after normalisation of data and REE profile plots in trace.elements script
# apatite



## Normalised values ##

# Eu anomaly
unknowns.norm <- mutate(unknowns.norm, Eua = ((Sm+Gd)/2))

unknowns.norm <- mutate(unknowns.norm, Eu.Eua = (Eu/Eua))

# light/heavy ratios
unknowns.norm <- mutate(unknowns.norm, Ce.Yb = (Ce/Yb))

unknowns.norm <- mutate(unknowns.norm, La.Nd = (La/Nd))


## Concentration values ##

# sum of rare earths

unknowns.conc <- mutate(unknowns.conc, sumREE = ((La_ppm_m139 + Ce_ppm_m140 + Pr_ppm_m141 + Nd_ppm_m146 + Sm_ppm_m147 + Eu_ppm_m153 + Gd_ppm_m157 + Tb_ppm_m159 + Dy_ppm_m163 + Y_ppm_m89 + Ho_ppm_m165 +Er_ppm_m166 + Tm_ppm_m169 + Yb_ppm_m172 + Lu_ppm_m175)/10000))

unknowns.conc <- mutate(unknowns.conc, lights.sumREE = ((La_ppm_m139 + Ce_ppm_m140 + Pr_ppm_m141)/((La_ppm_m139 + Ce_ppm_m140 + Pr_ppm_m141 + Nd_ppm_m146 + Sm_ppm_m147 + Eu_ppm_m153 + Gd_ppm_m157 + Tb_ppm_m159 + Dy_ppm_m163 + Y_ppm_m89 + Ho_ppm_m165 +Er_ppm_m166 + Tm_ppm_m169 + Yb_ppm_m172 + Lu_ppm_m175)/100)))



# extract and combine commonly plotted values

unknowns.norm.ratios <- select(unknowns.norm, age, Eu.Eua, Ce.Yb, La.Nd)

unknowns.conc.sum <- select(unknowns.conc, age, Mn_ppm_m55, Sr_ppm_m88, Y_ppm_m89, Th_ppm_m232, U_ppm_m238, sumREE, lights.sumREE)

unknowns.plot <- inner_join(unknowns.norm.ratios, unknowns.conc.sum, by = "age")



write.csv(unknowns.plot, paste0(unknown.name, "traceREE.csv"))




# discrimination plots

# Th vs U

sc <- scale_fill_viridis(option = 'plasma', direction = -1, limits=c(500,2500))
geompt <- geom_point(aes(fill=age), colour = "black", pch = 21, size = 3.5, stroke = 0.8)

unknowns.plot %>%
  ggplot(aes(x = Th_ppm_m232, y = U_ppm_m238)) +
  geompt +
  scale_y_log10(labels = comma, limits= c(0.01, 1000)) +
  scale_x_log10(labels = comma, limits= c(0.01, 1000)) +
  xlab("Th (ppm)") +
  ylab("U (ppm)") +
  labs(fill="Age (Ma)") +
  geom_abline(intercept = 0 , slope = 1) +
  geom_abline(intercept = 1 , slope = 1) +
  geom_abline(intercept = -1 , slope = 1) +
  sc

ggsave(paste0(unknown.name, "Th vs U.pdf"), width = 7, height = 7)



# Mn vs Sr

unknowns.plot %>%
  ggplot(aes(x = Mn_ppm_m55, y = Sr_ppm_m88)) +
  geompt +
  scale_y_log10(labels = comma,limits= c(10, 10000)) +
  scale_x_log10(labels = comma,limits= c(10, 10000)) +
  xlab("Mn (ppm)") +
  ylab("Sr (ppm)") +
  labs(fill="Age (Ma)") +
  sc

ggsave(paste0(unknown.name, "Mn vs Sr.pdf"), width = 7, height = 7)


# Y vs Sr

unknowns.plot %>%
  ggplot(aes(x = Y_ppm_m89, y = Sr_ppm_m88)) +
  geompt +
  scale_y_log10(labels = comma,limits= c(10, 10000)) +
  scale_x_log10(labels = comma,limits= c(10, 10000)) +
  xlab("Y (ppm)") +
  ylab("Sr (ppm)") +
  labs(fill="Age (Ma)") +
  sc

ggsave(paste0(unknown.name, "Y vs Sr.pdf"), width = 7, height = 7)


# Eu vs Y

unknowns.plot %>%
  ggplot(aes(x = Eu.Eua, y = Y_ppm_m89)) +
  geompt +
  scale_y_log10(labels = comma,limits= c(9, 10000)) +
  scale_x_log10(labels = comma,limits= c(0.01, 10)) +
  xlab("Eu/Eu*") +
  ylab("Y (ppm)") +
  labs(fill="Age (Ma)") +
  sc

ggsave(paste0(unknown.name, "Eu vs Y.pdf"), width = 7, height = 7)


# sumREE vs Ce/Yb

unknowns.plot %>%
  ggplot(aes(x = sumREE, y = Ce.Yb)) +
  geompt +
  scale_y_log10(labels = comma,limits= c(0.01, 10000)) +
  scale_x_continuous(limits= c(0, 3)) +
  xlab("Total REE (wt%)") +
  ylab("(Ce/Yb)cn") +
  labs(fill="Age (Ma)") +
  sc

ggsave(paste0(unknown.name, "ree vs ceyb.pdf"), width = 7, height = 7)


# lights/sumREE vs La/Nd

unknowns.plot %>%
  ggplot(aes(x = lights.sumREE, y = La.Nd)) +
  geompt +
  scale_y_continuous(limits= c(0, 6)) +
  scale_x_continuous(limits= c(0, 100)) +
  xlab("(La+Ce)/REE(tot)") +
  ylab("La/Nd") +
  labs(fill="Age (Ma)") +
  sc

ggsave(paste0(unknown.name, "l.ree vs land.pdf"), width = 7, height = 7)

