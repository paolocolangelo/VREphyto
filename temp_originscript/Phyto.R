phyto<-read.csv("PhytoplanktonWiserProject_2009.csv",sep=";")


### Indici Phyto ###


##### sinking #####

#Bacillariphyceae
Bacillariphyceae<-subset(phyto,class=="Bacillariophyceae")
sinking_Bacillariophyceae=0.036*(Bacillariphyceae$biovolume)^0.33
#hist(sinking_Bacillariophyceae,breaks = 30)

#Dinophyceae (Miozoa)
Dinophyceae<-subset(phyto,class=="Dinophyceae")
sinking_Dinophyceae=0.036*(Dinophyceae$biovolume)^0.33
#hist(sinking_Dinophyceae,breaks = 30)



####  Nutrient uptake ####

# Vmax = µmol nutrient cell-1day-1
# K= µmol nutrient L-1

Diatom<-subset(phyto,class=="Bacillariophyceae")
Vcell<-Diatom$biovolume
#Marine Diatom N(Nitrate)

log10VmaxN=-7.8+0.67*log(Vcell)
log10KN=-0.49+0.17*log(Vcell)

#Marine Diatom P
log10VmaxP=-7.95+0.81*log(Vcell)
log10KP=-1.16+0.39*log(Vcell)

#Uptake affinity

# Vmax/K
Uptake_affinityN=log10VmaxN/log10KN
Uptake_affinityP=log10VmaxP/log10KP

#Ligth absorpion

# log_a*(ʎ)=b+clog(Vcell)

## 440nm
b=-1.66
c=-0.067
log_a440=b+c*log(Vcell)


# 675nm
b=-1.82
c=-0.052
log_a675=b+c*log(Vcell)


## summary table

Diatom_summary_table<-cbind(Diatom[1:24],log10VmaxN,log10KN,log10VmaxP,log10KP,log_a440,log_a675)

par(mfrow=c(4,2))
hist(log(Diatom_summary_table$biovolume))
hist(log(Diatom_summary_table$cellcarboncontent))
hist(Diatom_summary_table$log10VmaxN)
hist(Diatom_summary_table$log10KN)
hist(Diatom_summary_table$log10VmaxP)
hist(Diatom_summary_table$log10KP)
hist(Diatom_summary_table$log_a440)
hist(Diatom_summary_table$log_a675)

# alternative with lattice
#count
lattice::histogram(~log(biovolume)+log(cellcarboncontent)+log10VmaxN+log10KN+log10VmaxP+log10KP+log_a440+log_a675,data=Diatom_summary_table,type="count")
#percent
lattice::histogram(~log(biovolume)+log(cellcarboncontent)+log10VmaxN+log10KN+log10VmaxP+log10KP+log_a440+log_a675,data=Diatom_summary_table)
#lattice::densityplot(~log(biovolume)+log(cellcarboncontent)+log10VmaxN+log10KN+log10VmaxP+log10KP+log_a440+log_a675,data=Diatom_summary_table)
#lattice::densityplot(~log(biovolume)+log(cellcarboncontent),data=Diatom_summary_table)
lattice::xyplot(log(biovolume)~log(cellcarboncontent)+log10VmaxN,data=Diatom_summary_table)
#### Community analysis ####

### Rarefaction curves
library(vegan)
tabella_specie_abundance<-table(phyto$parenteventid, phyto$scientificname)
rarecurve(tabella_specie_abundance)

### co-occurence with cooccur library
library(cooccur)
tabella_specie<-tabella_specie_abundance
tabella_specie[tabella_specie>1]<-1
cooccur.phyto <- cooccur(mat = t(tabella_specie), thresh = TRUE, spp_names = T)
summary(cooccur.phyto)
prob.table(cooccur.phyto)
#plot(cooccur.phyto)
pair.profile(cooccur.phyto)
obs.v.exp(cooccur.phyto)

### co-occurence with EcoSimR library
library(EcoSimR)
tabella_specie<-tabella_specie_abundance
tabella_specie[tabella_specie>1]<-1
fm<-cooc_null_model(tabella_specie, algo="sim9",nReps=10000,burn_in = 1000)
summary(fm)
plot(fm,type="burn_in")
plot(fm,type="hist")
