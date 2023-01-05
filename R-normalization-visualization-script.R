if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")
BiocManager::install("ape")
install.packages('SRS') #done
install.packages("vegan") #done
install.packages('ggplot2')
install.packages('ggstar')
BiocManager::install("treeio")
install.packages("reshape2")
install.packages("GGally")
install.packages("devtools")


library('SRS')
library('vegan')
library('ggplot2')
library('ggtree')
library('ggstar')
library('treeio')
library('ape')
library('reshape2')
library("GGally") 
library("devtools")
library("RColorBrewer")
library('dplyr')
library('tidyr')

n <- 32
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

family <- read.csv("minimap2_unite_family.txt",sep = "\t",header=TRUE )
head(family)
species<-read.csv("minimap_spp_species.txt",sep = "\t",header=TRUE )
head(species)

familyU <- family[c(2:9)] #This select columns 2-4 of table genusT and copies them to genusT2
head(familyU)


speciesU <- species[c(2:9)] #This select columns 2-4 of table genusT and copies them to genusT2
head(speciesU)


CminF<-min(colSums(familyU)) #sum up each individual column in genusT2 table and then store lowest value in CminG
CminF

CminS<-min(colSums(speciesU)) #sum up each individual column in genusT2 table and then store lowest value in CminG
CminS

species_df <- as.data.frame(t(speciesU)) #convert table to data.frame type object with t indicating switch of columns and rows
colnames(species_df)<-species[,1] # Add the genus names which where 1st column of genusT back in as column names in Genus_df table
head(species_df)

family_df <- as.data.frame(t(familyU)) #convert table to data.frame type object with t indicating switch of columns and rows
colnames(family_df)<-family[,1] # Add the genus names which where 1st column of genusT back in as column names in Genus_df table
head(family_df)

rarecurve(species_df, step = 20, sample = CminS, col = "blue", cex = 0.6) #Plot species vs sample size after rarefy
rarecurve(family_df, step = 20, sample = CminF, col = "blue", cex = 0.6) #Plot species vs sample size after rarefy

S <- specnumber(species_df) # observed number of species
raremax <- min(rowSums(species_df)) #calculate the smallest number of reads at any one of the sites 
Srare <- rarefy(species_df, raremax) #carry out rarefaction 
plot(S, Srare, xlab = "Observed No. of species", ylab = "Rarefied No. of species", ) # Create plot 
text(Srare~S, labels=names(S), font=2, pos=1) #label graph
abline(lm(Srare ~ S)) # add line to graph

F <- specnumber(family_df) # observed number of species
raremaxF <- min(rowSums(family_df)) #calculate the smallest number of reads at any one of the sites 
Frare <- rarefy(family_df, raremaxF) #carry out rarefaction 
plot(F, Frare, xlab = "Observed No. of families", ylab = "Rarefied No. of families", ) # Create plot 
text(Frare~F, labels=names(F), font=2, pos=1) #label graph
abline(lm(Frare ~ F)) # add line to graph

SRS_outputF2 <- SRS(familyU, CminF, set_seed=TRUE) #Carry out Scaling with Ranked Subsampling (SRS)
row.names(SRS_outputF2)<-family[,1] #Add genus names back to table 
SRS_outputF2

SRS_outputS2 <- SRS(speciesU, CminS, set_seed=TRUE) #Carry out Scaling with Ranked Subsampling (SRS)
row.names(SRS_outputS2)<-species[,1] #Add genus names back to table 
SRS_outputS2

write.csv(SRS_outputF2,"SRS_outputF-12-16.csv", row.names = TRUE)
write.csv(SRS_outputS2,"SRS_outputS-12-16.csv", row.names = TRUE)

SRSplotF<-SRScurve(familyU, metric = "richness", step = 500,
                  rarefy.comparison = TRUE ,ylab = "richness", rarefy.comparison.legend = TRUE,
                  col = c("red", "red","blue","blue","black","black"),
                  lty = c(1,2)) # From here to the top is the SRS function to create comparison with rarefy, col argument = line colour, lty = line type , step is sampling times 
segments(x0 = CminF, y0 = -55, x1 = CminF, y1 = 700,) #Add line at cMinG value colour the line dark green 


SRSplotS<-SRScurve(speciesU, metric = "richness", step = 500,
                  rarefy.comparison = TRUE ,ylab = "richness", rarefy.comparison.legend = TRUE,
                  ) # From here to the top is the SRS function to create comparison with rarefy, col argument = line colour, lty = line type , step is sampling times 
segments(x0 = CminS, y0 = -55, x1 = CminS, y1 = 700,)


## DIVERSITY METRICS ###

print("SRS Family richness")
Family_SRS <- as.data.frame(t(SRS_outputF2))  
specnumber(Family_SRS)

Species_SRS <- as.data.frame(t(SRS_outputS2))  
View(Species_SRS)
specnumber(Species_SRS)

shan <- diversity(Family_SRS, index = "shannon") #using vegan functions calculate diversity index
print(shan)
shan2<-as.data.frame(shan) #convert results from list to data frame 
xlab<-(rownames(shan2)) #label the x-axis 
print(shan2)
q<-ggplot(data=shan2, aes(x=xlab, y=shan))+geom_bar(stat="identity") # Make plot 
q+ ggtitle("Plot of shannon diversity at each site") +
  xlab("Site") + ylab("shannon diversity") #Add titles to plot 
ggsave("shan-div-indiv-family.pdf", plot = q)

Specshan <- diversity(Species_SRS, index = "shannon") #using vegan functions calculate diversity index
print(Specshan)
Specshan2<-as.data.frame(Specshan) #convert results from list to data frame 
xlab<-(rownames(Specshan2)) #label the x-axis 
print(Specshan2)
q<-ggplot(data=Specshan2, aes(x=xlab, y=Specshan))+geom_bar(stat="identity") # Make plot 
q+ ggtitle("Plot of shannon diversity at each site") +
  xlab("Site") + ylab("shannon diversity") #Add titles to plot 
ggsave("shan-div-unite-species.pdf", plot = q)

simpF <- diversity(Family_SRS, index = "simpson")
print(simpF)
simpF2<-as.data.frame(simpF)
xlab<-(rownames(simpF2))
p<-ggplot(data=simpF2, aes(x=xlab, y=simpF))+geom_bar(stat="identity")
p + ggtitle("Plot of simpson diversity at each site") +
  xlab("Site") + ylab(" Simpson diversity")
ggsave("simp-div-unite-fam.pdf", plot = p)

simpS <- diversity(Species_SRS, index = "simpson")
print(simpS)
simpS2<-as.data.frame(simpS)
xlab<-(rownames(simpS2))
p<-ggplot(data=simpS2, aes(x=xlab, y=simpS))+geom_bar(stat="identity")
p + ggtitle("Plot of simpson diversity at each site") +
  xlab("Site") + ylab(" Simpson diversity")
ggsave("simp-div-unite-species.pdf", plot = p)

invsimpF <- diversity(Family_SRS, index = "invsimpson")
print(invsimpF)
invsimpF2<-as.data.frame(invsimpF)
xlab<-(rownames(invsimpF2))
p<-ggplot(data=invsimpF2, aes(x=xlab, y=invsimpF))+geom_bar(stat="identity")
p + ggtitle("Plot of inverse simpson diversity at each site") +
  xlab("Site") + ylab("Inverse simpson diversity")
ggsave("Invsimp-div-indiv-fam.pdf", plot = p)

invsimpS <- diversity(Species_SRS, index = "invsimpson")
print(invsimpS)
invsimpS2<-as.data.frame(invsimpS)
xlab<-(rownames(invsimpS2))
p<-ggplot(data=invsimpS2, aes(x=xlab, y=invsimpS))+geom_bar(stat="identity")
p + ggtitle("Plot of inverse simpson diversity at each site") +
  xlab("Site") + ylab("Inverse simpson diversity")
ggsave("Invsimp-div-unite-species.pdf", plot = p)

unbias.simp <- simpson.unb(Family_SRS)
head(unbias.simp)
## Fisher alpha
alpha <- fisher.alpha(Family_SRS)
group2<-vector()
group2[row.names(Family_SRS)[1]] <- 1
group2[row.names(Family_SRS)[2]] <- 2
group2[row.names(Family_SRS)[3]] <- 3
group2[row.names(Family_SRS)[4]] <- 4
group2[row.names(Family_SRS)[5]] <- 5
group2[row.names(Family_SRS)[6]] <- 6
group2[row.names(Family_SRS)[7]] <- 7
group2[row.names(Family_SRS)[8]] <- 8

diversity<-cbind(shan, simpF, invsimpF, unbias.simp, alpha)

pairs(cbind(shan, simpF, invsimpF, unbias.simp, alpha), pch=c(1,1,1,1,8,8,5,5), col = hcl.colors(l, "viridis")[group2],oma=c(3,3,3,10))+
legend("bottomright", fill = hcl.colors(l, "viridis")[group2],legend = c( group2))

diversity2<-diversity
diversity$names<-rownames(diversity)
View(diversity2)
ggpairs(diversity, columns = 1:5, ggplot2::aes(colour=names)) #NEED TO WORK ON TWO PROBLEMS 1 not showing all daata? 2. Colour mapping


print("Pielous evenness Family")
SRF <- specnumber(Family_SRS)
shan/log(SRF)
print("Pielous evenness Species")
SRspec <- specnumber(Species_SRS)
Specshan/log(SRspec)


Family2.pca <- prcomp(SRS_outputF2, center = TRUE,scale. = FALSE)
Family1.pca <- prcomp(Family_SRS, center = TRUE,scale. = FALSE)
autoplot(Family1.pca, label=TRUE)
biplot(Family1.pca, label=TRUE)

pcs = data.frame(Family1.pca$x)
p = pcs %>% 
  ggplot(aes(x = PC1, 
             y = PC2)) +
  geom_point()
p
p+geom_text_repel(aes(label=row.names(Family_SRS)),size=3)

Species2.pca <- prcomp(SRS_outputS2, center = TRUE,scale. = FALSE)
Species1.pca <- prcomp(Species_SRS, center = TRUE,scale. = FALSE)
autoplot(Species1.pca, label=TRUE)
biplot(Species1.pca, label=TRUE)

pcsS = data.frame(Species1.pca$x)
pS = pcsS %>% 
  ggplot(aes(x = PC1, 
             y = PC2)) +
  geom_point()
pS
pS+geom_text_repel(aes(label=row.names(Species_SRS)),size=3)


ord <- princomp(Factor_SRS[, 1:4])


FamilyPlot<-as.data.frame(Family_SRS)
FamilyPlot$soiltype<-rownames(Family_SRS)
FamilyPlot2<-cSplit(FamilyPlot, "soiltype", ".")
FamilyPlot3 <- FamilyPlot2 %>% pivot_longer(cols=soiltype_1)
FamilyPlot4 = subset(FamilyPlot3, select = -c(soiltype_2,soiltype_3,soiltype_4,name) )
fit <- prcomp(FamilyPlot4[,c(1:97)])
autoplot( fit, data=FamilyPlot4, colour="value", loadings=TRUE, loadings.label=TRUE,loadings.label.size=2, loadings.label.repel=TRUE )

SpeciesPlot<-as.data.frame(Species_SRS)
SpeciesPlot$soiltype<-rownames(Species_SRS)
SpeciesPlot2<-cSplit(SpeciesPlot, "soiltype", ".")
SpeciesPlot3 <-SpeciesPlot2 %>% pivot_longer(cols=soiltype_1)
SpeciesPlot4 = subset(SpeciesPlot3, select = -c(soiltype_2,soiltype_3,soiltype_4,name) )
fitS <- prcomp(SpeciesPlot4[,c(1:208)])
autoplot( fitS, data=SpeciesPlot4, colour="value", loadings=TRUE, loadings.label=TRUE,loadings.label.size=2, loadings.label.repel=TRUE )

var_explained_df <- data.frame(PC= paste0("PC",1:8), var_explained=(fitS$sdev)^2/sum((fitS$sdev)^2))
ggplot(aes(x=PC, y=var_explained))+
geom_point(size=4)+
geom_line()+
labs(title="Scree plot: PCA on scaled data")

var_explained_dfF <- data.frame(PC= paste0("PC",1:8), var_explained=(fit$sdev)^2/sum((fit$sdev)^2))
var_explained_dfF %>%
ggplot(aes(x=PC, y=var_explained))+
  geom_point(size=4)+
  geom_line()+
  labs(title="Scree plot Family: PCA on scaled data")

AbundF2<-SRS_outputF2
iF1 <- order(decreasing = T,rowSums(AbundF2[c("SM1.R1","IV3","BM2.R2","BM2.D3.P8","BM2.D1.P7","BM2.D1.P7.1","SM1.R2","IV2")]))
AbundF2new <- AbundF2[iF1, ]
head(AbundF2new)
AbundF2new$X <- rownames(AbundF2new)

row.names.remove <- c("unidentified","Saxifragales")
AbundF2new<-AbundF2new[!(row.names(AbundF2new) %in% row.names.remove), ]



AbundF4<- head(AbundF2new, n=10)
AbundF4$X<- rownames(AbundF4)
meltedF <- melt(AbundF4, id.vars = "X")
meltedF$ranks <- factor(meltedF$X,levels=unique(meltedF$X))

ggplot(meltedF, aes(x = variable, y = value, fill = X)) +
geom_bar(stat = "identity")+scale_fill_manual(values=col_vector)



##PERCENTAGE ### 
AbundF2newp<-AbundF2new[!(row.names(AbundF2new) %in% row.names.remove), ]

#Convert percentage to loop 
AbundF2newp$SM1.R1= 100*(AbundF2newp$SM1.R1/sum(AbundF2newp$SM1.R1))
AbundF2newp$IV3= 100*(AbundF2newp$IV3/sum(AbundF2newp$IV3))
AbundF2newp$BM2.R2= 100*(AbundF2newp$BM2.R2/sum(AbundF2newp$BM2.R2))
AbundF2newp$BM2.D3.P8= 100*(AbundF2newp$BM2.D3.P8/sum(AbundF2newp$BM2.D3.P8))
AbundF2newp$BM2.D1.P7= 100*(AbundF2newp$BM2.D1.P7/sum(AbundF2newp$BM2.D1.P7))
AbundF2newp$BM2.D1.P7.1= 100*(AbundF2newp$BM2.D1.P7.1/sum(AbundF2newp$BM2.D1.P7.1))
AbundF2newp$SM1.R2= 100*(AbundF2newp$SM1.R2/sum(AbundF2newp$SM1.R2))
AbundF2newp$IV2= 100*(AbundF2newp$IV2/sum(AbundF2newp$IV2))
AbundP<-AbundF2newp
AbundP$X<- rownames(AbundF2newp)
meltedFP <- melt(AbundP, id.vars = "X")
meltedFP$ranks <- factor(meltedFP$X,levels=unique(meltedFP$X))
ggplot(meltedFP, aes(x = variable, y = value, fill = X)) +
geom_bar(stat = "identity")

#filter atleast one at 1 percent 
dat.query <- filter(AbundP, SM1.R1 >= 1 |IV3 >= 1 | BM2.R2 >= 1 |BM2.D3.P8>= 1 | BM2.D1.P7.1>= 1 |BM2.D1.P7 >= 1|SM1.R2 >= 1 |IV2 >= 1 )
meltedFP1 <- melt(dat.query, id.vars = c("X"))
ggplot(meltedFP1, aes(x = variable, y = value, fill = X)) +
geom_bar(stat = "identity")

write.csv(AbundP,"SRS_outputP-10-3-family-unite-percent-saxifragelles-unidentified-rm.csv", row.names = TRUE)


##SPECIES BARPLOTS ### 

AbundS2<-SRS_outputS2
iS1 <- order(decreasing = T,rowSums(AbundS2[c("SM1.R1","IV3","BM2.R2","BM2.D3.P8","BM2.D1.P7","BM2.D1.P7.1","SM1.R2","IV2")]))
AbundS2new <- AbundS2[iS1, ]
head(AbundS2new)
AbundS2new$X <- rownames(AbundS2new)

row.names.removeS <- c("unidentified","Mitella_pauciflora")
AbundS2new<-AbundS2new[!(row.names(AbundS2new) %in% row.names.removeS), ]



AbundS4<- head(AbundS2new, n=10)
AbundS4$X<- rownames(AbundS4)
meltedS <- melt(AbundS4, id.vars = "X")
meltedS$ranks <- factor(meltedS$X,levels=unique(meltedS$X))

ggplot(meltedS, aes(x = variable, y = value, fill = X)) +
  geom_bar(stat = "identity")



##PERCENTAGE ### 
AbundS2newp<-AbundS2new[!(row.names(AbundS2new) %in% row.names.remove), ]

#Convert percentage to loop 
AbundS2newp$SM1.R1= 100*(AbundS2newp$SM1.R1/sum(AbundS2newp$SM1.R1))
AbundS2newp$IV3= 100*(AbundS2newp$IV3/sum(AbundS2newp$IV3))
AbundS2newp$BM2.R2= 100*(AbundS2newp$BM2.R2/sum(AbundS2newp$BM2.R2))
AbundS2newp$BM2.D3.P8= 100*(AbundS2newp$BM2.D3.P8/sum(AbundS2newp$BM2.D3.P8))
AbundS2newp$BM2.D1.P7= 100*(AbundS2newp$BM2.D1.P7/sum(AbundS2newp$BM2.D1.P7))
AbundS2newp$BM2.D1.P7.1= 100*(AbundS2newp$BM2.D1.P7.1/sum(AbundS2newp$BM2.D1.P7.1))
AbundS2newp$SM1.R2= 100*(AbundS2newp$SM1.R2/sum(AbundS2newp$SM1.R2))
AbundS2newp$IV2= 100*(AbundS2newp$IV2/sum(AbundS2newp$IV2))
AbundSP<-AbundS2newp
AbundSP$X<- rownames(AbundS2newp)
meltedSP <- melt(AbundSP, id.vars = "X")
meltedSP$ranks <- factor(meltedSP$X,levels=unique(meltedSP$X))
ggplot(meltedSP, aes(x = variable, y = value, fill = X)) +
geom_bar(stat = "identity")

write.csv(AbundSP,"SRS_outputP-10-3-species-unite-percent-unidentified-rm.csv", row.names = TRUE)


#filter atleast one at 1 percent 
dat.queryS <- filter(AbundSP, SM1.R1 >= 1 |IV3 >= 1 | BM2.R2 >= 1 |BM2.D3.P8>= 1 | BM2.D1.P7.1>= 1 |BM2.D1.P7 >= 1|SM1.R2 >= 1 |IV2 >= 1 )
meltedSP1 <- melt(dat.queryS, id.vars = c("X"))
ggplot(meltedSP1, aes(x = variable, y = value, fill = X)) +
geom_bar(stat = "identity")

dat.queryMycoP <- AbundSP %>% filter(grepl('Glomus|Septoglomus|Rhizophagus|Sclerocystis|Claroideoglomus|Pacispora|Racocetra|Gigaspora|Acaulospora|Redeckera|Diversispora|Geosiphon|Ambispora|Archaeospora|Paraglomus|Scutellospora_|Funneliformis_',X))

dat.queryMyco <- AbundS2new%>% filter(grepl('Glomus|Septoglomus|Rhizophagus|Sclerocystis|Claroideoglomus|Pacispora|Racocetra|Gigaspora|Acaulospora|Redeckera|Diversispora|Geosiphon|Ambispora|Archaeospora|Paraglomus|Scutellospora_|Funneliformis_|Glomeromyc|Dominikia|Microdominikia|Kamienskia|Rhizoglomus|Nanoglomus|Oehlia|Microkamienskia|Funneliglomus|Sclerocarpum|Orientoglomus',X))
View(dat.queryMyco)
meltedSPM <- melt(dat.queryMyco, id.vars = c("X"))
bubbleplot <- ggplot(meltedSPM, aes(x = variable, y = ranks)) + 
  geom_point(aes(size = value, fill = ranks), alpha = 0.75, shape = 21) + 
  scale_size_continuous(range = c(1,20),breaks=c(50,100,200)) + 
  labs( x= "Site", y = "Species", size = "Reads per Genera", fill = "Species")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right") #+   
bubbleplot


query.P.lit<-AbundSP %>% filter(grepl('Yarrowia|Alternaria|Arthrobotrys|Aspergillus|Cephalosporium|Cladosporium|Curvularia|Cunninghamella|Chaetomium|Fusarium|Glomus|Helminthosporium|Micromonospora|Mortierella|Myrothecium|Oidiodendron|Paecilomyces|Penicillium|Phoma|Pichia|Populospora|Pythium|Rhizoctonia|Rhizopus|Saccharomyces|Schizosaccharomyces|Schwanniomyces|Sclerotium|Torula|Trichoderma',X))
meltedPlit <- melt(query.P.lit, id.vars = c("X"))
bubbleplot <- ggplot(meltedPlit, aes(x = variable, y = X)) + 
  geom_point(aes(size = value, fill = X), alpha = 0.75, shape = 21) + 
  scale_size_continuous(range = c(1,20),breaks=c(50,100,200)) + 
  labs( x= "Site", y = "Genera", size = "Reads per Genera", fill = "Genera")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right") #+   
bubbleplot


meltedSP1 <- melt(dat.queryS, id.vars = c("X"))