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
install.packages("data.table")  


library("data.table")
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
library('ggrepel')
library('splitstackshape')
library('tidyr')
library(devtools)
install_github('sinhrks/ggfortify')
library(ggfortify); library(ggplot2)

n <- 32
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

genus <- read.csv("minimap2_k12_genus.txt",sep = "\t",header=TRUE )
head(genus)


genusU <- genus[c(2:5)] #This select columns 2-4 of table genusT and copies them to genusT2
head(genusU)

CminG<-min(colSums(genusU)) #sum up each individual column in genusT2 table and then store lowest value in CminG
CminG

genus_df <- as.data.frame(t(genusU)) #convert table to data.frame type object with t indicating switch of columns and rows
colnames(genus_df)<-genus[,1] # Add the genus names which where 1st column of genusT back in as column names in Genus_df table
head(genus_df)

rarecurve(genus_df, step = 500, sample = CminG, col = "blue", cex = 0.6) #Plot species vs sample size after rarefy

g <- specnumber(genus_df) # observed number of species
raremax <- min(rowSums(genus_df)) #calculate the smallest number of reads at any one of the sites 
Grare <- rarefy(genus_df, raremax) #carry out rarefaction 
plot(g, Grare, xlab = "Observed No. of genera", ylab = "Rarefied No. of genera", ) # Create plot 
text(Grare~g, labels=names(g), font=2, pos=1) #label graph
abline(lm(Grare ~ g)) # add line to graph


SRS_outputG2 <- SRS(genusU, CminG, set_seed=TRUE) #Carry out Scaling with Ranked Subsampling (SRS)
row.names(SRS_outputG2)<-genus[,1] #Add genus names back to table 
SRS_outputG2


write.csv(SRS_outputG2,"SRS_outputG-12-27.csv", row.names = TRUE)

SRSplotG<-SRScurve(genusU, metric = "richness", step = 500,
                  rarefy.comparison = TRUE ,ylab = "richness", rarefy.comparison.legend = TRUE,
                  col = c("red", "red","blue","blue","black","black"),
                  lty = c(1,2)) # From here to the top is the SRS function to create comparison with rarefy, col argument = line colour, lty = line type , step is sampling times 
segments(x0 = CminF, y0 = -55, x1 = CminF, y1 = 700,) #Add line at cMinG value colour the line dark green 





## DIVERSITY METRICS ###

print("SRS Genus richness")
genus_SRS <- as.data.frame(t(SRS_outputG2))  
specnumber(genus_SRS)


shan <- diversity(genus_SRS, index = "shannon") #using vegan functions calculate diversity index
print(shan)
shan2<-as.data.frame(shan) #convert results from list to data frame 
xlab<-(rownames(shan2)) #label the x-axis 
print(shan2)
q<-ggplot(data=shan2, aes(x=xlab, y=shan))+geom_bar(stat="identity") # Make plot 
q+ ggtitle("Plot of shannon diversity at each site") +
  xlab("Site") + ylab("shannon diversity") #Add titles to plot 
ggsave("shan-16s-genus.pdf", plot = q)

simpG <- diversity(genus_SRS, index = "simpson")
print(simpG)
simpG2<-as.data.frame(simpG)
xlab<-(rownames(simpG2))
p<-ggplot(data=simpG2, aes(x=xlab, y=simpG))+geom_bar(stat="identity")
p + ggtitle("Plot of simpson diversity at each site") +
  xlab("Site") + ylab(" Simpson diversity")
ggsave("simp-div-16s-genus.pdf", plot = p)


invsimpG <- diversity(genus_SRS, index = "invsimpson")
print(invsimpG)
invsimpG2<-as.data.frame(invsimpG)
xlab<-(rownames(invsimpG2))
p<-ggplot(data=invsimpG2, aes(x=xlab, y=invsimpG))+geom_bar(stat="identity")
p + ggtitle("Plot of inverse simpson diversity at each site") +
  xlab("Site") + ylab("Inverse simpson diversity")
ggsave("Invsimp-div-indiv-genus.pdf", plot = p)



unbias.simp <- simpson.unb(genus_SRS)
head(unbias.simp)
## Fisher alpha
alpha <- fisher.alpha(genus_SRS)
group2<-vector()
group2[row.names(genus_SRS)[1]] <- 1
group2[row.names(genus_SRS)[2]] <- 2
group2[row.names(genus_SRS)[3]] <- 3
group2[row.names(genus_SRS)[4]] <- 4

diversity<-as.dataframe(cbind(shan, simpG, invsimpG, unbias.simp, alpha))

pairs(cbind(shan, simpG, invsimpG, unbias.simp, alpha), pch=c(1,4,5,3,6),oma=c(3,3,3,10))+
legend("bottomright",legend = c(group2))

diversity2<-diversity
diversity$names<-rownames(diversity)
View(diversity2)
ggpairs(diversity, columns = 1:5, ggplot2::aes(colour=names)) #NEED TO WORK ON TWO PROBLEMS 1 not showing all daata? 2. Colour mapping


print("Genus evenness Family")
GRF <- specnumber(genus_SRS)
shan/log(GRF)



genus2.pca <- prcomp(SRS_outputG2, center = TRUE,scale. = FALSE)
genus1.pca <- prcomp(genus_SRS, center = TRUE,scale. = FALSE)
autoplot(genus1.pca, label=TRUE)
autoplot(genus2.pca, label=TRUE)
biplot(genus1.pca, label=TRUE)

pcs = data.frame(genus1.pca$x)
p = pcs %>% 
  ggplot(aes(x = PC1, 
             y = PC2)) +
  geom_point()
p
p+geom_text_repel(aes(label=row.names(genus_SRS)),size=3)


ord <- princomp(genus_SRS[, 1:4])


genusPlot<-as.data.frame(genus_SRS)
fit <- prcomp(genusPlot)
autoplot( fit, data=genusPlot, loadings=TRUE, ntopLoadings=3, loadings.label=TRUE,loadings.label.size=2, loadings.label.repel=TRUE )

***DOESNT WORK YET 

fit <- svd(genus_SRS)
fit <- factanal(genusPlot4[,c(1:782)],3,rotation = 'varimax')
print(fit, digits=2, cutoff=.3, sort=TRUE)

var_explained_df <- data.frame(PC= paste0("PC",1:8), var_explained=(fitS$sdev)^2/sum((fitS$sdev)^2))
var_explained_df %>%
ggplot(aes(x=PC, y=var_explained))+
geom_point(size=4)+
geom_line()+
labs(title="Scree plot: PCA on scaled data")


pca1 <- rda(SRS_outputG2, scale = TRUE)
biplot(pca1, scaling=3, type = c("text", "text"))
autoplot( fit, data=genusPlot, loadings=TRUE, ntopLoadings=3, loadings.label=TRUE,loadings.label.size=2, loadings.label.repel=TRUE)


****

AbundG2<-SRS_outputG2
iG1 <- order(decreasing = T,rowSums(AbundG2[c("BM2.R2","SM1.R2","SM1.R1","BM2.R1")]))
AbundG2new <- AbundG2[iG1, ]
head(AbundG2new)
AbundG2new$X <- rownames(AbundG2new)

row.names.remove <- c("uncultured")
AbundG2new<-AbundG2new[!(row.names(AbundG2new) %in% row.names.remove), ]




AbundG4<- head(AbundG2new, n=10)
AbundG4$X<- rownames(AbundG4)
meltedG <- melt(AbundG4, id.vars = "X")
meltedG$ranks <- factor(meltedG$X,levels=unique(meltedG$X))

ggplot(meltedG, aes(x = variable, y = value, fill = X)) +
geom_bar(stat = "identity")+scale_fill_manual(values=col_vector)



##PERCENTAGE ### 
AbundG2newp<-AbundG2new[!(row.names(AbundG2new) %in% row.names.remove), ]

#Convert percentage to loop 
AbundG2newp$BM2.R2= 100*(AbundG2newp$BM2.R2/sum(AbundG2newp$BM2.R2))
AbundG2newp$SM1.R2= 100*(AbundG2newp$SM1.R2/sum(AbundG2newp$SM1.R2))
AbundG2newp$SM1.R1= 100*(AbundG2newp$SM1.R1/sum(AbundG2newp$SM1.R1))
AbundG2newp$BM2.R1= 100*(AbundG2newp$BM2.R1/sum(AbundG2newp$BM2.R1))

AbundP<-AbundG2newp
AbundP$X<- rownames(AbundG2newp)
meltedGP <- melt(AbundP, id.vars = "X")
meltedGP$ranks <- factor(meltedGP$X,levels=unique(meltedGP$X))
ggplot(meltedGP, aes(x = variable, y = value, fill = X)) +
geom_bar(stat = "identity")

#filter atleast one at 1 percent 
dat.query <- filter(AbundP, BM2.R2 >= 1 |SM1.R2 >= 1 | SM1.R1 >= 1 |BM2.R1 >= 1 )
dat.querys <- filter(AbundP, BM2.R2 >= 0.1 |SM1.R2 >= 0.1 | SM1.R1 >= 0.1 |BM2.R1 >= 0.1 )
meltedGP1 <- melt(dat.query, id.vars = c("X"))
meltedGP1s <- melt(dat.querys, id.vars = c("X"))
ggplot(meltedGP1, aes(x = variable, y = value, fill = X)) +
  geom_bar(stat = "identity", linewidth= 1, colour = "black")

ggplot(meltedGP1s, aes(x = variable, y = value, fill = X)) +
  geom_bar(stat = "identity", linewidth= 0.1, colour = "black")+
  scale_colour_viridis(D)

write.csv(AbundP,"SRS_outputP-10-3-genus-SILVA-1percent-uncultured-rm.csv", row.names = TRUE)

bubbleplot <- ggplot(meltedGP1, aes(x = variable, y = X)) +
geom_point(aes(size = value, fill = X),shape = 21) +
scale_size_continuous(range = c(1,10)) +
labs( x= "Site", y = "Genera", size = "Percent of reads assigned each Genera > 1%", fill = "X")
bubbleplot

View(AbundP)

dat.queryCollected <- filter(AbundP, X == 'Bacillus' |X == 'Rahnella1' | X =='Bradyrhizobium' | X=='Pseudomonas' | X=='Ensifer' |X=='Azotobacter' | X == 'Micromonospora' |X == 'Frankia' |X == 'Aspergillus' | X =='Azospirillum' | X=='Streptomyces' | X=='Gluconacetobacter' |X=='Hyphomicrobium')
meltedCollected <- melt(dat.queryCollected, id.vars = c("X"))
meltedCollected[meltedCollected == 0]<-NA
bubbleplot <- ggplot(meltedCollected, aes(x = variable, y = X)) +
  geom_point(aes(size = value)) +
  scale_size_continuous(range = c(1,10)) +
  labs( x= "Site", y = "Genera", size = "Percent of reads assigned each Genera", fill = "X")
  
bubbleplot



dat.queryRhizobium <- AbundP[AbundP$X %like% "rhizobium", ]        # Extract matching rows with %like%
head(dat.queryRhizobium )
meltedRhizo <- melt(dat.queryRhizobium, id.vars = c("X"))
meltedRhizo[meltedRhizo == 0]<-NA
bubbleplot <- ggplot(meltedRhizo, aes(x = variable, y = X)) +
  geom_point(aes(size = value)) +
  scale_size_continuous(range = c(1,10)) +
  labs( x= "Site", y = "Genera", size = "Percent of reads assigned each Genera", fill = "X")

bubbleplot
                                        
dat.queryRhizobium <- AbundP[AbundP$X %like% "Frankia", ]        # Extract matching rows with %like%
head(dat.queryRhizobium )
meltedRhizo <- melt(dat.queryRhizobium, id.vars = c("X"))
meltedRhizo[meltedRhizo == 0]<-NA
bubbleplot <- ggplot(meltedRhizo, aes(x = variable, y = X)) +
  geom_point(aes(size = value)) +
  scale_size_continuous(range = c(1,10)) +
  labs( x= "Site", y = "Genera", size = "Percent of reads assigned each Genera", fill = "X")

bubbleplot                               
                               
dat.queryP <- filter(AbundP,X=='Pseudomonas' |X=='Micromonospora'|X =='Bacillus'|X=='Aspergillus'|X=='Streptomyces'|X=='Bacillus'|X=='Burkholderia'|X=='Enterobacter'|X=='Erwinia'|X=='Kushneria'|X=='Paenibacillus'|X=='Ralstonia'|X=='Rhizobium'|X=='Rhodococcus'|X=='Serratia'|X=='Bradyrhizobium'|X=='Salmonella'|X=='Sinomona'|X=='Thiobacillus')
meltedP <- melt(dat.queryP, id.vars = c("X"))
meltedP[meltedP == 0]<-NA
bubbleplot <- ggplot(meltedP, aes(x = variable, y = X)) +
  geom_point(aes(size = value)) +
  scale_size_continuous(range = c(1,10)) +
  labs( x= "Site", y = "Genera", size = "Percent of reads assigned each Genera Phosphorous", fill = "X")
bubbleplot


dat.queryN <- filter(AbundP,X=='Frankia' |X=='Paramesorhizobium'|X =='Mesorhizobium'|X=='Bradyrhizobium'|X=='Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium '|X=='Rhanella1'|X=='Hyphomicrobiales'|X=='Azospirillum'|X=='Ensifer'|X=='Azotobacter'|X=='Paenibacillus'|X=='Derxia'|X=='Herbapirillum')
meltedN <- melt(dat.queryN, id.vars = c("X"))
meltedN[meltedN == 0]<-NA
bubbleplot <- ggplot(meltedN, aes(x = variable, y = X)) +
  geom_point(aes(size = value)) +
  scale_size_continuous(range = c(1,10)) +
  labs( x= "Site", y = "Genera", size = "Percent of reads assigned each Genera Nitrogen", fill = "X")
bubbleplot
