### R code for Research Methods of Microbial Biogeography by Peng 

### 1. Upload
env <- read.csv("env.csv", row.names = 1)
geo <- read.csv("geo.csv", row.names = 1) 
otu <- read.csv("otu.csv", row.names = 1)
taxon <- read.csv("taxon.csv", row.names = 1) 
clim <- read.csv("clim.csv", row.names = 1) 
bnti <- read.csv("nti.csv", row.names = 1) 

library(pacman)
p_load(ggplot2,reshape2,patchwork,vegan,geosphere,psych,corrplot)


### 2. Species composition
otu_abun <- otu/colSums(otu) 
p<- aggregate(otu_abun, by=list(taxon$Phylum), sum) 
rownames(p) <- p$Group.1
p <- p[,-1]
top <- names(head(sort(rowSums(p), decreasing = T), 57)) 
top_10<- c("Proteobacteria", "Actinobacteria", "Acidobacteria", "Chloroflexi", "Gemmatimonadetes","Bacteroidetes", "Planctomycetes","Firmicutes", "Thermomicrobia", "Nitrospirae")
taxon$Phylum <- as.character(taxon$Phylum)
taxon$Phylum[!(taxon$Phylum)%in%top_10] <- "Others"
p_top <- aggregate(otu_abun, by=list(taxon$Phylum), sum) #获取top10
rownames(p_top) <- p_top[,1] 
p_top <- p_top[,-1]
p_top <- p_top[order(rowSums(p_top)),] 
q_top <- t(p_top)
q_top <- as.data.frame(q_top)
q_top$sample <- rownames(q_top)
q <- melt(q_top,ID="names")
colnames(q)[names(q)=="variable"]<-"Taxa" 
q$factor <- "forest"

colors<-c("grey50","darkolivegreen3","gold","dodgerblue4","darkseagreen",
          "chartreuse4","darkorange","burlywood2","brown3","#984EA3","cyan3")
ggplot(q, aes( x = sample, y = value, fill = Taxa))+
  geom_bar(position = "fill", stat = "identity")+ 
  theme_bw()+
  scale_fill_manual(values=colors)+  
  scale_y_continuous(expand = c(0,0))+ 
  labs(x="",y="Relative Abundance",fill="Phylum")+
  theme(text=element_text(size=12),
        axis.text.y=element_text(size=12,color = "black"), 
        axis.text.x=element_text(size=12,color = "black",angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12))+ 
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent')) + 
  guides(fill=guide_legend(keywidth = 1.3, keyheight = 2))


### 3. Latitude diversity pattern
otu <- as.data.frame(t(otu))
richness <- specnumber(otu) #计算丰富度
shannon <- diversity(otu,index = "shannon") #计算shannon多样性
diversity <- data.frame(richness,shannon) 
aa <- cbind(geo, diversity)
summary(lm(aa$richness ~ aa$lat))
summary(lm(aa$shannon ~ aa$lat))

p1 <- ggplot(aa,aes(x = lat, y = richness))+
          geom_point()+
          geom_smooth(method = "lm",alpha = 0.2)
p2 <- ggplot(aa,aes(x = lat, y = shannon))+
          geom_point()+
          geom_smooth(method = "lm",alpha = 0.2)
plot <- p1 + p2
plot


### 4. Factors affecting bacterial alpha-diversity
diversity <- data.frame(richness,shannon)
cor_soil <- corr.test(diversity, env[,1:14], method = "spearman", adjust = "none")
cor_soil_r <- cor_soil$r
cor_soil_p <- cor_soil$p
cor_soil_p[cor_soil_p >= 0.05] <- -1
cor_soil_p[cor_soil_p < 0.05 & cor_soil_p >= 0] <- 1
cor_soil_p[cor_soil_p == -1] <- 0
corr_soil <- cor_soil_r*cor_soil_p
cor_clim <- corr.test(diversity, clim[,1:19],method = "spearman",adjust = "none")
cor_clim_r <- cor_clim$r
cor_clim_p <- cor_clim$p
cor_clim_p[cor_clim_p >= 0.05] <- -1
cor_clim_p[cor_clim_p < 0.05 & cor_clim_p >= 0] <- 1
cor_clim_p[cor_clim_p == -1] <- 0
corr_clim <- cor_clim_r*cor_clim_p

col<-colorRampPalette(c("#77AADD","#4477AA","#FFFFFF", "#EE9988","#BB4444"))
corrplot(corr_soil, method = 'color', addCoef.col = 'black',  number.cex = 0.8, rect.col = "black",addgrid.col = "black",col = col(200),tl.col = 'black', cl.pos = "n")
corrplot(corr_clim, method = 'color', addCoef.col = 'black', number.cex = 0.8, rect.col = "black",addgrid.col = "black",col = col(200),tl.col = 'black', cl.pos = "n")


### 5. Distance decay relationship
d.geo <- distm(geo, fun = distHaversine) 
dist.geo <- as.dist(d.geo)
dist_geo <- as.data.frame(as.vector(dist.geo))
dist_geo <- dist_geo/1000 
colnames(dist_geo) <- "dist_geo"
dist.otu <- vegdist(otu, method = "bray")
dist_otu <- as.data.frame(as.vector(dist.otu))
colnames(dist_otu) <- "dist_otu"
data <- data.frame(dist_geo, dist_otu)
data$dist_otu <- 1-data$dist_otu
data$dist_otu <- data$dist_otu * 100 
summary(lm(data$dist_otu ~ data$dist_geo))

ggplot(data, aes(x = dist_geo,y = dist_otu)) + 
geom_point() + 
geom_smooth(method = "lm",alpha = 0.2) + 
labs(x = "Geographic Distance (Km)", y = "Community Similarity (%) ")


### 6. Factors affecting bacterial community structure
dist.otu <- vegdist(otu, method = "bray")
soil <- NULL
for (i in 1:ncol(env)) {
aa <- mantel(dist.otu, dist(env[,i], method = "euclidean"), method = "pearson", permutations = 9999, na.rm = TRUE)
soil <- rbind(soil,c(colnames(env)[i],aa$statistic, aa$signif))
}
climate <- NULL
for (i in 1:ncol(clim)) {
 aa <- mantel(dist.otu, dist(clim[,i],method = "euclidean"), method = "pearson", permutations = 9999, na.rm = TRUE)
 climate <- rbind(climate,c(colnames(clim)[i],aa$statistic, aa$signif))
}


### 7. Bacterial community β diversity
pcoa <- cmdscale(dist.otu, k = 3, eig = TRUE)
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
site <- data.frame(pcoa$points)[1:2]
site <- cbind(site,env$CEC)
colnames(site) <- c('PCoA1', 'PCoA2',"CEC")

ggplot(site, aes(PCoA1, PCoA2)) +
  geom_point(aes(color = CEC),size = 3)+  
  scale_color_gradient(low = "blue",high = "red") + 
  labs(x = paste('PCoA1: ', round(100 * pcoa_eig[1], 2), '%'), 
       y = paste('PCoA2: ', round(100 * pcoa_eig[2], 2), '%'))


### 8. Community assembly
## Calculation method of Beta_NTI 
library(picante)
commu <- read.csv("otu.csv",fileEncoding = "UCS-2LE",row.names = 1)
phylo <- read.tree("OTUs.tre") # Phylogenetic tree of each OTU
Beta_NTI<-function(phylo,comun,beta.reps=999){
   comun=t(comun) 
   match.phylo.comun = match.phylo.data(phylo, t(comun)) 
   beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.comun$data),cophenetic(match.phylo.comun$phy),abundance.weighted=T))
   rand.weighted.bMNTD.comp = array(c(-999), dim=c(ncol(match.phylo.comun$data), ncol(match.phylo.comun$data),beta.reps))
   for (rep in 1:beta.reps) {
       rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.comun$data),taxaShuffle(cophenetic(match.phylo.comun$phy)),abundance.weighted=T,exclude.conspecifics = F))
       print(c(date(),rep))
   }
   weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.comun$data),ncol=ncol(match.phylo.comun$data))
   for(columns in 1:(ncol(match.phylo.comun$data)-1)) {
          for(rows in (columns+1):ncol(match.phylo.comun$data)) {
                  rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
                  weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals)
                  rm("rand.vals")
          }
   }
  rownames(weighted.bNTI) = colnames(match.phylo.comun$data);
  colnames(weighted.bNTI) = colnames(match.phylo.comun$data);
  return(as.dist(weighted.bNTI))
}
bnti <- Beta_NTI(phylo,comun,beta.reps=999)

bnti <- bnti[lower.tri(bnti)]
sto <- length(which(abs(bnti)<2))
det <- length(which(abs(bnti)>2))
assem <- c(sto,det)
pie(assem, labels=c("stochastic process","deterministic process"),col = c("#F8766D","#00BFC4"),radius = 1)
