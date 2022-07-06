##### Full Script for exploring the flexibility in social structure of East African chimpanzees. 
#The main part of the dataset is repeated including only adult individuals and using subset randomisation. The step required for subset randomisation can be seen in the ESM script at the bottom

#####Packages required#####
library(ape)
library(asnipe)
library(igraph)
library(vegan)
library(effsize)
library(ggplot2)
library(data.table)
library(adehabitatHR)

####Data upluad and creation of social network graphs####
#every step is repeated for Waibira dataset
#upload gbi data for Sonso
sgbi1<-as.matrix(read.csv(file.choose(),header=FALSE))
sgbi2<-as.matrix(read.csv(file.choose(),header=FALSE))
sgbi3<-as.matrix(read.csv(file.choose(),header=FALSE))
sgbi<-cbind(sgbi1,sgbi2,sgbi3)
sgbi<-as.matrix(sgbi)
Sonso<-t(sgbi)

#upload attributes for each individual
Sattributes<-read.csv(file.choose(),header = TRUE)
#define columns of gbi data by indiviual ID
colnames(Sonso)<-Sattributes$ID

#create network object
networkS<-get_network(association_data=Sonso, data_format="GBI")

#build adjacency graph for simple ratio network
netS<-graph.adjacency(networkS, mode="undirected", weighted=TRUE, diag=FALSE)

#create random networks through data stream permutation
permutations<-10000 #use 1000 permutations
randomS <- network_permutation(Sonso, association_matrix=networkS, permutations=permutations)


#calculate mean weight in random networks
randomES<-rep(NA, 10000)
for (i in 1:10000){
  net.r<-graph.adjacency(randomS[i,,], mode="undirected", weighted=TRUE, diag=FALSE)
  edges<-get.data.frame(net.r)
  randomES[i] <- mean(edges$weight)} 

mean(randomES)

#create network with only above random edges
Scut.off<-mean(randomES)
dil_netS <- delete_edges(netS, E(netS)[weight<Scut.off])
plot(dil_netS,edge.width= E(netS)$weight*5,
     vertex.size=10,
     vertex.label.cex=0.5) 
     
#calculate network measures Sonso, check if they are different from random networks
#transitivity 
obsSt<-transitivity(netS)

	#create object that holds results for all the random networks
random.St<-rep(NA, 10000)
for (i in 1:10000){
  net.r<-graph.adjacency(randomS[i,,], mode="undirected", weighted=TRUE, diag=FALSE)
  random.St[i] <- transitivity(net.r)} 
  
  #plot histogram of results from random networks with observed result in as red line
hist(random.St, breaks=1,col="grey", main = "Sonso", xlab = "Transitivity")
abline(v=obsSt, col="red")#add line for observed measure

	#compare results from random and observed networks
sum(random.St >= obsSt)/10000


#mean strength 
obsSstr<-mean(strength(netS))
range(strength(netS))
mean(strength(netS))
sd(strength(netS))

	#create object that holds results for all the random networks
random.Sstr<-rep(NA, 10000)
for (i in 1:10000){
  net.r<-graph.adjacency(randomS[i,,], mode="undirected", weighted=TRUE, diag=FALSE)
  random.Sstr[i] <- mean(strength(net.r))} 
  
	 #plot histogram of results from random networks with observed result in as red line
hist(random.Sstr, breaks=50,col="grey", main = "Sonso", xlab = "Mean Strength")
abline(v=obsSstr, col="red")

	#compare results from random and observed networks
sum(random.Sstr >= obsSstr)/10000


#modularity: optimal clustering
clusS_O<- cluster_optimal(netS) #calculate clusterng algorithm (Louavain's)
obsmodS_O<-modularity(clusS_O) #calculate modularity

	#create object that holds results for all the random networks
random.modS_O<-rep(NA, 10000)
for (i in 1:10000){
  net.r<-graph.adjacency(randomS_O[i,,], mode="undirected", weighted=TRUE, diag=FALSE)
  clus<- cluster_optimal(net.r)
  random.modS_O[i] <- modularity(clus)} 
  
	#plot histogram of results from random networks with observed result in as red line
hist(random.modS_O, breaks=50,col="grey", main = "Sonso", xlab = "Modularity")
abline(v=obsmodS_O, col="red")

	#compare results from random and observed networks
sum(random.modS_O >= obsmodS_O)/10000
  
#modularity: Louvain's clustering
clusS_L<- cluster_louvain(netS) #calculate clusterng algorithm (Louavain's)
obsmodS_L<-modularity(clusS_L) #calculate modularity

#create object that holds results for all the random networks
random.modS_L<-rep(NA, 10000)
for (i in 1:10000){
  net.r<-graph.adjacency(randomS_L[i,,], mode="undirected", weighted=TRUE, diag=FALSE)
  clus<- cluster_louvain(net.r)
  random.modS_L[i] <- modularity(clus)} 

#plot histogram of results from random networks with observed result in as red line
hist(random.modS_L, breaks=50,col="grey", main = "Sonso", xlab = "Modularity")
abline(v=obsmodS_L, col="red")

#compare results from random and observed networks
sum(random.modS_L >= obsmodS_L)/10000

####Create networks for Waibira Core (Waibira_A) and Waibira periphery (Waibira_B) clusters####
#create Waibira_A network
Waibira_A<-subset(Waibira, select=c("ALF", "BEN","DAU", "DOU", "FID", "GER","MAC", "MAS","SAM", "TAL", "TRS"))

#create network object
networkWA<- get_network(association_data=Waibira_A, data_format="GBI")

#build adjacency graph for simple ratio network
netWA<-graph.adjacency(networkWA, mode="undirected", weighted=TRUE, diag=FALSE)


#create random networks through data stream permutation
randomWA <- network_permutation(Waibira_A, association_matrix=networkWA, permutations=permutations)

#calculate mean random edge weight 
randomEWA<-rep(NA, 10000)
for (i in 1:10000){
  net.r<-graph.adjacency(randomWA[i,,], mode="undirected", weighted=TRUE, diag=FALSE)
  edges<-get.data.frame(net.r)
  randomEWA[i] <- mean(edges$weight)} 

mean(randomEWA)

#create above chance edges network
WAcut.off<-mean(randomEWA)
dil_netWA <- delete_edges(netWA, E(netWA)[weight<WAcut.off])
plot(dil_netWA,edge.width= E(netWA)$weight*5,
     vertex.size=10,
     vertex.label.cex=0.5) 
     

#create Waibira_B network
Waibira_B<-subset(Waibira, select=c("ABO", "ARD", "CHN", "ILA", "KAS", "KEV","LAF", "LAN", "MOR", "MUG", "URS"))

#create network object
networkWB<- get_network(association_data=Waibira_B, data_format="GBI")

#build adjacency graph for simple ratio network
netWB<-graph.adjacency(networkWB, mode="undirected", weighted=TRUE, diag=FALSE)

#create random networks through data stream permutation
randomWB <- network_permutation(Waibira_B, association_matrix=networkWB, permutations=permutations)

#calculate mean random edge weight 
randomEWB<-rep(NA, 10000)
for (i in 1:10000){
  net.r<-graph.adjacency(randomWB[i,,], mode="undirected", weighted=TRUE, diag=FALSE)
  edges<-get.data.frame(net.r)
  randomEWB[i] <- mean(edges$weight)} 

mean(randomEWB)

#create above chance edges network
WBcut.off<-mean(randomEWB)
#choose Waibira full mean randomE so that we keep unconnected nodes unconnected (only one dyad within Waibira B are connected in the full Waibira network)
dil_netWB <- delete_edges(netWB, E(netWB)[weight<Wcut.off])
plot(dil_netWB,edge.width= E(netWB)$weight*5,
     vertex.size=10,
     vertex.label.cex=0.5) 

####Yearly sociograms####
#Create yearly sociograms - same for all years and both communities, here shown example for year 2015-2016 in Waibira
wgbi15<-as.matrix(read.csv(file.choose(),header=FALSE))
Waibira15<-t(wgbi15)

#define columns of gbi data by indiviual ID
colnames(Waibira15)<-Wattributes$ID

#create network object
networkW15<- get_network(association_data=Waibira15, data_format="GBI")

#build adjacency graph for simple ratio network
netW15<-graph.adjacency(networkW15, mode="undirected", weighted=TRUE, diag=FALSE)

#create random networks through data stream permutation
randomW15 <- network_permutation(Waibira15, association_matrix=networkW15, permutations=permutations)

#calculate mean random edge weight 
randomEW15<-rep(NA, 10000)
for (i in 1:10000){
  net.r<-graph.adjacency(randomW15[i,,], mode="undirected", weighted=TRUE, diag=FALSE)
  edges<-get.data.frame(net.r)
  randomEW15[i] <- mean(edges$weight)} 

mean(randomEW15)

#create above chance edges network
Wcut.off15<-mean(randomEW15)
dil_netW15 <- delete_edges(netW15, E(netW15)[weight<Wcut.off15])
plot(dil_netW15,edge.width= E(netW15)$weight*5,
     vertex.size=10,
     vertex.label.cex=0.5) 
#modularity across years - repeated for all years and both communities
#Waibira
clusW15<- cluster_optimal(netW_15)
mod_W15<-modularity(clusW15)

	#create object that holds results for all the random networks
random.modW_15<-rep(NA, 10000)
for (i in 1:10000){
  net.r<-graph.adjacency(randomW_15[i,,], mode="undirected", weighted=TRUE, diag=FALSE)
  clus<- cluster_optimal(net.r)
  random.modW_15[i] <- modularity(clus)} 

	#compare results from random and observed networks
sum(random.modW_15 >= mod_W15)/10000

####Core-Periphery Analysis####
#plot distribution of degree strength for Sonso and Waibira - this is the step Alex suggested but it doens't work very well... 
hist(strength(netS))
hist(strength(dnetW))


#Here I want to calculate the probability of edges being present within core-core (Pcc), core-periphery (Pcp), and periphery-periphery (Ppp) dyads. core-pheryphery networks should show that Pcc > Pcp > Ppp

#create subsets of individuals from the Waibira graph that includes only edges with above chance SRI
W_da<-V(dil_netW)[c("ALF", "BEN", "DOU", "FID", "GER","MAC", "MAS", "TAL", "TRS")]
W_db<-V(dil_netW)[c("ABO", "ARD", "CHN", "DAU","ILA", "KAS", "KEV","LAF", "LAN", "MOR", "MUG", "SAM", "URS")]

#create subsets of individuals from the Waibira graph that includes oall possible edges (saturated graph)
W_a<-V(netW)[c("ALF", "BEN", "DOU", "FID", "GER","MAC", "MAS", "TAL", "TRS")]
W_b<-V(netW)[c("ABO", "ARD", "CHN", "DAU","ILA", "KAS", "KEV","LAF", "LAN", "MOR", "MUG", "SAM", "URS")]

#induce a subgraph of Waibira_A (core) from unsaturated graph that includes only edges with above mean SRI
net_da<-induced.subgraph(graph=dil_netW,vids = W_da)

#induce a subgraph of Waibira_B (periphery) from unsaturated graph that includes only edges with above mean SRI
net_db<-induced.subgraph(graph=dil_netW,vids = W_db)

#induce a subgraph of Waibira_A (core) from saturated graph
net_a<-induced.subgraph(graph=netW,vids = W_a)

#induce a subgraph of Waibira_B (periphery) from saturated graph
net_b<-induced.subgraph(graph=netW,vids = W_b)

#count number of edges in unsaturated network between individuals in:
nA<-gsize(net_da)#core-core
nB<-gsize(net_db)#periphery-periphery
nAB<-(gsize(dil_netW))-(nA+nB)#core-periphery

#count number of possible edges between individuals in:
tA<-gsize(net_a)#core-core
tB<-gsize(net_b)#periphery-periphery
tAB<-(gsize(netW))-(tA+tB)#core-periphery


#calculate edge probability between:
Pcc<-nA/tA#core-core
Ppp<-nB/tB#periphery-periphery
Pcp<-nAB/tAB#core-periphery

####distance between graphs using D dissimilarity####
setwd("/Users/gb64/Desktop/")
write_graph(dil_netS, "/Users/gb64/Desktop/elS.txt")
write_graph(dil_netW, "/Users/gb64/Desktop/elW.txt")
write_graph(dil_netWA, "/Users/gb64/Desktop/elWA.txt")
write_graph(dil_netWB, "/Users/gb64/Desktop/elWB.txt")


D_full<-D("elS.txt","elW.txt",0.45,0.45,0.1)
D_S_WA<-D("elS.txt","elWA.txt",0.45,0.45,0.1)
D_S_WB<-D("elS.txt","elWB.txt",0.45,0.45,0.1)

#comparing degree between subgroup
Wattributes$Group<-c("B","A","B","A","B","A","A","A","A","B","B","B","B","B","A", "A", "B", "B", "A", "A", "A", "B" ) #assign group membership to each individual based on Louvain's clustering algorithm (above) and add to attribute table
Wattributes$Colour<-c("blue","red","blue","red","blue","red","red","red","red","blue","blue","blue","blue","blue","red","red","blue","blue","red","red","red","blue") #assign colour to each individual based on their subgroup membership and add to attribute table

#### Linear Model using Double Permutation ####
# Specify metric
metric <- "DEGREE"
# calculate observed network
network <- get_network(Waibira, data_format="GBI",
                       association_index="SRI")
# Calculate observed metric (degree)
degrees <- rowSums(network)
# Do randomisation (as above, permutations should be >=1000)
networks.perm <- network_permutation(Waibira, data_format="GBI",
                                     association_matrix=network,  permutations=10000)
# Now calculate the same metric on all the random networks
degrees.rand <- apply(networks.perm,1,function(x) { rowSums(x)})
# Now substract each individual's median from the observed
degree.controlled <- degrees - apply(degrees.rand,1,median)

#### Now use degree.controlled for any later test. For example, to related against a trait:
# Make a trait
trait <- Wattributes$Group
# get the coefficient of this:
coef <- summary(lm(degree.controlled~trait))$coefficients[2,3]
# Compare this to a node permutation
# (here just randomising the trait values)
# note this should be done >= 1000 times
n.node.perm <- 100000
coefs.random <- rep(NA, n.node.perm)
for (i in 1:n.node.perm) {
  trait.random <- sample(trait)
  coefs.random[i] <- summary(lm(degree.controlled~trait.random))$coefficients[2,3]
}
# calculate P value (note this is only one sided)
P <- sum(coef <= coefs.random)/n.node.perm

#####Waibira Range calculations of range overlap between Waibira core and Waibira periphery clusters#####
#use either All_ranging_CBlock.csv or FirstDaily_Ranging_CBlock.csv depending on whether calculating only with first daily observation or all observations
wbr<-read.csv(file.choose(),header = TRUE) 
head(wbr)

#subset full dataset for the parties in which each individual is present
abo<-subset(wbr,wbr$ABO=="1")
alf<-subset(wbr,wbr$ALF=="1")
ben<-subset(wbr,wbr$BEN=="1")
chn<-subset(wbr,wbr$CHN=="1")
dou<-subset(wbr,wbr$DOU=="1")
ger<-subset(wbr,wbr$GER=="1")
ila<-subset(wbr,wbr$ILA=="1")
kas<-subset(wbr,wbr$KAS=="1")
kev<-subset(wbr,wbr$KEV=="1")
mac<-subset(wbr,wbr$MAC=="1")
mor<-subset(wbr,wbr$MOR=="1")
sam<-subset(wbr,wbr$SAM=="1")
tal<-subset(wbr,wbr$TAL=="1")
trs<-subset(wbr,wbr$TRS=="1")
urs<-subset(wbr,wbr$URS=="1")

##add a column includig the ID of the individual 
abo$ID<-ifelse(abo$ABO=="1","ABO")
alf$ID<-ifelse(alf$ALF=="1","ALF")
ben$ID<-ifelse(ben$BEN=="1","BEN")
chn$ID<-ifelse(chn$CHN=="1","CHN")
dou$ID<-ifelse(dou$DOU=="1","DOU")
ger$ID<-ifelse(ger$GER=="1","GER")
ila$ID<-ifelse(ila$ILA=="1","ILA")
kas$ID<-ifelse(kas$KAS=="1","KAS")
kev$ID<-ifelse(kev$KEV=="1","KEV")
mac$ID<-ifelse(mac$MAC=="1","MAC")
mor$ID<-ifelse(mor$MOR=="1","MOR")
sam$ID<-ifelse(sam$SAM=="1","SAM")
tal$ID<-ifelse(tal$TAL=="1","TAL")
trs$ID<-ifelse(trs$TRS=="1","TRS")
urs$ID<-ifelse(urs$URS=="1","URS")

#create dataframe with everyone's ranging patterns by individual instead of by party
all<-rbind(abo,alf,ben,chn,dou,ger,ila,kas,kev,mac,mor,sam,trs,tal,urs)
head(all)

#create smaller data frame from the one above in which only the lat, lon and ID columns are present
comp<-cbind.data.frame(all$lat,all$lon,all$ID)
colnames(comp) <- c("lat", "lon","ID")
head(comp)

#create file of coordinates only
xy<-cbind(comp$lat,comp$lon)
xy<-SpatialPoints(xy)
#create spatial point data frame
spdf <- SpatialPointsDataFrame(coords = xy, data = comp,bbox=NULL,match.ID=TRUE,
                               proj4string = CRS(as.character(NA)))
head(spdf)
spdf@data <- spdf@data %>% 
  dplyr::select(3) #keeps column 2 in the spdf object.


#kernel overlap function: creates matrix reporting the percentage of home range over lap between each dyad. This function was ran on each of the kernals reported (here 95% of individual home range)
k<-kerneloverlap(spdf[,1],method = "HR",percent = 95,conditional = FALSE)
k#this is the matrix, we export it from R and rearrange this matrix in excel to include 5 columns. The first two columns include the ID of each individual in the dyad, the 3rd includeds the percent overlap between their home ranges, the 4th and 5th column indicate the subgroup each individual belongs to and the final column indicates if the two individuals are in the same or different subgroups 

write.csv(k, "/Users/gb64/Desktop/1block_a_5k.csv")#to create dataframe for t.tests in next step

#t.test comparing overlap between dyads in the same or different subgroups: this t-test was repeated for each kernal measures
#use any/all of overlap csvs separately 
#(e.g. All_overlap5csv = overlap of 5% core area between males with all observations included/ FirstDaily_overlap5.csv = same but with just first daily observaton)
perc<-read.csv(file.choose(), header=TRUE)
t.test(perc$K~perc$match)

####subgroup stability across the years Waibira using Levenshtein distance ####
clusW15<- cluster_optimal(netW_15) 
clusW16<- cluster_optimal(netW_16) 
clusW17<- cluster_optimal(netW_17) 
clusW18<- cluster_optimal(netW_18)

#create list of individuals in Waibira core for each year
Wattributes$ID 
WA_15<-c("BDFGHIOPSU")
WA_16<-c("BDGHIOPSTU")
WA_17a<-c("ACDFHJMNRST")
WA_17b<-c("BGIKOPU")
WA_18<-c("BCFGHIJMNOPRSU")

#compare membership in the core between years
library(stringdist)
stringdist(WA_15, WA_16)
stringdist(WA_15, WA_17a)
stringdist(WA_15, WA_17b)
stringdist(WA_15, WA_18)
stringdist(WA_16, WA_17a)
stringdist(WA_16, WA_17b)
stringdist(WA_17a, WA_18)
stringdist(WA_17b, WA_18)


#####Figures#####

#Figure 1: Histogram of Waibira and Sonso network measures from null models and observed values (also created for ESM figures S2 and S5 to show results from subset randomisation and adult date respectively)
par(mfrow=c(2,4))
hist(random.St, breaks=1,col="grey", main = "A)", xlab = "Transitivity")
par(xpd = TRUE)
lines(x = c(obsSt,obsSt), y = c(0, par('usr')[4]), col="red")

hist(random.Sstr, breaks=50,col="grey", main = "B)", xlab = "Mean Strength")
par(xpd = TRUE)
lines(x = c(obsSstr,obsSstr), y = c(0, par('usr')[4]), col="red")

hist(random.modS_O, breaks=50,col="grey", main = "C)", xlab = "Modularity (Optimal)")
par(xpd = TRUE)
lines(x = c(obsmodS_O,obsmodS_O), y = c(0, par('usr')[4]), col="red")

hist(random.modS_L, breaks=50,col="grey", main = "C)", xlab = "Modularity (Louvain's)")
par(xpd = TRUE)
lines(x = c(obsmodS_O,obsmodS_O), y = c(0, par('usr')[4]), col="red")


hist(random.Wt, breaks=1,col="grey", main = "D)", xlab = "Transitivity")
par(xpd = TRUE)
lines(x = c(obsWt,obsWt), y = c(0, par('usr')[4]), col="red")

hist(random.Wstr, breaks=50,col="grey", main = "E)", xlab = "Mean Strength")
par(xpd = TRUE)
lines(x = c(obsWstr,obsWstr), y = c(0, par('usr')[4]), col="red")

hist(random.modW_O, breaks=50,col="grey", main = "F)", xlab = "Modularity (Optimal)")
par(xpd = TRUE)
lines(x = c(obsmodW_O,obsmodW_O), y = c(0, par('usr')[4]), col="red")

hist(random.modW_L, breaks=50,col="grey", main = "F)", xlab = "Modularity (Louvain's)")
par(xpd = TRUE)
lines(x = c(obsmodW_L,obsmodW_L), y = c(0, par('usr')[4]), col="red")

dev.copy2pdf(file= "/Users/gb64/Desktop/Flexibility in Chimpanzee Social Network/BEH_revisions/Figures/Figure1ms_20220517.pdf")

#Figure 2: Sociograms for Waibira and Sonso
#Waibira

V(dil_netW)$community <- clusW$membership#set group membership
V(dil_netW)$shape <- c("square","circle")[V(dil_netW)$community]#set group membership by shap

#set plot space
par(mfrow=c(1,2),mar=c(1,1,1,1))
#plot Waibira sociogram
plot(dil_netW, vertex.color="dimgrey",#plot sociograph
     edge.width=E(netW)$weight*7,edge.color="grey"(alpha=0.30, level=0.25),vertex.size=(strength(netW)*3),
     vertex.label=NA, vertex.frame.color="white")
title("Waibira", cex.main=1)

     
#Sonso

#plot Sonso sociogram
plot(dil_netS, vertex.color="dimgrey",#plot sociograph
     edge.width=E(netS)$weight*7,edge.color="grey"(alpha=0.30, level=0.25),vertex.size=(strength(netS)*3),
     vertex.label=NA, vertex.frame.color="white")
title("Sonso", cex.main=1)

dev.copy2pdf(file= "C:/Users/gb64/Desktop/SNA/Figure2ms_20220517.pdf")

#Figure 3 (and ESM fig S3 and S4): box plot illustrating range of individual strenght in Waibira core and periphery clusters (repeated with subset randomisation)
par(mfrow=c(1,2),cex.lab=1.5,mar=c(4,4,1,1), tck=0.01, mgp=c(2.5,0.5,0),las=1)

# Plot original boxplot
plot(Wattributes$DEGREE.r~factor(Wattributes$Group),xlab="Cluster",ylab="Individual Strength",cex.axis=1, tck=0.01)
text(par('usr')[1] - (par('usr')[2]-par('usr')[1])/5,par('usr')[4] - (par('usr')[4]-par('usr')[3])/15, "A", cex=1.5) 
dev.copy2pdf(file= "/Users/gb64/Desktop/R_temp/Figure3ms_20220520.pdf")

# Plot resulting distribution histogram for double permutation model (ESM S3 and S4b)
hist <- hist(coefs.random,xlim=c(min(coefs.random),max(coefs.random)),col="black",xlab="Coefficient value",ylab="Frequency",breaks=100,cex.axis=1,main="", tck=0.01)
segments(coef,0,coef,max(hist$counts),col="red",lwd=3)
box()
text(par('usr')[1] - (par('usr')[2]-par('usr')[1])/5,par('usr')[4] - (par('usr')[4]-par('usr')[3])/15, "B", cex=1.5) 

dev.copy2pdf(file= "/Users/gb64/Desktop/R_temp/hist_lm1_20220519.pdf")

#Figure 4: ranging visualisation for Waibira Core and Periphery (also repeated to include all parties observed and not just 1st party of each day for ESM Fig S7)
library(dplyr)
library(adehabitatHR)
library(gplots)

gapa<-read.csv(file.choose(),check.names = FALSE,sep = ",",comment.char = "#")
rnames<-gapa[,1]
mat_gapa<-data.matrix(gapa[,2:ncol(gapa)])
rownames(mat_gapa)<-rnames

my_pallete1<-colorRampPalette(c("skyblue","yellow","tomato"))(n=299)


heatmap.2(mat_gapa,
          main = "Waibira_B",
          notecol="black",
          density.info = "none",
          trace = "none",
          margins = c(12,9),
          col=my_pallete1,
          dendrogram="row",
          srtCol=360,
          Colv = "NA",
          Rowv=FALSE)

dev.copy2pdf(file= "/Users/gb64/Desktop/Flexibility in Chimpanzee Social Network/BEH_revisions/Figures/Figure3bms_1block.pdf")

#Figure 5
par(mfrow=c(2,2))#create plot space to include 4 plots

V(netW_15)$community <- clusW15$membership#create object including individual by group
V(netW_16)$community <- clusW16$membership
V(netW_17)$community <- clusW17$membership
V(netW_18)$community <- clusW18$membership

#create a triangle vertex shape option
mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }

  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}

add_shape("triangle",
                 plot=mytriangle)


#membership each year by shape
V(netw15.aa)$shape <- c("square", "circle", "triangle")[V(netW_15)$community ]
V(netw16.aa)$shape <- c("square", "circle", "triangle")[V(netW_16)$community]
V(netw17.aa)$shape <- c("circle", "triangle","square")[V(netW_17)$community]
V(netw18.aa)$shape <- c("square","circle", "triangle")[V(netW_18)$community]

#plot all sociograms
par(mfrow=c(2,2), mar=c(3, 3, 2, 0.5), mgp=c(2, 0.8, 0))

plot(dil_netW_15, vertex.color="dimgrey",#plot sociogram and repeat for each year
     edge.width=E(netW_15)$weight*7,edge.color="grey"(alpha=0.3, level=0.3),vertex.size=(strength(netW_15)*3),vertex.shape=V(netw15.aa)$shape ,
     vertex.label=NA, vertex.frame.color="white")
title("A", cex.main=1.5)


plot(dil_netW_16, vertex.color="dimgrey",#plot sociogram and repeat for each year
     edge.width=E(netW_16)$weight*7,edge.color="grey"(alpha=0.3, level=0.3),vertex.size=(strength(netW_16)*3),vertex.shape=V(netw16.aa)$shape ,
     vertex.label=NA, vertex.frame.color="white")
title("B", cex.main=1.5)


plot(dil_netW_17, vertex.color="dimgrey",#plot sociogram and repeat for each year
     edge.width=E(netW_17)$weight*7,edge.color="grey"(alpha=0.3, level=0.3),vertex.size=(strength(netW_17)*3),vertex.shape=V(netw17.aa)$shape ,
     vertex.label=NA, vertex.frame.color="white")
title("C", cex.main=1.5)


plot(dil_netW_18, vertex.color="dimgrey",#plot sociogram and repeat for each year
     edge.width=E(netW_18)$weight*7,edge.color="grey"(alpha=0.3, level=0.3),vertex.size=(strength(netW_18)*3),vertex.shape=V(netw18.aa)$shape ,
     vertex.label=NA, vertex.frame.color="white")
title("D", cex.main=1.5)

dev.copy2pdf(file= "C:/Users/msstt/Desktop/SNA/FigureS4_20220517.pdf")

#####ESM Analyses and Figures#####

#####Subset randomisation##### - all analysis other than this step is the same as individual randomisation, this step is repeated for each community, but this example uses Sonso

#add observation number added to dataframe once raw data are added to workspace and individual attributes are added
Sonso_work <- cbind(Sonso, "observation"=1:nrow(Sonso)) 
Sonso_work<-as.data.frame(Sonso_work)

#Gives every 6 consecutive parties a unique identifier/number
Sonso_work<-Sonso_work %>%
  mutate(N_split = floor(observation/6))

#splits data into groups based on their unique number
Split_S<-split(Sonso_work, f = Sonso_work$N_split)

#creates new data-frame including one randomly sampled party from each group defined above
Sub_S<- data.frame(matrix(ncol = 11, nrow = 0))
colnames(Sub_S)<-Sattributes$ID
head(Sub_S)
for(i in 1:length(Split_S)){
  row<-sample_n(Split_S[[i]], 1)
  Sub_S[i,]<-row[1,1:11]
}

####mantel test Figure S1####
par(mfrow=c(4,3))
mants1<-mantel.test(networkS_15,networkS_16, nperm=1000, graph = T, alternative = "greater",main="a) 2015/2016")
mants2<-mantel.test(networkS_15,networkS_17, nperm=1000, graph = T, alternative = "greater",main="b) 2015/2017")
mants3<-mantel.test(networkS_15,networkS_18, nperm=1000, graph = T, alternative = "greater",main="c) 2015/2018")
mants4<-mantel.test(networkS_16, networkS_17, nperm=1000, graph = T, alternative = "greater",main="d) 2016/2017")
mants5<-mantel.test(networkS_16, networkS_18, nperm=1000, graph = T, alternative = "greater",main="e) 2016/2018")
mants6<-mantel.test(networkS_17, networkS_18, nperm=1000, graph = T, alternative = "greater",main="f) 2017/2018")

mantw1<-mantel.test(networkW_15,networkW_16, nperm=1000, graph = T, alternative = "greater",main="g) 2015/2016")
mantw2<-mantel.test(networkW_15,networkW_17, nperm=1000, graph = T, alternative = "greater",main="h) 2015/2017")
mantw3<-mantel.test(networkW_15,networkW_18, nperm=1000, graph = T, alternative = "greater",main="i) 2015/2018")
mantw4<-mantel.test(networkW_16, networkW_17, nperm=1000, graph = T, alternative = "greater",main="j) 2016/2017")
mantw5<-mantel.test(networkW_16, networkW_18, nperm=1000, graph = T, alternative = "greater",main="k) 2016/2018")
mantw6<-mantel.test(networkW_17, networkW_18, nperm=1000, graph = T, alternative = "greater",main="l) 2017/2018")

dev.copy2pdf(file= "/Users/gb64/Desktop/Flexibility in Chimpanzee Social Network/BEH_revisions/Figures/SuppFig1_20210520.pdf")

####Hierarchy calculations Figure S6####
library(EloRating)
pgs<-read.csv(file.choose(),header = TRUE)#add Waibira_mm_pg.csv (dataset with male pant grunts)

#create elo-rating seqence
elo<-elo.seq(winner = pgs$Winner, loser = pgs$Loser, Date = pgs$Date,draw=NULL, presence = NULL,
             startvalue = 1000, k=200, normprob = TRUE, init = "average",
             intensity = NULL, iterate = 0, runcheck = FALSE, progressbar = FALSE)
#create object with individuals who are included in the study            
m<-c("ABO", "ALF", "ARD", "BEN","CHN", "DAU", "DOU", "FID", "GER", "ILA", "KAS", "KEV", "LAF", "LAN", "MAC", "MAS", "MOR", "MUG", "SAM", "TAL", "TRS", "URS" )
#extract the elo-rating of the individuals included in the study
extract_elo(elo, "2019-09-30",IDs = m)

eloplot(elo, ids=m,from="start", to = "2019-09-30", color = TRUE )#get plot of rank progression 
dev.copy2pdf(file= "/Users/gb64/Desktop/SuppFig6_20220517.pdf")