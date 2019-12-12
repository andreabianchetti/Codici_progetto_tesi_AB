
#calculate a dendrogram of plasmids based on the similarities in their signatures
hc <- hclust(dist(plasprofile_evalue), "ave")
dend <- as.dendrogram(hc)

library(colorspace)
library(colormap)
#cut the dendrogram at this height to put together
#plasmids with extremely similar signatures
#THe objective of this analysis is to reduce the redundancy in the rows
#and then columns
mycut<-cutree(dend, k= NULL,h=3)
dend1 <- color_branches(dend, h = 1)  

plot(dend1, main = "Plasmids clusters before modeling")

#just as an exercise, cut the dendrogram at many increasing thresholds
#and count the number of clusters
ncl<-matrix(0,nrow=50,ncol=2)
for(i in 1:50){
  mycut<-cutree(dend,h=(i/10))
  ncl[i,]<-c(length(unique(mycut)),i/10)
}

plot(ncl[,2],ncl[,1],xlab="Cutree height threshold",ylab="Number of clusters")

#we want to put together only plasmids that are very similar or identical
#after visual inspection
#we choose to cut at 3
#this indicates there are very similar plasmids in the dataset, we should also remove them 
#from the columns (or the database), as they provide redundant information
#now we sample 1 representant per cluster and build the new reduced matrix for
#model training

mycut<-cutree(dend,h=3)
#declare a vector to record the plasmids to be taken for continuing the approach
totake<-matrix(0,nrow=length(query),ncol=1)

#define the unique cluster names
ucl<-unique(mycut)

#for each cluster, identify plasmids it is composed of
#and get one random plasmid
for(i in 1:length(ucl)){
  w<-as.numeric(which(mycut==ucl[i]))
  if(length(w)>1){
    #if the cluster contains more than one plasmid
    #sample only one plasmid, randomly
    s<-sample(w,1)
  }else{
    #if there is only one plasmid in this cluster, then take it
    s=w
  }
  
  #add 1 for plasmid to take
  totake[s]<-1
  
  
}


#simple check. the number of ones in totake must equal the number of clusters
#in ucl
if(colSums(totake)!=length(ucl)){
  
  print("There is a problem")
}

#identify plasmids to be taken
tt<-which(totake==1)

#store original version of the plasmid_profile matrix
plasprofile_evalue_original<-plasprofile_evalue

#get only lines for plasmids with 1 in totake
plasprofile_evalue<-plasprofile_evalue[tt,]
#do the same with the query name vector
query<-query[tt]

######now do the same for the columns of the matrix
##exactly the same thing on columns
hc <- hclust(dist(t(plasprofile_evalue)), "ave")
dend <- as.dendrogram(hc)

ncl<-matrix(0,nrow=50,ncol=2)
for(i in 1:50){
  mycut<-cutree(dend,h=(i/10))
  ncl[i,]<-c(length(unique(mycut)),i/10)
}

plot(ncl[,2],ncl[,1],xlab="Cutree height threshold",ylab="Number of clusters")


mycut<-cutree(dend,h=3)

totake<-matrix(0,nrow=length(subjects),ncol=1)

ucl<-unique(mycut)

for(i in 1:length(ucl)){
  w<-as.numeric(which(mycut==ucl[i]))
  if(length(w)>1){
    s<-sample(w,1)
  }else{
    s=w
  }
  totake[s]<-1
  
}


if(colSums(totake)!=length(ucl)){
  
  print("There is a problem")
}

tt<-which(totake==1)



plasprofile_evalue<-plasprofile_evalue[,tt]
dim(plasprofile_evalue)


subjects<-subjects[tt]

#this script therefore will change the plasprofile_evalue matrix
#by removing rows that have very similar (if not identical) signatures
#this for instance because the db contains the same plasmid more than once
#e.g. should be quite common for clinical strains


#save(list=ls(),file="plasmidClassification_v3.RData")