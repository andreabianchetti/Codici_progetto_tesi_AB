library(randomcoloR)
library(hashmap)
library(Matrix)
library(tictoc)
library(gplots)
library(cluster)
library(clValid) 
library(colormap)
library(flexclust)
library(randomcoloR)
library(hashmap)
library(dendextend)
library(randomForest)
library(uwot)
library(e1071)
library(colorspace)
library(doParallel)
library(foreach)
library(iterators)
library(modeltools)
library(parallel)


setwd("/home/andrea/Documents/andrea/klebsiella_patric/blast")
tfS <- "hashmapS.RData"
modelRF <- "/home/andrea/Documents/bioinformatica_condivisa/plasConfig/plasprofile_evalue_UMAP7bbb328730de"
#load("/home/andrea/Documents/bioinformatica_condivisa/plasConfig/Tot_environment.RData")
#pth<-file.path(cat("/home/andrea/Documents/bioinformatica_condivisa/plasConfig/",tfS,sep=""))
S<-load_hashmap("/home/andrea/Documents/bioinformatica_condivisa/plasConfig/hashmapS.RData")


#setwd("/Users/andreabianchetti/Desktop/bioinformatica_condivisa/klebsiella_patric/blast")

evalue_tau <- 1E-10
id_tau <- 90
ali_length_tau <- 500


#evalue == 0 are set to a minimum of 1E-190
minEval <- 1E-190


listaBlast <- read.table("lista_all_blasts.txt",stringsAsFactors = F, header=F)


pp <- Matrix(0,nrow = length(listaBlast$V1), ncol = S$size(), byrow = FALSE,sparse=TRUE)


for( i in 1:length(listaBlast$V1)){
  tic()
  evalue_row <- Matrix(0,nrow = 1, ncol = S$size(), byrow = FALSE,sparse=TRUE)
  
  #genomesID[i,1] <- listaBlast$V1[i]
  G <- read.table(listaBlast$V1[i],stringsAsFactors = F, header=F)
  
  
  #remove lines of the blast with an evalue>evalue_tau
  rem <- which(G$V11>evalue_tau)
  if(length(rem)>0){
    G <- G[-rem,]
  }
  
  #similarly, remove lines for alignments that are too short
  rem <- which(G$V4<ali_length_tau)
  if(length(rem)>0){
    G <- G[-rem,]
  }
  
  w <- which(G$V11<minEval)
  if(length(rem)>0){
    G$V11[w] <- minEval
  }
  
  #these will be the rows of the preliminary data matrix
  query <- sort(unique(G$V1))
  
  #and these will be the columns'
  subjects <- sort(unique(G$V2))
  
  
  
  cid <- S$find(G$V2)
  rem<-which(is.na(cid))
  cid <- cid[-rem]
  #get evalues for the HSP above
  ev <- G$V11
  ev<-ev[-rem]
  #and their alignment lengths
  aliL<-G$V4
  aliL<-aliL[-rem]
  #now sort the evalues
  sortev<-sort(ev,index.return=TRUE)
  #and sort accordingly the columns
  cid<-cid[sortev$ix]
  #and the alignment lengths
  sortAliL<-aliL[sortev$ix]
  
  #for each column we only store the first HSP
  #this is done by taking the unique columns
  #unique returns the first index, therefore, as we sorte the alignments by evalue,
  #we are sure that for a certain query and a certain subject, this is the most significant HSP
  cidx <- unique(cid);
  
  
  cidxx <- match(cidx,cid);
  
  #set the values by using the indices
  evalue_row[cidx] <- -log10(sortev$x[cidxx])
  pp[i,] <- evalue_row
  toc()
  
    
    }

###########################################################################################################################################################################################

    # umeval_model<-load_uwot(rf_plasmids)   #
  umeval_model<-load_uwot(file = p_ev_UMAP_rid)

  umeval_new_tmp <- uwot::umap_transform(as.matrix(pp),umeval_model)
 
  
  RF_genvspl <- readRDS("model_RF2C_plasmids.rds")
  
 # pred_new <- predict(rf_plasmids, as.matrix(pp)) #
  pred_new <- as.numeric(predict(RF_genvspl, umeval_new_tmp))
  #prepare colors for plotting later
  C2 = randomColor(max(pred_new))
  plot(jitter(umeval_new_tmp,amount=.25),col=C2[pred_new])
  
  plot((umeval_new_tmp),col=C2[pred_new])
  
  #################
  
  nrep<-100
  allcl<-matrix(NA,nrow=nrow(umeval_new_tmp),ncol=nrep)
  
  
  for(i in 1:nrep){
    tic()
    umeval_noisy <- jitter(umeval_new_tmp,amount=.25)
    pred_noisy <- as.numeric(predict(RF_genvspl, umeval_noisy))
    allcl[,i] <- (pred_noisy)
    toc()
  }
  
  ###############################################

  
  cl_inclusion_p <- matrix(0,nrow=nrow(umeval_new_tmp),ncol=max(pred_new))
  
  for(i in 1:nrow(umeval_new_tmp)){
    tmp<-as.matrix(table(allcl[i,]))
    c<-as.numeric(rownames(tmp))
      for(n in 1:length(c)){
        cl_inclusion_p[i,c[n]] <- tmp[n]
        
      }
    
        
    
  }
  
  
  finalCL <- matrix(NA,nrow=nrow(umeval_new_tmp),ncol=2)
    for(i in 1:nrow(umeval_new_tmp)){
      mx <- max(cl_inclusion_p[i,]);
      mx_id <- which(cl_inclusion_p[i,]==mx)
      finalCL[i,]<- c(mx_id[1],mx[1]/nrep)
  }

  
  
  
  perfect<-which(finalCL[,2]>=.95)
  library(rgl)
  plot3d(umeval_new_tmp[,1],umeval_new_tmp[,2],jitter(1,amount=.25),col="black",size=1.5)
  plot3d(umeval_new_tmp[perfect,1],umeval_new_tmp[perfect,2],jitter(2,amount=.25),col=C2[finalCL[perfect,1]],size=5,add=TRUE)
  



  plot(umeval_new_tmp[,1],umeval_new_tmp[,2],col=C2[finalCL[,1]],cex=finalCL[,2],pch=16)
  
######################
  ########
  #####################




####################**********************************************************######################################################***********************************************************
####################**********************************************************######################################################***********************************************************
metadata<-read.csv("/home/andrea/Documents/bioinformatica_condivisa/klebsiella_patric/PATRIC_genome.csv",stringsAsFactors = FALSE,header=TRUE,fill=TRUE)
#metadata<-read.csv("/Users/andreabianchetti/Desktop/bioinformatica_condivisa/klebsiella_patric/PATRIC_genome.csv")


rem<-which(is.na(as.numeric(metadata$Collection.Date)))
metadata<-metadata[-rem,]

genomes<- gsub(".fna_blast_vs_ncbi_plasmids.txt","",listaBlast$V1)

m<-match(genomes,metadata$Genome.ID)

rem <- which(is.na(m))
m <- m[-rem]
genomes<-genomes[-rem]
umeval_new_tmpSEL <- umeval_new_tmp[-rem,]
metadataSel <- metadata[m,]
finalCLSEL <- finalCL[-rem,]
years<-as.numeric(metadataSel$Collection.Date)



uY<-sort(unique(years))
w<-which(metadataSel$Collection.Date == uY[1])

anno<-rep(1,length(w))

rgl.bg(color = "grey")
plot3d(umeval_new_tmpSEL[w,1],anno,umeval_new_tmpSEL[w,2],col = C2[finalCLSEL[w,1]])



for(i in 2:length(uY)-1){
  
  w<-which(metadataSel$Collection.Date == uY[i])
  anno<-rep(i,length(i))
  
  plot3d(umeval_new_tmpSEL[w,1],anno,umeval_new_tmpSEL[w,2],col = C2[finalCLSEL[w,1]],add=TRUE)
  
  
  
  
}
rgl.postscript("plot3d.pdf","pdf") 



############

save(list=ls(),file="model_PlasConfig.RData")

############

save(list=ls(),file="Classificator_BLAST_new_genomes_vs_olddata.RData")
