library(randomcoloR)
library(hashmap)
library(Matrix)
library(tictoc)
library(gplots)
library(cluster)
library(clValid) 
library(flexclust)
library(randomcoloR)
library(hashmap)
library(dendextend)
library(randomForest)
library(uwot)
library(e1071)



setwd("/home/andrea/Documents/bioinformatica_condivisa/plasConfig/")
#setwd("/Users/andreabianchetti/Desktop/bioinformatica_condivisa/plasConfig")
evalue_tau <- 1E-10
id_tau <- 90
ali_length_tau <- 500

#command below would load the blast into variable G
#run if the database changes (and comment out the next load command)
#G <-read.table("Klebsiella_plasmids_ncbi_10012019_AAA.txt", header = FALSE,sep = "\t", dec = ".", quote = "",stringsAsFactors = FALSE)
#G was previously loaded and saved in the below session data file
load("original_plasmid_AAA_blast.RData")

#evalue == 0 are set to a minimum of 1E-190
G$V11[which(G$V11==0)] <- 1E-190

#store original variable for record
G_orig <- G

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
#these will be the rows of the preliminary data matrix
query <- sort(unique(G$V1))

#and these will be the columns'
subjects <- sort(unique(G$V2))

#we associate query labels to progressive numbers in a hashmap
#in this way, every time we read a certain query identifier
#we know the number of the row with this query identifier as label
tmp<-c(1:length(query))
#hashmap query name --> row number in the matrix
Q <-hashmap(query,tmp)
tfQ <- "hashmapQ.RData"


save_hashmap(Q, tfQ)

tmp<-c(1:length(subjects))
#hashmap subject name --> column number in matrix
S <-hashmap(subjects,tmp)

tfS <- "hashmapS.RData"


save_hashmap(S, tfS)
#build the matrix of plasmid profiles
#storing evalues
plasprofile_evalue <- Matrix(0,nrow = length(query), ncol = length(subjects), byrow = FALSE,sparse=TRUE)
#storing alignment length
plasprofile_alilength <- Matrix(0,nrow = length(query), ncol = length(subjects), byrow = FALSE,sparse=TRUE)


    #process one query at a time
    for(j in 1:length(query)) {
          tic()
          #prepare the row to store the values for this query
          test_evalue_row <- Matrix(0,nrow = 1, ncol = length(subjects), byrow = FALSE,sparse=TRUE)
          test_alilength_row <- Matrix(0,nrow = 1, ncol = length(subjects), byrow = FALSE,sparse=TRUE)
          
          #take the subset of the blast where query[j] was the query
          B<-G[which(G$V1 == query[j]),]
          
          #now extract all the columns corresponding to the HSPs
          #in this portion of blast
          cid <- S$find(B$V2)
          
          #get evalues for the HSP above
          ev <- B$V11
          #and their alignment lengths
          aliL<-B$V4
          
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
          
          #now we ask for the indices in cid of the values from above
          #this is necessary to associate the evalues in sortev$x to the right places
          cidxx <- match(cidx,cid);
          
          #set the values by using the indices
          test_evalue_row[cidx] <- -log10(sortev$x[cidxx])
          test_alilength_row[cidx] <- -log10(sortAliL[cidxx])
          
          #and put the rows in the full matrix
          plasprofile_evalue[j,]<-test_evalue_row
          plasprofile_alilength[j,]<-test_alilength_row
          toc()
    }



#as we are filtering at a minimum HSP length, there may be
#columns that are all zeros if no alignment
#for a certain plasmid survived
#basically this could only happen if some of the plasmids are shorter than
#the alignment length threshold defined at the beginning of the script

#check if columns with no non-zero value exist
#and eventually remove them
rem <- which(colSums(plasprofile_evalue)==0)
    if(length(rem)>0){
      plasprofile_evalue<-plasprofile_evalue[,-rem]
      subjects<-subjects[-rem]
      plasprofile_alilength <- plasprofile_alilength[,-rem]
      
    }


#the following script will cluster the plasprofile_evalue matrix
#to identify plasmids with an almost identical signature
#only one representative per pattern will be maintained in the
#final matrix to avoid redundant patterns
#these do not pose additional difficulty to their own classification
#but they distort the model and its performances
dim(plasprofile_evalue)
tic()
source("select_plasmids_from_plasprofile_evalue.R")
toc()
dim(plasprofile_evalue)

tmp<-c(1:length(subjects))
#hashmap subject name --> column number in matrix
S <-hashmap(subjects,tmp)
tfS <- "hashmapS.RData"


save_hashmap(S,tfS)
#make a preliminary umap transformation
umeval_tmp <- uwot::umap(as.matrix(plasprofile_evalue),n_neighbors=4, verbose = TRUE, n_threads = 8,ret_model = TRUE)

#p_ev_UMAP <- tempfile("plasprofile_evalue_UMAP",tmpdir=getwd())

#save_uwot(umeval_tmp, file = p_ev_UMAP)

#first perform a kmeans with clvalid to detect the best number of clusters
Kfound <- clValid(umeval_tmp$embedding, nClust = (20:80), clMethods = "kmeans", validation = "internal", maxitems = 2000)

#the following command finds the best number of clusters based on the Silhouette score
OS <- optimalScores(Kfound,measures="Silhouette")$Clusters
nOptimClusters<-as.numeric(attributes(OS)$levels)

#now perform a clustering with the discovered number of clusters
#maybe this command can be removed.
#the clusters discovered here will only be used to train a rf model
#on this very same data
#and it will be the rf model that is used to predict new plasmids
#so the cl1 model will not be used in this approach as it was before
cl1 = kcca(umeval_tmp$embedding, k=nOptimClusters, kccaFamily("kmeans"))

#Knew store the cluster number for each of the plasmids in the row of the profile matrix
Knew <- predict(cl1)
Knew <- data.frame(Knew, row.names = query, stringsAsFactors = F)

#now re-do the umap transformation by giving the cluster labels to guide the manifold approximation
#this is a feature of the uwot package
#the crucial parameter here is y=Knew$Knew because it tells the algorithm to
#use it to guide the transformation
umeval <- uwot::umap(as.matrix(plasprofile_evalue),n_neighbors=4, verbose = TRUE, n_threads = 8, y = Knew$Knew, target_weight = 0.5,ret_model = TRUE)

p_ev_UMAP_rid <- tempfile("plasprofile_evalue_UMAP_rid")

save_uwot(umeval, file = p_ev_UMAP_rid )


#do the same for the umap transformation with 5 components
umeval_tmp5 <- uwot::umap(as.matrix(plasprofile_evalue),n_neighbors=4, verbose = TRUE, n_threads = 8,ret_model = TRUE, n_components = 5)

#p_ev_UMAP5 <- tempfile("plasprofile_evalue_UMAP5")

#save_uwot(umeval_tmp5, file = p_ev_UMAP5)

Kfound5 <- clValid(umeval_tmp$embedding, nClust = (20:80), clMethods = "kmeans", validation = "internal", maxitems = 2000)
OS5 <- optimalScores(Kfound5,measures="Silhouette")$Clusters
nOptimClusters5<-as.numeric(attributes(OS5)$levels)
cl5 = kcca(umeval_tmp$embedding, k=nOptimClusters5, kccaFamily("kmeans"))
Knew5 <- predict(cl5)
Knew5 <- data.frame(Knew5, row.names = query, stringsAsFactors = F)
umeval5 <- uwot::umap(as.matrix(plasprofile_evalue),n_neighbors=4, verbose = TRUE, n_threads = 8, y = Knew5$Knew5, target_weight = 0.5,ret_model = TRUE)

#p_ev_UMAP_rid5 <- tempfile("plasprofile_evalue_UMAP_rid5")

#save_uwot(umeval5, file = p_ev_UMAP_rid5 )

#prepare colors for plotting later
C2 = randomColor(max(Knew$Knew))
plot(jitter(umeval$embedding,amount=.25),col=C2[Knew$Knew])

C5 = randomColor(max(Knew5$Knew5))
plot(jitter(umeval5$embedding,amount=.25),col=C5[Knew5$Knew5])


# ################NOT USED FOR THE MOMENT BUT TO TEST##############################################
# #####and also build a clustering based on the original data matrix (no umap)
# #first perform a kmeans with clvalid to detect the best number of clusters
# Kfound_origdata <- clValid(as.matrix(plasprofile_evalue), nClust = (20:80), clMethods = "kmeans", validation = "internal", maxitems = 2000)
# #the following command finds the best number of clusters based on the Silhouette score
# OS_origdata <- optimalScores(Kfound_origdata,measures="Silhouette")$Clusters
# nOptimClusters_origdata<-as.numeric(attributes(OS)$levels)
# 
# #now perform a clustering with the discovered number of clusters
# #maybe this command can be removed.
# #the clusters discovered here will only be used to train a rf model
# #on this very same data
# #and it will be the rf model that is used to predict new plasmids
# #so the cl1 model will not be used in this approach as it was before
# kcca_origdata = kcca(plasprofile_evalue, k=nOptimClusters_origdata, kccaFamily("kmeans"))
# Knew_origdata <- predict(kcca_origdata)
# Knew_origdata <- data.frame(Knew_origdata, row.names = query, stringsAsFactors = F)
# #################################################################################################

tic()

source("strategy_N2_vBeta_parpar_sel_uwot.R")
toc()




#p_ev_UMAP_rid  name of the umap transform file to load
#file containing the hashmap to load
#hashmapQ.RData                                    
#hashmapS.RData
save(list=ls(),file="plasmidClassification_v3_agg7.RData")
