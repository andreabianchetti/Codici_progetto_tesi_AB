#rm(list=ls())
#setwd("/media/X/andrea_simul/")



#########
#The objective of this script is to take the plasprofile_evalue matrix and the clustering results obtained
#when clustering the blast profiles of plasmids against plasmids
#to cluster novel plasmid signatures obtained when a chromosome is added to the plasmid.
#the genome indeed can distort the original plasmid blast signature and this might therefore preclude
#the obtaining of good predictive models
#Therefore, in this script, we build several matrices simulating the profiles of
#plasmid-chromosome pairs for many associations among a plasmid and a chromosome
#and by using the classification made above for pure plasmid signature, we train a model
#able to remove the effect of the chromosomes therefore still able to classify correctly the plasmids
#




#number of times the procedure is repeated
nrep <- 100

#the plasprofile_evalue matrix can be made bigger, by concatenating several plasprofile_evalue
#matrices after each plasmid has been combined with a different chromosome
#nagg=1 corresponds to using exaclty the plasprofile_evalue matrix after adding the chromosome signals
nagg <- 1

#this is the proportion of rows of the matrix that we leave out
#of the training, and it is used at the end to test the performances of the model on data 
#that when nagg=1, it has never used for training

p_train<-0.8

#the following is a function that make the performance matrix squared
complete_matrix<-function(cm){
  
  missing_cols<-as.numeric(setdiff(rownames(cm),colnames(cm)))
  missing_rows<-as.numeric(setdiff(colnames(cm),rownames(cm)))
  
  if(length(missing_cols)>0){
  for(k in 1:length(missing_cols)){
    if(missing_cols[k]==1){
      cnames=colnames(cm)
      cnames=c(1,cnames)
      cm<-cbind(matrix(0,nrow=nrow(cm),ncol=1),cm)
      colnames(cm)<-cnames
    }else if(missing_cols[k]>ncol(cm)){
      cnames=colnames(cm)
      cnames<-c(cnames,missing_cols[k])
      cm<-cbind(cm,matrix(0,nrow=nrow(cm),ncol=1))
      colnames(cm)<-cnames
    }else{
      cnames=colnames(cm)
      cnames=c(cnames[1:missing_cols[k]-1],missing_cols[k],cnames[missing_cols[k]:ncol(cm)])
      cm<-cbind(cm[,1:missing_cols[k]-1],matrix(0,nrow=nrow(cm),ncol=1),cm[,missing_cols[k]:ncol(cm)])
      colnames(cm)<-cnames
      
    }
    
  }
  }
  
  
  
  
  
  if(length(missing_rows)>0){
    for(k in 1:length(missing_rows)){
      if(missing_rows[k]==1){
        rnames=rownames(cm)
        rnames=c(1,rnames)
        cm<-rbind(matrix(0,nrow=1,ncol=ncol(cm)),cm)
        rownames(cm)<-rnames
      }else if(missing_rows[k]>nrow(cm)){
        rnames=rownames(cm)
        rnames<-c(rnames,missing_rows[k])
        cm<-rbind(cm,matrix(0,nrow=1,ncol=ncol(cm)))
        rownames(cm)<-rnames
      }else{
        rnames=rownames(cm)
        rnames=c(rnames[1:missing_rows[k]-1],missing_rows[k],rnames[missing_rows[k]:nrow(cm)])
        cm<-rbind(cm[1:missing_rows[k]-1,],matrix(0,nrow=1,ncol=ncol(cm)),cm[missing_rows[k]:nrow(cm),])
        rownames(cm)<-rnames
        
      }
      
    }
  }
  if(dim(as.matrix(cm))[1]!=dim(as.matrix(cm))[2]){
    print("Problem with the size of the confusion matrix.")
    
  }
  return(cm)
}


setwd("/home/andrea/Documents/bioinformatica_condivisa/plasConfig/")
#load the blast of the chromosomes against the plasmids
load("Blast_genomes.RData")



num_rows= length(query)#550#  default value, reduce for testing purposes#
#old plasmid vs plasmids blast was called G, rename it

cat("Number of replications:",nrep,"\nSize of the final training matrix:",num_rows*nagg,"\nNumber of rows sampled:",num_rows,"\n")

blastP<-G
tmp<-c(1:length(query))
#carica hashmap usata per associare identificativo query a riga corrispondente nella matrice originale

load_hashmap(tfS)


#filter the blasts as we did before
rem <- which(blastG$V11>evalue_tau)
if(length(rem)>0){
  blastG <- blastG[-rem,]
}

rem <- which(blastP$V11>evalue_tau)
if(length(rem)>0){
  blastP <- blastP[-rem,]
}
rem <- which(blastG$V4<ali_length_tau)
if(length(rem)>0){
  blastG <- blastG[-rem,]
}
rem <- which(blastP$V4<ali_length_tau)
if(length(rem)>0){
  blastP <- blastP[-rem,]
}


#build matrix of genome profiles
#this will be used to merge the plasmid and genome profiles
genprofile_evalue <- Matrix(0,nrow = length(genomes), ncol = length(subjects), byrow = FALSE,sparse=TRUE)

#similarly to the building of the plasprofile_evalue matrix,
#we build the matreix with chromosome profiles
for(j in 1:length(genomes)) {
  tic()
  test_evalue_row <- Matrix(0,nrow = 1, ncol = length(subjects), byrow = FALSE,sparse=TRUE)
  B<-blastG[which(blastG$V1 == genomes[j]),]
  cid <- S$find(B$V2)
  ev <- B$V11
  nanne<-which(is.na(cid))
  cid<-cid[-nanne]
  ev<-ev[-nanne]
  sortev<-sort(ev,index.return=TRUE)
  cid<-cid[sortev$x]
  #unique columns to fill in row rid
  cidx <- unique(cid);
  #now we ask for the indices in cid of the values from above
  #this is necessary to associate the evalues in sortev$x to the right places
  cidxx <- match(cidx,cid);
  test_evalue_row[cidx] <- -log10(sortev$x[cidxx])
  genprofile_evalue[j,]<-test_evalue_row
  toc()
}


library(parallel)
library(iterators)
library(foreach)
library(doParallel)
registerDoParallel(cores=2)

#declare matrices to store the performances
perfSVM<-matrix(NaN,nrow=nrep,ncol=8)
perfRF<-matrix(NaN,nrow=nrep,ncol=8)
perf_testSVM<-matrix(NaN,nrow=nrep,ncol=8)
perf_testRF<-matrix(NaN,nrow=nrep,ncol=8)

#now we want to mix plasmid and genome blasts artificially
#to find models able to clasify the plasmids irrespectively of the
#noise in the signatures introduced by the homologies of the chromosomes with 
#plasmids in the database
test_evalue<-NULL

for(h in 1:nrep){
  tic()
      for(j in 1:nagg){
          
          #rearrange randomly the rows of the genprofile_evalue matrix
          #and resize it into a matrix with nrow=nrow(plasprofile_evalue)
          gp<-genprofile_evalue[sample(1:nrow(genprofile_evalue),nrow(plasprofile_evalue),replace=TRUE),]
          #now take the difference among the two matrices
          Zgp<-which(gp==0)
          Zpl<-which(plasprofile_evalue==0)
          ZERI<-intersect(Zgp,Zpl)
          #set to zero cells that are zero in BOTH matrices
          #in this way cells that are zero in D are necessarily 
          
          #different from zero in the original matrices
          gp[ZERI]<-NaN
          plasprofile_evalue[ZERI]<-NaN
          D<-plasprofile_evalue-gp
          #values that are larger than zero here
          #correspond to evalues that were larger for the plasmid than for the genome
          fromPL<-which(D>0,arr.ind= FALSE)
          valFromPL<-plasprofile_evalue[fromPL]
          fromPL<-which(D>0,arr.ind = TRUE)
          
          #viceversa, the evalue was larger in the genome profiles
          fromGP<-which(D<0,arr.ind = FALSE)
          valFromGP<-genprofile_evalue[fromGP]
          
          fromGP<-which(D<0,arr.ind = TRUE)
          
          #those that are different from zero and take the same value, 
          #can be taken from either matrix  
          fromBoths<-which(D==0,arr.ind = FALSE)
          valFromBoth<-plasprofile_evalue[fromBoths]
          fromBoth<-as.data.frame(which(D==0,arr.ind = TRUE))
          tmp<-sparseMatrix(i=c(fromPL[,1],fromGP[,1],fromBoth[,1]),j=c(fromPL[,2],fromGP[,2],fromBoth[,2]),x=c(valFromPL,valFromGP,valFromBoth),dims=c(nrow(plasprofile_evalue),ncol(plasprofile_evalue)))
              if(j==1){
                test_evalue<-as.matrix(tmp)
              }else{
                test_evalue<-rbind(as.matrix(test_evalue),as.matrix(tmp))
              }
          
        
        
      }
print("Done growing the train matrix")


#also concatenate the cluster id vectors
KNN<-NULL
    for(k in 1:nagg){
        KNN<-c(KNN,Knew$Knew)
    }



#remove the lines for testing
trainidx<-sample(1:nrow(test_evalue),p_train*nrow(test_evalue))
test_evalue_test<-test_evalue[-trainidx,]
test_evalue<-test_evalue[trainidx,]
KNN_test<-KNN[-trainidx]
KNN<-KNN[trainidx]

  #do the umap transform with the uwot package
#and by using the previous transformation
#in parallel with both 2 and 5 components
  umap_train_tmp<-foreach(K=1:2) %dopar%{
    
    if(K==1){  
      #do it with 2 components
      um_tmp<-umap_transform(as.matrix(test_evalue),umeval, verbose = TRUE, n_threads = 6)
      
    }else{
      um_tmp<-umap_transform(as.matrix(test_evalue),umeval5, verbose = TRUE, n_threads = 6)
    }
  }  
     umeval_new<-umap_train_tmp[[1]]
     
      
     # p_ev_UMAP_new <- tempfile("plasprofile_evalue_UMAP_new")
      
      #save_uwot(umeval_new, file = p_ev_UMAP_new)
     
      
    umeval_new5<-umap_train_tmp[[2]]
    
    #p_ev_UMAP_new5 <- tempfile("plasprofile_evalue_UMAP_new5")
    
    #save_uwot(umeval_new5, file = p_ev_UMAP_new_5)
    
      
  print("End umap transformation for this sample")
  
  
  #build the svm and the rf model with the umap transformation
  #with 2 components
  svmrf_2comp<-foreach(K=1:2) %dopar%{

    if(K==1){    
       res<- svm(umeval_new,as.factor(KNN),scale=F,cross=10)
    }else{
      res<- randomForest(umeval_new,as.factor(KNN))
      
      
    }
    
  } 
  
  svm_plasmids<-svmrf_2comp[[1]]
  rf_plasmids<-svmrf_2comp[[2]]
  

  file_nameRFmodel <- "model_RF2C_plasmids.rds"
  saveRDS(rf_plasmids, file = "file_nameRFmodel")
  
  file_nameSVMmodel <- "model_SVM2C_plasmids.rds"
  saveRDS(svm_plasmids, file_nameSVMmodel)
  
  
  
  
  
  #build the svm and the rf model with the uman transformation
  #with 5 components
  svmrf_5comp<-foreach(K=1:2) %dopar%{
    
    if(K==1){    
      res<- svm(umeval_new5,as.factor(KNN),scale=F,cross=10)
    }else{
      res<-randomForest(umeval_new5,as.factor(KNN))
      
      
    }
    
  } 
  svm_plasmids5<-svmrf_5comp[[1]]
  rf_plasmids5<-svmrf_5comp[[2]]
  
  
  
  
  
  
  
  
  
  print("End model training")
  
  #confusion matrices for all models built
  cmSVM <- as.matrix(table(Actual = as.factor(KNN), Predicted = as.factor(svm_plasmids$fitted)))
  cmSVM <- complete_matrix(cmSVM)
  cmRF <- as.matrix(table(Actual = as.factor(KNN), Predicted = as.factor(rf_plasmids$predicted)))
  cmRF <- complete_matrix(cmRF)
  cmSVM5 <- as.matrix(table(Actual = as.factor(KNN), Predicted = as.factor(svm_plasmids5$fitted)))
  cmSVM5 <- complete_matrix(cmSVM5)
  cmRF5 <- as.matrix(table(Actual = as.factor(KNN), Predicted = as.factor(rf_plasmids5$predicted)))
  cmRF5 <- complete_matrix(cmRF5)
  
  
  ######DO IT FOR THE SVM MODEL and umap 2 components
  ns = sum(cmSVM) # number of instances
  nc = nrow(cmSVM) # number of classes
  diags = diag(cmSVM) # number of correctly classified instances per class 
  rowsums = apply(cmSVM, 1, sum) # number of instances per class
  colsums = apply(cmSVM, 2, sum) # number of predictions per class
  p = rowsums / ns # distribution of instances over the actual classes
  q = colsums / ns # distribution of instances over the predicted classes
  accuracySVM = sum(diags) / ns
  precisionSVM = diags / colsums 
  precisionSVM[is.nan(precisionSVM)]<-0
  recallSVM = diags / rowsums 
  
  macroPrecisionSVM = mean(precisionSVM,na.rm=T)
  macroRecallSVM = mean(recallSVM,na.rm=T)
  
  ######DO IT FOR THE SVM MODEL and umap 5 components
  ns = sum(cmSVM5) # number of instances
  nc = nrow(cmSVM5) # number of classes
  diags = diag(cmSVM5) # number of correctly classified instances per class 
  rowsums = apply(cmSVM5, 1, sum) # number of instances per class
  colsums = apply(cmSVM5, 2, sum) # number of predictions per class
  p = rowsums / ns # distribution of instances over the actual classes
  q = colsums / ns # distribution of instances over the predicted classes
  accuracySVM5 = sum(diags) / ns
  precisionSVM5 = diags / colsums 
  precisionSVM5[is.nan(precisionSVM5)]<-0
  recallSVM5 = diags / rowsums 
  
  macroPrecisionSVM5 = mean(precisionSVM5,na.rm=T)
  macroRecallSVM5 = mean(recallSVM5,na.rm=T)
  
  perfSVM[h,1:8]<-c(accuracySVM,macroPrecisionSVM,macroRecallSVM,1/((1/macroPrecisionSVM+1/macroRecallSVM)/2),accuracySVM5,macroPrecisionSVM5,macroRecallSVM5,1/((1/macroPrecisionSVM5+1/macroRecallSVM5)/2))
  
  ######DO IT FOR THE RF MODEL and umap 2 components
  ns = sum(cmRF) # number of instances
  nc = nrow(cmRF) # number of classes
  diags = diag(cmRF) # number of correctly classified instances per class 
  rowsums = apply(cmRF, 1, sum) # number of instances per class
  colsums = apply(cmRF, 2, sum) # number of predictions per class
  p = rowsums / ns # distribution of instances over the actual classes
  q = colsums / ns # distribution of instances over the predicted classes
  accuracyRF = sum(diags) / ns
  precisionRF = diags / colsums 
  precisionRF[is.nan(precisionRF)]<-0
  recallRF = diags / rowsums 
  
  macroPrecisionRF = mean(precisionRF,na.rm=T)
  macroRecallRF = mean(recallRF,na.rm=T)
  
  
  ######DO IT FOR THE RF MODEL and umap 5 components
  
  ns = sum(cmRF5) # number of instances
  nc = nrow(cmRF5) # number of classes
  diags = diag(cmRF5) # number of correctly classified instances per class 
  rowsums = apply(cmRF5, 1, sum) # number of instances per class
  colsums = apply(cmRF5, 2, sum) # number of predictions per class
  p = rowsums / ns # distribution of instances over the actual classes
  q = colsums / ns # distribution of instances over the predicted classes
  accuracyRF5 = sum(diags) / ns
  precisionRF5 = diags / colsums 
  precisionRF5[is.nan(precisionRF5)]<-0
  recallRF5 = diags / rowsums

  macroPrecisionRF5 = mean(precisionRF5,na.rm=T)
  macroRecallRF5 = mean(recallRF5,na.rm=T)
  
  perfRF[h,1:8]<-c(accuracyRF,macroPrecisionRF,macroRecallRF,1/((1/macroPrecisionRF+1/macroRecallRF)/2),accuracyRF5,macroPrecisionRF5,macroRecallRF5,1/((1/macroPrecisionRF5+1/macroRecallRF5)/2))
  
  
  print("End performance calculation for the training set")
  
  
  
  #now focus on the test set and do what we did above
  umap_test_tmp<-foreach(K=1:2) %dopar%{
    
    if(K==1){  
      umtmp<-umap_transform(as.matrix(test_evalue_test),umeval, verbose = TRUE, n_threads = 6)
      
    }else{
      umtmp<-umap_transform(as.matrix(test_evalue_test),umeval5, verbose = TRUE, n_threads = 6)
    }
  }  
      um_test<-umap_test_tmp[[1]]
      um_test5<-umap_test_tmp[[2]]
      
  
  
      print("End umap transformation test sample")
      
  
  svm_pred_test_tmp<-foreach(K=1:2) %dopar%{
    
    if(K==1){    
      pred_classes<-predict(svm_plasmids,um_test)
    }else{
      pred_classes<-predict(svm_plasmids5,um_test5)
      
      
    }
    
  } 
  
  svm_pred_test<-svm_pred_test_tmp[[1]]
  svm_pred_test5<-svm_pred_test_tmp[[2]]
  
  
  
  rf_pred_test_tmp<-foreach(K=1:2) %dopar%{
    
    if(K==1){    
      pred_classes<-predict(rf_plasmids,um_test)
    }else{
      pred_classes<-predict(rf_plasmids5,um_test5)
      
      
    }
    
  } 
  rf_pred_test<-rf_pred_test_tmp[[1]]
  rf_pred_test5<-rf_pred_test_tmp[[2]]
  
  
  print("End Classification test sample")
  cmSVM <- as.matrix(table(Actual = as.factor(KNN_test), Predicted = as.factor(svm_pred_test)))
  cmSVM <- complete_matrix(cmSVM)
  cmRF <- as.matrix(table(Actual = as.factor(KNN_test), Predicted = as.factor(rf_pred_test)))
  cmRF <- complete_matrix(cmRF)
  cmSVM5 <- as.matrix(table(Actual = as.factor(KNN_test), Predicted = as.factor(svm_pred_test5)))
  cmSVM5 <- complete_matrix(cmSVM5)
  cmRF5 <- as.matrix(table(Actual = as.factor(KNN_test), Predicted = as.factor(rf_pred_test5)))
  cmRF5 <- complete_matrix(cmRF5)
  
  
  ######DO IT FOR THE SVM MODEL 2C
  ns = sum(cmSVM) # number of instances
  nc = nrow(cmSVM) # number of classes
  diags = diag(cmSVM) # number of correctly classified instances per class 
  rowsums = apply(cmSVM, 1, sum) # number of instances per class
  colsums = apply(cmSVM, 2, sum) # number of predictions per class
  p = rowsums / ns # distribution of instances over the actual classes
  q = colsums / ns # distribution of instances over the predicted classes
  accuracySVM = sum(diags) / ns
  precisionSVM = diags / colsums 
  precisionSVM[is.nan(precisionSVM)]<-0
  recallSVM = diags / rowsums 
  
  macroPrecisionSVM = mean(precisionSVM,na.rm=T)
  macroRecallSVM = mean(recallSVM,na.rm=T)
 
  
   ######DO IT FOR THE SVM MODEL 5C
  ns = sum(cmSVM5) # number of instances
  nc = nrow(cmSVM5) # number of classes
  diags = diag(cmSVM5) # number of correctly classified instances per class 
  rowsums = apply(cmSVM5, 1, sum) # number of instances per class
  colsums = apply(cmSVM5, 2, sum) # number of predictions per class
  p = rowsums / ns # distribution of instances over the actual classes
  q = colsums / ns # distribution of instances over the predicted classes
  accuracySVM5 = sum(diags) / ns
  precisionSVM5 = diags / colsums 
  precisionSVM5[is.nan(precisionSVM5)]<-0
  recallSVM5 = diags / rowsums 
  
  macroPrecisionSVM5 = mean(precisionSVM5,na.rm=T)
  macroRecallSVM5 = mean(recallSVM5,na.rm=T)
  
  #performances of the test set, SVM first 4 columns for umap transformation 2 comonents
  #last 4 columns for the umap with 5 components
  perf_testSVM[h,1:8]<-c(accuracySVM,macroPrecisionSVM,macroRecallSVM,1/((1/macroPrecisionSVM+1/macroRecallSVM)/2),accuracySVM5,macroPrecisionSVM5,macroRecallSVM5,1/((1/macroPrecisionSVM5+1/macroRecallSVM5)/2))
  
  ######DO IT FOR THE RF MODEL 2C
  ns = sum(cmRF) # number of instances
  nc = nrow(cmRF) # number of classes
  diags = diag(cmRF) # number of correctly classified instances per class 
  rowsums = apply(cmRF, 1, sum) # number of instances per class
  colsums = apply(cmRF, 2, sum) # number of predictions per class
  p = rowsums / ns # distribution of instances over the actual classes
  q = colsums / ns # distribution of instances over the predicted classes
  accuracyRF = sum(diags) / ns
  precisionRF = diags / colsums 
  precisionRF[is.nan(precisionRF)]<-0
  recallRF = diags / rowsums 
  
  macroPrecisionRF = mean(precisionRF,na.rm=T)
  macroRecallRF = mean(recallRF,na.rm=T)
  
  
  ######DO IT FOR THE RF MODEL 5C
  ns = sum(cmRF5) # number of instances
  nc = nrow(cmRF5) # number of classes
  diags = diag(cmRF5) # number of correctly classified instances per class 
  rowsums = apply(cmRF5, 1, sum) # number of instances per class
  colsums = apply(cmRF5, 2, sum) # number of predictions per class
  p = rowsums / ns # distribution of instances over the actual classes
  q = colsums / ns # distribution of instances over the predicted classes
  accuracyRF5 = sum(diags) / ns
  precisionRF5 = diags / colsums 
  precisionRF5[is.nan(precisionRF5)]<-0
  recallRF5 = diags / rowsums
  
  macroPrecisionRF5 = mean(precisionRF5,na.rm=T)
  macroRecallRF5 = mean(recallRF5,na.rm=T)
  print("End performance calculation for the test set")
  
  perf_testRF[h,1:8]<-c(accuracyRF,macroPrecisionRF,macroRecallRF,1/((1/macroPrecisionRF+1/macroRecallRF)/2),accuracyRF5,macroPrecisionRF5,macroRecallRF5,1/((1/macroPrecisionRF5+1/macroRecallRF5)/2))
  
  print(paste("end rep N",h))
  toc()
  if(perf_testRF[h,8]<.4){
    print(paste("end rep N",h))
    
    
  }
  
  
  
}

  colnames(perfSVM) <- c("Accuracy-2c","Precision-2c","Recall-2c","F1-2c","Accuracy-5c","Precision-5c","Recall-5c","F1-5c")
  colnames(perfRF) <- c("Accuracy-2c","Precision-2c","Recall-2c","F1-2c","Accuracy-5c","Precision-5c","Recall-5c","F1-5c")
  
colnames(perf_testSVM) <- c("Accuracy-2c","Precision-2c","Recall-2c","F1-2c","Accuracy-5c","Precision-5c","Recall-5c","F1-5c")
colnames(perf_testRF) <- c("Accuracy-2c","Precision-2c","Recall-2c","F1-2c","Accuracy-5c","Precision-5c","Recall-5c","F1-5c")

par(mfrow=c(2,2))
boxplot(perfSVM,main="Performance score SVM internal",las=2, ylim=c(.6,1))
boxplot(perfRF,main="Performance score RF internal",las=2, ylim=c(.6,1))
boxplot(perf_testSVM,main="Performance score SVM test",las=2, ylim=c(.6,1))
boxplot(perf_testRF,main="Performance score RF test",las=2, ylim=c(.6,1))







lista<-ls()

  save(list=lista[grep("perf",lista)],file="performances_strategy_N2_vBeta_agg3.RData")

  #save(list=lista[grep("rf_plasmid",lista)],file="FINAL_classificationModels_RF.RData")
  
  
  
  
  
  
