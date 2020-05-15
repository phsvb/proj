#make sure ms and ys are good..
##change these as see fit

Xtrains=rbind(M1,M2,M3,M4,M5,M6,M7,M8)
Xtests=M9
Ytrains=rbind(T1Ys,T2Ys,T3Ys,T4Ys,T5Ys,T6Ys,T7Ys,T8Ys)
Ytests=T9Ys
for (i in 1:dim(M1)[2]){ #i in 1:22. Center and scale all features.
  Xtrains[,i]=(Xtrains[,i]-mean(Xtrains[,i]))/sd(Xtrains[,i])  
  Xtests[,i]=(Xtests[,i]-mean(Xtests[,i]))/sd(Xtests[,i])  
}
colnames(Xtrains)=1:22
colnames(Xtests)=1:22
#maybe just do it once
##

resprinter=function(n=1){
  start_time <- Sys.time()
  res=matrix(NA,ncol=3,nrow=n)
  for (resi in 1:n){
    start_time2=Sys.time()
    MP=preProcess(Xtrains,method='pca',thresh=0.95)
    MPa=predict(MP,Xtrains) #PCAcomp for training set
    Ndims=dim(MPa)[2]
    MP2=preProcess(Xtests,method='pca',pcaComp=Ndims)
    MP2a=predict(MP2,Xtests) #PCAcomp for test set
    Str=fastICA(Xtrains,n.comp=Ndims) #train
    StrS=Str$S #S_train signals
    SteS=MP2a%*%Str$W #S_test from W_train
#    SteS2=t(Str$W%*%t(MP2a))
    #alternatively try #SteS=t(Str$W%*%t(Mpa2))
    Mtrains=matrix(NA,nrow=dim(Ytrains)[1],ncol=Ndims)
    Mtests=matrix(NA,nrow=dim(Ytests)[1],ncol=Ndims)
 #   Mtests2=matrix(NA,nrow=dim(Ytests)[1],ncol=Ndims)
    for(i in 1:dim(Ytrains)[1]){
      Mtrains[i,]=log(diag(var(StrS[1:750+(i-1)*750,]))) #band power
    }
    for(i in 1:dim(Ytests)[1]){
      Mtests[i,]=log(diag(var(SteS[1:750+(i-1)*750,]))) #band power
    }
  #  for(i in 1:dim(Ytests)[1]){
   #   Mtests2[i,]=log(diag(var(SteS2[1:750+(i-1)*750,]))) #band power
    #}
    Trainset=cbind(Ytrains[,2],Mtrains)
    Testset=cbind(Ytests[,2],Mtests)
    rf1=randomForest(y=as.factor(Trainset[,1]),x=Trainset[,-1])
  ## fastICA results
    ficate1=mean(predict(rf1,newdata=Testset[,-1])==Testset[,1]) 
  ##COROICA
    Ctr=coroICA(MPa,partitionsize=48*750,groupsize=750*48*6,n_components=Ndims,max_iter = 10000,pairing = 'neighbouring')
    CtrS=Ctr$Shat
    MStrains=matrix(NA,nrow=dim(Ytrains)[1],ncol=Ndims)
    MStests=matrix(NA,nrow=dim(Ytests)[1],ncol=Ndims)
    MStests2=matrix(NA,nrow=dim(Ytests)[1],ncol=Ndims)
    for(i in 1:dim(Ytrains)[1]){
      MStrains[i,]=log(diag(var(CtrS[1:750+(i-1)*750,]))) #band power
    }
    test2=t(Ctr$V%*%t(MP2a))  #looks silly! but works. should be doing this i imagine!!
    for(i in 1:dim(Ytests)[1]){
      MStests[i,]=log(diag(var(test2[1:750+(i-1)*750,]))) #band power
    }
    temp1=cbind(Ytrains[,2],MStrains)
    temp2=cbind(Ytests[,2],MStests)
    rfc1=randomForest(y=as.factor(temp1[,1]),x=temp1[,-1]) #OOB error rate 62%
    cicate=mean(predict(rfc1,newdata=temp2[,-1])==temp2[,1]) #test accuracy 32. despite non-convergence. issa okay?
    cat('i=',resi,'conver=',Ctr$converged,'\n')
    cat('runtime=',Sys.time()-start_time2,'\n')
    res[resi,]=c(ficate1,cicate,Ctr$converged)
  }
  end_time <- Sys.time()
  cat('runningtime=',end_time - start_time,'\n')
  colnames(res)=c('fastICA','coroICA','coroica conver')
  return(res) 
}
#0 if false. 1 if true
n20=resprinter(n=20) #testsubject is T9
#residentsleeper timezzz

n10=resprinter(n=10)
cat(n,Ctr$converged)

cat(1-2,'\n',2)

n10
#didnt converge once..
n10
mean(n10[,1])
mean(n10[,2])
mean(n10[,3])

var(n10[,3])
var(n10[,1])


c(1,2,FALSE,TRUE,FALSE)


n10
#each try takes about a minute?
mean(n10[,2])

colMeans(n20)
#31
#n100=resprinter(100)