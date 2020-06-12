#load(TMsRAW.RData) #RAW
#load(TYs.RData)
Mlist=list(M1,M2,M3,M4,M5,M6,M7,M8,M9)
Ylist=list(T1Ys,T2Ys,T3Ys,T4Ys,T5Ys,T6Ys,T7Ys,T8Ys,T9Ys)
resmatrix=matrix(NA,nrow=9,ncol=3)

main4=function(N=22){ #N=dims
#  Mlist=list(M1,M2,M3,M4,M5,M6,M7,M8,M9)
#  Ylist=list(T1Ys,T2Ys,T3Ys,T4Ys,T5Ys,T6Ys,T7Ys,T8Ys,T9Ys)
  Mlist=list(rbind(M1,N1),rbind(M2,N2),rbind(M3,N3),rbind(M4,N5),rbind(M5,N5),rbind(M6,N6),rbind(M7,N7),rbind(M8,N8),rbind(M9,N9))
  Ylist=list(rbind(T1Ys,E1Ys),rbind(T2Ys,E2Ys),rbind(T3Ys,E3Ys),rbind(T4Ys,E4Ys),rbind(T5Ys,E5Ys),rbind(T6Ys,E6Ys),rbind(T7Ys,E7Ys),rbind(T8Ys,E8Ys),rbind(T9Ys,E9Ys))
  resmatrix=matrix(NA,nrow=9,ncol=3)
  for(inde in 1:9){
    LM=Mlist
    LY=Ylist
    Xtests=LM[[inde]]
    Ytests=LY[[inde]] 
    LM[[inde]]=NULL #removing the ith, 
    LY[[inde]]=NULL
    Xtrains=do.call(rbind,LM)
    Ytrains=do.call(rbind,LY)
    #
    #removing artifacts for training set.
    #may screw some partition/group-size for coroICA
    #    indy=rep(NA,length(which(Ytrains[,3]==1)))
    #    for(i in 1:length(which(Ytrains[,3]==1))){
    #      indy[(1:750)+(i-1)*750]=(which(Ytrains[,3]==1)[i]*750-749):(which(Ytrains[,3]==1)[i]*750)
    #    }
    #    Xtrains=Xtrains[-indy,]
    #    Ytrains=Ytrains[-which(Ytrains[,3]==1),]
    ##############^######maybe unsafe
    #CAR
    #   for(i in 1:dim(Xtests)[1]){
    #      Xtests[i,]=Xtests[i,]-mean(Xtests[i,])
    #    }##
#############center+scale:::
#    for (i in 1:22){ #i in 1:22. Center and scale all features.
#      Xtrains[,i]=(Xtrains[,i]-mean(Xtrains[,i]))/sd(Xtrains[,i])  
#      Xtests[,i]=(Xtests[,i]-mean(Xtests[,i]))/sd(Xtests[,i])  
#    }
    colnames(Xtrains)=1:22
    colnames(Xtests)=1:22
    Ndims=N
#    Mtrains=matrix(NA,nrow=dim(Ytrains)[1],ncol=Ndims)
#    Mtests=matrix(NA,nrow=dim(Ytests)[1],ncol=Ndims)
#    MP=preProcess(Xtrains,method="pca",pcaComp = Ndims)
#    pcaTrainComp=predict(MP,Xtrains)
    #   MP2=preProcess(Xtests,method="pca",pcaComp = Ndims)
#    pcaTestsComp=predict(MP,Xtests)

    #    Str=fastICA(Xtrains,n.comp=Ndims) #train
    #   StrS=Str$S                        #S_train signals
    #  SteS=(Xtests%*%Str$K)%*%Str$W     #test signals
    ##CAR
    for(i in 1:dim(Xtrains)[1]){
        Xtrains[i,]=Xtrains[i,]-mean(Xtrains[i,])
    }
    for(i in 1:dim(Xtests)[1]){
        Xtests[i,]=Xtests[i,]-mean(Xtests[i,])
    }
    Mtrains=matrix(NA,nrow=dim(Ytrains)[1],ncol=Ndims)
    Mtests=matrix(NA,nrow=dim(Ytests)[1],ncol=Ndims)
    MP=preProcess(Xtrains,method="pca",pcaComp = Ndims)
    pcaTrainComp=predict(MP,Xtrains)
 #   MP2=preProcess(Xtests,method="pca",pcaComp = Ndims)
    pcaTestsComp=predict(MP,Xtests)
#    for(j in 1:Ndims){
#      pcaTrainComp[,j]=signal::filter(coefs$b, coefs$a, pcaTrainComp[,j])
#      pcaTestsComp[,j]=signal::filter(coefs$b, coefs$a, pcaTestsComp[,j])
#    }
    #CAR
#    for(i in 1:dim(Xtrains)[1]){
#        Xtrains[i,]=Xtrains[i,]-mean(Xtrains[i,])
#    }##
    #CAR
#    for(i in 1:dim(Xtests)[1]){
#        Xtests[i,]=Xtests[i,]-mean(Xtests[i,])
#    }##
    #######filter before?
#    
#    
    ##try PCA
    ##
    ##
#    for(j in 1:22){
#      Xtrains[,j]=signal::filter(coefs$b, coefs$a, Xtrains[,j])
#      Xtests[,j]=signal::filter(coefs$b, coefs$a, Xtests[,j])
#    }
 ##   Ctr=coroICA(Xtrains,partitionsize=48*750/10,groupsize=750*48*6,n_components=Ndims,max_iter = 10000,pairing = 'neighbouring',minimize_loss=T)
    #?
    Ctr=coroICA(pcaTrainComp,partitionsize=48*750/10,groupsize=750*48*6,n_components=Ndims,max_iter = 10000,pairing = 'neighbouring')
  #  Ctr$converged
    StrS=Ctr$Shat
    SteS=pcaTestsComp%*%t(Ctr$V)
 #   SteS=Xtests%*%t(Ctr$V)
    ###FILTER
    #  for(i in 1:dim(Xtrains)[1]){
    #    Xtrains[i,]=Xtrains[i,]-mean(Xtrains[i,])
    #  }
    
    for(j in 1:Ndims){
      StrS[,j]=signal::filter(coefs$b, coefs$a, StrS[,j])
      SteS[,j]=signal::filter(coefs$b, coefs$a, SteS[,j])
    }
    ##CAR
#    for(i in 1:dim(StrS)[1]){
#      StrS[i,]=StrS[i,]-mean(StrS[i,])
#    }
#    for(i in 1:dim(SteS)[1]){
#      SteS[i,]=SteS[i,]-mean(SteS[i,])
#    }
    ###band power
    for(i in 1:dim(Ytrains)[1]){
      Mtrains[i,]=log(diag(var(StrS[1:750+(i-1)*750,]))) #band power
    }
    
    for(i in 1:dim(Ytests)[1]){
      Mtests[i,]=log(diag(var(SteS[1:750+(i-1)*750,]))) #band power
    }
    
    Trainset=cbind(Ytrains[,2],Mtrains)
    Testset=cbind(Ytests[,2],Mtests)
    B=200            #bootstrap: is it tho? ye i think
    n=nrow(Trainset)
    L=numeric(B)
    cat('conv=',Ctr$converged)
    
    #    L=rep(0,1) #10 for now.. change to bigger eventually
    for (k in 1:B){
      i= sample(n,n,replace=T)
      S1=sda(Trainset[i,-1],as.factor(Trainset[i,1]),verbose=F)
      L[k]=mean(predict(S1,Testset[,-1])$class==Testset[,1])
      #      S1=sda(Trainset[,-1],L=as.factor(Trainset[,1]),verbose=F)
      #     L[k]=mean(predict(S1,Testset[,-1])$class==Testset[,1]) #decent!
    }
    cat('conv=',Ctr$converged,inde)
    resmatrix[inde,]=c(mean(L),mean(L)+c(-1,1)*1.96*sd(L)/(sqrt(length(L))))
  }
  return(resmatrix)
}
#no pca
#test2=main4()
#colMeans(test2)

#no pca
test1=main4()
colMeans(test1)

#yes pca #no car tho xD
test21=main4(21)
colMeans(test21)

#yes pca #yes car
test21c=main4(21)
colMeans(test21c)

#pca no good! no need to whitten
minl=main4()
minl
colMeans(minl)

ffirst=main4()
ffirst
colMeans(ffirst)#32.7%  

faf=main4()
faf
colMeans(faf) #35% 

faf14=main4(14)
faf14
colMeans(faf14) #31

faf21=main4(21)
faf21
colMeans(faf21) #35.2%


#only car, no scale center
t10306=main4() #bis 25k
t10306
colMeans(t10306) #35.4

#scale center car later
t103061=main4() #bis 100k
t103061
colMeans(t103061) #35.5

#########
t1030621=main4(21) #bis 100k
t1030621
colMeans(t1030621)#36

preProcess()

#/1
coro1CAR=main4()#w. car
coro1CAR
mean(coro1CAR[,1]) #33.6

coro1=main4()#no car
mean(coro1[,1]) #33.6

##/10
coro10car=main4() #no car 
coro10car
mean(coro10car[,1])#35.3..


coro10=main4() #w. car
coro10
mean(coro10[,1])#35.9

##basically the same no matter what?
#/2
coro2=main4()
coro2
mean(coro2[,1])#33.8%

coro2CAR=main4()
coro2CAR
mean(coro2CAR[,1])#33.2%
#/4
coronopca=main4(22) #no pca no car
coronopca
mean(coronopca[,1]) #34.4% USING /4
coronopca2=main4(22) #no pca no car
coronopca2
mean(coronopca2[,1]) #34.4% USING /4


coroCAR=main4() #no pca +car
coroCAR[,1]
mean(coroCAR[,1]) #34.9% USING /4


#/8
coro8=main4()
coro8
mean(coro8[,1])#35%


#/16
coro16=main4()
coro16
mean(coro16[,1])#34.6%



corosda14nocar=main4(14)
corosda14nocar
mean(corosda14nocar[,1])
##
corosda5nocar=main4(5)
corosda5nocar
mean(corosda5nocar[,1])

qplot(1:9,corosda5nocar[,1])

qplot(1:9,coronopca[,1],geom = 'line')

#with V instead of t(V)
corosda14v=main4(14)
mean(corosda14v[,1]) 

corosda14[,1]
corosda14v[,1]

#   resmatrix[inde,]=c(mean(L),mean(L)+c(-1,1)*1.96*sd(L)/(sqrt(length(L))))
#not converging with 22, try 9

nopca1
mean(nopca1[,1])

########seems not tobad^
corosda22=main4()
mean(corosda22[,1]) #28

corosda9=main4(9)
mean(corosda9[,1])

corosda5=main4(5)
mean(corosda5[,1])


corosda14=main4(14)
mean(corosda14[,1])


nopca1=main4(22) #28% bad!
noicaica=main4(22) #28% bad!
mean(noicaica[,1])
corosda5=main4(5) #28% bad!
mean(corosda5[,1])

corosda9=main4(9) #28% bad!
mean(corosda9[,1])


corosda16=main4(16)
mean(corosda16[,1])

corosda14=main4(14)
mean(corosda14[,1])

noicaica=main4(22) #28% bad!
mean(noicaica[,1])

coroSDA=main3(22) #22
coroSDA #26
coroSDA2=main3(20) #22
coroSDA2#26
mean(coroSDA[,1])

Ndims=10

head(Xtrains)
head(Xtests)
prcomp(Xtrains)

t2=prcomp(Xtrains,center=F,scale.=F,)
?prcomp

t1=preProcess(Xtrains,method = 'pca',pcaComp = Ndims)
t1x=predict(t1,Xtrains)
head(t1x)
dim(t2$x)

t2x=preProcess(Xtrains,method = 'pca',pcaComp = 22)
t2x=predict(t2x,Xtrains)

mean(t2x=t1x)
dim(t2x)
dim(t1x)
dim(t2$x[,1:10])
mean(abs(t2$x-t2x)<1e-10)
mean(abs(t2$x[,1:10]-t1x)<1e-16)
