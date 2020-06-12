#load(TMsRAW.RData) #RAW
#load(TYs.RData)
Mlist=list(M1,M2,M3,M4,M5,M6,M7,M8,M9)
Ylist=list(T1Ys,T2Ys,T3Ys,T4Ys,T5Ys,T6Ys,T7Ys,T8Ys,T9Ys)
resmatrix=matrix(NA,nrow=9,ncol=3)

main3=function(N=22){ #N=dims
  Mlist=list(M1,M2,M3,M4,M5,M6,M7,M8,M9)
  Ylist=list(T1Ys,T2Ys,T3Ys,T4Ys,T5Ys,T6Ys,T7Ys,T8Ys,T9Ys)
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
    for (i in 1:22){ #i in 1:22. Center and scale all features.
      Xtrains[,i]=(Xtrains[,i]-mean(Xtrains[,i]))/sd(Xtrains[,i])  
      Xtests[,i]=(Xtests[,i]-mean(Xtests[,i]))/sd(Xtests[,i])  
    }
    colnames(Xtrains)=1:22
    colnames(Xtests)=1:22
    Ndims=N
#    Str=fastICA(Xtrains,n.comp=Ndims) #train
 #   StrS=Str$S                        #S_train signals
  #  SteS=(Xtests%*%Str$K)%*%Str$W     #test signals
##tryICA
    Mtrains=matrix(NA,nrow=dim(Ytrains)[1],ncol=Ndims)
    Mtests=matrix(NA,nrow=dim(Ytests)[1],ncol=Ndims)
    MP=preProcess(Xtrains,method="pca",pcaComp = Ndims)
    pcaTrainComp=predict(MP,Xtrains)
    #   MP2=preProcess(Xtests,method="pca",pcaComp = Ndims)
    pcaTestsComp=predict(MP,Xtests)

##    
    Mtrains=matrix(NA,nrow=dim(Ytrains)[1],ncol=Ndims)
    Mtests=matrix(NA,nrow=dim(Ytests)[1],ncol=Ndims)
    StrS=pcaTrainComp
    SteS=pcaTestsComp
    #StrS=Xtrains
    #SteS=Xtests
    ###FILTER
    #for(i in 1:dim(Xtrains)[1]){
    #   Xtrains[i,]=Xtrains[i,]-mean(Xtrains[i,])
    #}
    
    for(j in 1:Ndims){
      StrS[,j]=signal::filter(coefs$b, coefs$a, StrS[,j])
      SteS[,j]=signal::filter(coefs$b, coefs$a, SteS[,j])
    }
    ##CAR
    for(i in 1:dim(StrS)[1]){
      StrS[i,]=StrS[i,]-mean(StrS[i,])
    }
    for(i in 1:dim(SteS)[1]){
      SteS[i,]=SteS[i,]-mean(SteS[i,])
    }
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
    #    L=rep(0,1) #10 for now.. change to bigger eventually
    for (k in 1:B){
      i= sample(n,n,replace=T)
      S1=sda(Trainset[i,-1],as.factor(Trainset[i,1]),verbose=F)
      L[k]=mean(predict(S1,Testset[,-1])$class==Testset[,1])
      #      S1=sda(Trainset[,-1],L=as.factor(Trainset[,1]),verbose=F)
      #     L[k]=mean(predict(S1,Testset[,-1])$class==Testset[,1]) #decent!
    }
    resmatrix[inde,]=c(mean(L),mean(L)+c(-1,1)*1.96*sd(L)/(sqrt(length(L))))
  }
  return(resmatrix)
}

#   resmatrix[inde,]=c(mean(L),mean(L)+c(-1,1)*1.96*sd(L)/(sqrt(length(L))))
tepc=main3(22)
colMeans(tepc)

tepc21=main3(21)
colMeans(tepc21)

tepc14=main3(14)
colMeans(tepc14)

tete=main3()
mean(tete[,1])

tete21=main3(21) #w. pca
mean(tete21[,1])
tete21

tete222=main3(22)
tete222
colMeans(tete222)

ress=rep(0,22)
for(i in 2:22){
  tmp=main(i)
  ress[i]=mean(tmp[,1])
}
noica=main3()
noica
mean(noica[,1])

noica2=main3()
noica2
mean(noica2[,1])

noicanof=main3()
noicanof
mean(noicanof[,1])

noicanofnocar=main3()
noicanofnocar
mean(noicanofnocar[,1])

nothing2=main3()
nothing2
mean(nothing[,1])



nothing3=main3()
nothing3
mean(nothing3[,1])
mean(nothing[,1])
mean(nothing2[,1])
mean(nothing3[,1])


sda5=main()
sda5
colMeans(sda5)

sda9=main(9)
sda9
colMeans(sda9)
#invsible()
sda20a=main(20)
sda20a
colMeans(sda20a)

sda22a=main(22)
sda22a
colMeans(sda22a)

sda3=main(3)
sda3
colMeans(sda3)

sda2=main(2)
sda2
colMeans(sda2)
sda4=main(4)
sda6=main(6)
sda7=main(7)
sda8=main(8)
sda10=main(10)

sda4
mean(sda4[,1])
colMeans(sda2)

mean(sda2[,1])
mean(sda3[,1])
mean(sda4[,1])
mean(sda5[,1])
mean(sda9[,1])
mean(sda20a[,1])
mean(sda22a[,1])


#sda22=main()#ndims=22 - no dimension reduction
sda22#with artifacts
colMeans(sda22) #WOW! 36.6% - no reduction

mean(colMeans(sda20k))


B=10
n=nrow(Trainset)
beta=numeric(B)
for(b in 1:B){
  i= sample(n,n,replace=T)
  S1=sda(Trainset[i,-1],as.factor(Trainset[i,1]),verbose=F)
  beta[b]=mean(predict(S1,Testset[,-1])$class==Testset[,1])
}
beta
mean(beta)
i1= sample(n,n,replace=T)
sample(10,10,replace=F)

numeric(10)
