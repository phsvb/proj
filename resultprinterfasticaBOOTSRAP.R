load('proj/TMsRAW.RData') #RAW
load('proj/TYs.RData') #RAW
load('proj/EMsRAW.RData') #RAW
load('proj/EYs.RData') #RAW
#load(TYs.RData)
Mlist=list(M1,M2,M3,M4,M5,M6,M7,M8,M9)
Ylist=list(T1Ys,T2Ys,T3Ys,T4Ys,T5Ys,T6Ys,T7Ys,T8Ys,T9Ys)
resmatrix=matrix(NA,nrow=9,ncol=3)
dim(M1)
dim(N1)
dim(E1Ys)

main=function(N=22){ #N=dims
 # Mlist=list(M1,M2,M3,M4,M5,M6,M7,M8,M9)
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
#    for(i in 1:dim(Xtrains)[1]){
 #       Xtrains[i,]=Xtrains[i,]-mean(Xtrains[i,])
#    }##
#    #CAR
#    for(i in 1:dim(Xtests)[1]){
 #     Xtests[i,]=Xtests[i,]-mean(Xtests[i,])
 #   }##
    for (i in 1:22){ #i in 1:22. Center and scale all features.
      Xtrains[,i]=(Xtrains[,i]-mean(Xtrains[,i]))/sd(Xtrains[,i])  
      Xtests[,i]=(Xtests[,i]-mean(Xtests[,i]))/sd(Xtests[,i])  
    }
    #CAR
    #   for(i in 1:dim(Xtrains)[1]){
    #      Xtrains[i,]=Xtrains[i,]-mean(Xtrains[i,])
    #    }##
    #CAR
    #   for(i in 1:dim(Xtests)[1]){
    #      Xtests[i,]=Xtests[i,]-mean(Xtests[i,])
    #    }##
    
    ###filter
    
    #for(j in 1:22){
    #        Xtrains[,j]=signal::filter(coefs$b, coefs$a, Xtrains[,j])
    #        Xtests[,j]=signal::filter(coefs$b, coefs$a, Xtests[,j])
    #}
      
      
    colnames(Xtrains)=1:22
    colnames(Xtests)=1:22
    Ndims=N
    Str=fastICA(Xtrains,n.comp=Ndims) #train
    StrS=Str$S                        #S_train signals
    SteS=(Xtests%*%Str$K)%*%Str$W     #test signals
    Mtrains=matrix(NA,nrow=dim(Ytrains)[1],ncol=Ndims)
    Mtests=matrix(NA,nrow=dim(Ytests)[1],ncol=Ndims)
    ###FILTER
    #  for(i in 1:dim(Xtrains)[1]){
    #     Xtrains[i,]=Xtrains[i,]-mean(Xtrains[i,])
    #    }
    
    for(j in 1:Ndims){
      StrS[,j]=signal::filter(coefs$b, coefs$a, StrS[,j])
      SteS[,j]=signal::filter(coefs$b, coefs$a, SteS[,j])
    }
    ##CAR
    #for(i in 1:dim(StrS)[1]){
    #  StrS[i,]=StrS[i,]-mean(StrS[i,])
    #}
    #for(i in 1:dim(SteS)[1]){
    #  SteS[i,]=SteS[i,]-mean(SteS[i,])
    #}
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

test14=main(14)
colMeans(test14)
test22=main(22)
colMeans(test22)

all14f=main(14)#no car
colMeans(all14f) # using all data 35%. drop of 1.6 percent

toda142=main(14)#no car
colMeans(toda142) # 36.6

toda14=main(14)#no car
colMeans(toda14)

#   resmatrix[inde,]=c(mean(L),mean(L)+c(-1,1)*1.96*sd(L)/(sqrt(length(L))))
tod14c=main(14)
colMeans(tod14c) #no car at all

tod14cc=main(14) #car and center early
colMeans(tod14cc)

colMeans(tod21CAR)
tod14CAR=main(14)
tod21CAR=main(21)
colMeans(tod14CAR)
colMeans(tod21CAR)

tod14NCNS=main(14) #no car
colMeans(tod14NCNS)


tod14NCfaf=main(14)
colMeans(tod14NCfaf)
tod14=main(14)
tod22=main(22)

tod14NC=main(14) #no car
tod22NC=main(22) #no car

colMeans(tod14)
colMeans(tod22)

colMeans(tod14NC)
colMeans(tod22NC)
##
##
##
ress=rep(0,21)
save(ress,file='ficancomp.RData')
for(i in 2:22){
  tmp=main(i)
  ress[i-1]=mean(tmp[,1])
}

qplot(1:21,ress,geom='line')+
  labs(x='amount of comps',y='test-accuracy')#+
#  title('fastICA test-accuracy given amount of compositions')

#######
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
