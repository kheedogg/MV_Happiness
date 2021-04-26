

######HW2_for Multivariate Statistics I
######201711505_통계학과_김도희

######Chapter 2. Principal Component Analysis
#Consider the  Exercise 2.8 and Exercise 2.9 in page 126.

##1 [Data 2.8.2](airpollution.txt)

setwd('C:\\Users\\kheed\\Desktop\\data')
airpollution<-read.table('airpollution.txt')
head(airpollution)

##1-(1)

install.packages('MVT')
library(MVT)
#공분산행렬사용
S<-round(cov(airpollution),3)
S
eigenS=eigen(S)
round(eigenS$values,3)
V<-round(eigenS$vectors, 3)
V
gofS<-eigenS$values/sum(eigenS$values)*100 # Goodness-of fit
round(gofS,2)
plot(eigenS$values, type='b', main='Scree Graph', xlab='Component Number', ylab='Eigenvalue')
V2<-V[,1:2]  # 주성분 2개 선택
V2
Y<-scale(airpollution, scale=F) # Centred Data Matrix
Y
P.S<-Y%*%V2  # PCs scores
P.S
plot(P.S[,1], P.S[,2], main='Plot of PCs Scores', xlab='1st PC', ylab='2nd PC')
text(P.S[,1], P.S[,2]+2, labels=rownames(P.S),cex=.6, col='blue')
abline(v=0, h=0)
#상관행렬사용
R<-round(cor(airpollution),3)
R
eigenR<-eigen(R)
round(eigenR$values, 2)  # Eigenvalues
V=round(eigenR$vectors,3)
V
gofR<-eigenR$values/sum(eigenR$values)*100  # Goodness-of fit
round(gofR, 2)
plot(eigenR$values, type='b', main='Scree Graph', xlab='Component Number', ylab='Eigenvalue')
#주성분 3개 선택
V3<-V[,1:3]
V3
Z<-scale(airpollution, scale=T)  # Standardized Data Matrix
Z
P.R<-Z%*%V3
P.R
par(mfrow=c(2,2))
plot(P.R[,1], P.R[,2], main='Plot of PCs Scores', xlab='1st PC', ylab='2nd PC')
text(P.R[,1], P.R[,2], labels=rownames(airpollution), cex=0.8, col='blue', pos=3)
abline(v=0, h=0)

plot(P.R[,2], P.R[,3], main='Plot of PCs Scores', xlab='2nd PC', ylab='3rd PC')
text(P.R[,2], P.R[,3], labels=rownames(airpollution), cex=0.8, col='blue', pos=3)
abline(v=0, h=0)

plot(P.R[,1], P.R[,3], main='Plot of PCs Scores', xlab='1st PC', ylab='3rd PC')
text(P.R[,1], P.R[,3], labels=rownames(airpollution), cex=0.8, col='blue', pos=3)
abline(v=0, h=0)





##1-(2)  Find the appropriate principal components with the goodness-of-fit and interpret them.

V3
round(gofR, 2)

##1-(3),(4) 

par(mfrow=c(1,3))
biplot(prcomp(airpollution,scale=T), choices=c(1,2))
abline(v=0, h=0)
biplot(prcomp(airpollution, scale=T), choices=c(2,3))
abline(v=0, h=0)
biplot(prcomp(airpollution, scale=T), choices=c(1,3))
abline(v=0, h=0)





##2 [Data 2.8.3](trackrecord2005-men.txt)

data<-read.table('trackrecord2005-men.txt')
library(dplyr)
glimpse(data)


##2-(1)

#install.packages('MVT')
library(MVT)
SS<-cov(data)
eigenSS=eigen(SS)
round(eigenSS$values,3)
VV<-eigenSS$vectors
gofSS<-eigenSS$values/sum(eigenSS$values)*100 # Goodness-of fit
plot(eigenSS$values, type='b', main='Scree Graph', 
xlab='Component Number', ylab='Eigenvalue')

# 주성분의 점수를 보여줄 2차원 공간을 얻기 위해 편의상 주성분 2개 선택
VV2<-VV[,1:2]
YY<-scale(data, scale=F) # Centred Data Matrix
P.SS<-YY%*%VV2  # PCs scores
P.SS
par(pty='s')  # pty=비율설정
lim<-range(pretty(P.SS))
plot(P.SS[,1],P.SS[,2], main='Plot of PCs Scores', xlab='1st PC', 
ylab='2nd PC', xlim=lim,ylim=c(-5,5))
text(P.SS[,1], P.SS[,2]+.5, labels=rownames(P.SS),cex=.6, col='blue')
abline(v=0, h=0) 
 
RR<-cor(data)
RR
eigenRR<-eigen(RR)
round(eigenRR$values, 2)  # Eigenvalues
VV=eigenRR$vectors
gofRR<-eigenRR$values/sum(eigenRR$values)*100  # Goodness-of fit
round(gofRR, 2)
plot(eigenRR$values, type='b', main='Scree Graph', 
xlab='Component Number', ylab='Eigenvalue')
VV.2<-VV[,1:2]
ZZ<-scale(data, scale=T)  # Standardized Data Matrix
P.RR<-ZZ%*%VV.2
P.RR
par(pty='s')
lim<-range(pretty(P.RR))
plot(P.RR[,1], P.RR[,2], main='Plot of PCs Scores', 
xlab='1st PC', ylab='2nd PC', xlim=lim, ylim=c(-5,5))
text(P.RR[,1], P.RR[,2], labels=rownames(data), cex=0.6, col='blue', pos=3)
abline(v=0, h=0)

##2-(2)

SS
RR
round(VV.2, 3)
round(gofRR, 3)


##2-(3)

round(gofRR, 2)


##2-(4)

sort(P.RR[,1], decreasing=T)

##2-(5)
#
par(mfrow=c(1,1))
biplot(prcomp(data,scale=T), choices=c(1,2), main='Biplot for trackrecord2005')
abline(v=0, h=0)
#
n<-nrow(data)
rownames(data)
colnames(data)
YY<-scale(data, scale=T)
svd.YY<-svd(YY)
UU<-svd.YY$u
VV<-svd.YY$v
DD<-diag(svd.YY$d)
GG<-(sqrt(n-1)*UU)[,1:2]
HH<-(sqrt(1/(n-1))*VV%*%DD)[,1:2]
CC<-rbind(GG,HH)
rownames(GG)<-rownames(data)
rownames(HH)<-colnames(data)
eig2<-(svd.YY$d)^2
per2<-eig2/sum(eig2)*100
gof2<-sum(per2[1:2])
round(per2, 2)
round(gof2, 2)
lim<-range(pretty(GG))
biplot(GG,HH, main='Biplot for trackrecord2005', xlab='1st PC', ylab='2nd PC',
       xlim=lim, ylim=lim, cex=.5, pch=12)
abline(v=0,h=0)

##2-(6)
biplot(GG,HH, main='Biplot for trackrecord2005', xlab='1st PC', ylab='2nd PC',
       xlim=lim, ylim=c(-1,1), cex=.5, pch=12)
abline(v=0,h=0)











