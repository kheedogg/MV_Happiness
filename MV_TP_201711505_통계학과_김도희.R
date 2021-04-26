
######################## Term Project_MVN
######################## 201711505_통계학과_김도희


setwd('C:\\Users\\kheed\\Desktop')

####install.packages

install.packages('readxl')
install.packages('dplyr')
install.packages('MVN')
install.packages('psych')
install.packages('biotools')
install.packages('cluster')


### 데이터 불러오기

library(readxl)
hap<-read_excel('Chapter2OnlineData.xls',sheet=2)
head(hap)
hap<-as.matrix(hap)
country<-hap[,1]
rownames(hap)<-country
hap<-hap[,-1]
head(hap)
View(hap)
dim(hap)
n=dim(hap)[1]
p=dim(hap)[2]
colnames(hap)
X<-matrix(NA, n,p )
for ( i in 1:p ) {
  X[,i]<-as.numeric(hap[,i])
}


colnames(X)<-colnames(hap)
rownames(X)<-rownames(hap)
X<-as.data.frame(X)

library(dplyr)
glimpse(X)
str(X)
plot(X)
colnames(X)[5:10]<-substr(colnames(X)[5:10],15,100)
X<-X[,-c(2,3)]
n=dim(X)[1]
p=dim(X)[2]

##### 3-3) 정규성 검토

xbar=colMeans(X)
S=cov(X)
m=mahalanobis(X, xbar, S)
m=sort(m)
id=seq(1, n)
pt=(id-0.5)/n
q=qchisq(pt, p)
plot(q, m, pch="*", xlab="Quantile", ylab="Ordered Squared Distance", main='Chi-Sqaure Q-Q Plot')
abline(0, 1)

rq=cor(cbind(q, m))[1,2]
rq


##### 4-1) PCA

# S 선택
S<-cov(X)
eigen.S=eigen(S)
round(eigen.S$values, 3) # Eigenvalus
V=round(eigen.S$vectors, 3) # Eigenvaectors
V
gof=eigen.S$values/sum(eigen.S$values)*100 # Goodness-of fit
round(gof, 2)

# R 선택
R=round(cor(X),3)
R
eigen.R=eigen(R)
round(eigen.R$values, 2) # Eigenvalues
V=round(eigen.R$vectors, 2) # Eigenvectors
gof=eigen.R$values/sum(eigen.R$values)*100 # Goodness-of fit
round(gof, 2)
plot(eigen.R$values, type="b", main="Scree Graph", xlab="Component Number", ylab="Eigenvalue")
V3=V[,1:3]
rownames(V3)<-colnames(X)
colnames(V3)<-c('P1','P2','P3')
V3
Z=scale(X, scale=T) # Standardized Data Matrix
Z
P=Z%*%V3            # PCs Scores
round(P, 3)
par(mfrow=c(1,3))
plot(P[,1], P[, 2], main="Plot of PCs Scores", xlab="1st PC", ylab="2nd PC")
text(P[,1], P[, 2], labels=rownames(X), cex=0.8, col="blue", pos=3)
abline(v=0, h=0)
plot(P[,2], P[, 3], main="Plot of PCs Scores", xlab="2nd PC", ylab="3rd PC")
text(P[,2], P[, 3], labels=rownames(X), cex=0.8, col="blue", pos=3)
abline(v=0, h=0)
plot(P[,1], P[, 3], main="Plot of PCs Scores", xlab="1st PC", ylab="3rd PC")
text(P[,1], P[, 3], labels=rownames(X), cex=0.8, col="blue", pos=3)
abline(v=0, h=0)


# Biplot
Y <- scale(X,scale=T)
svd.Y <- svd(Y) 
U <- svd.Y$u    
V <- svd.Y$v 
D <- diag(svd.Y$d)
G <- (sqrt(n-1)*U)[,1:3]
H <- (sqrt(1/(n-1))*V%*%D)[,1:3]
rownames(G)<-rownames(X)
rownames(H)<-colnames(X) 
par(mfrow=c(1,1))
lim<-range(pretty(G))
biplot(G[,1:2],H[,1:2], xlab="1st PC(47.92%)",ylab="2nd PC(17.87%)", main="Biplot for Happiness Scores Data",
       xlim=lim,ylim=lim,cex=0.8,pch=16)
abline(v=0,h=0)
biplot(G[,2:3],H[,2:3], xlab="2nd PC(47.92%)",ylab="3rd PC(14.53%)", main="Biplot for Happiness Scores Data",
       xlim=lim,ylim=lim,cex=0.8,pch=16)
abline(v=0,h=0)
biplot(G[,c(1,3)],H[,c(1,3)], xlab="1st PC(47.92%)",ylab="3rd PC(14.53%)", main="Biplot for Happiness Scores Data",
       xlim=lim,ylim=lim,cex=0.8,pch=16)
abline(v=0,h=0)





##### 4-2) FA

###PCFA
library(psych)
R<-cor(X)
pcfa<-principal(R, nfactors=3, rotate="varimax")
pcfa$loadings
Psi=pcfa$uniquenesses
Rm = R-(L%*%t(L) + diag(Psi))
colnames(Rm)<-c('X1','X2','X3','X4','X5','X6','X7','X8')
round(Rm, 3)


###MLFA
library(psych)
mlfa<-factanal(covmat=R, factors = 3,  rotation="varimax") # rotation="none" 
mlfa$loadings
Psi=mlfa$uniquenesses # specific variance
Rm = R-(L%*%t(L) + diag(Psi)) 
colnames(Rm)<-c('X1','X2','X3','X4','X5','X6','X7','X8')
round(Rm, 3)





## 결과비교
par(mfrow=c(1,3))
fpc<-pcfa$loadings[,1:3]
fml<-mlfa$loadings[,1:3]
lim<-range(pretty(fml))
plot(fml[,1], fpc[,1],main="(a) Factor Scores : ml f1 & pc f1",  xlab="f1", ylab="f1",
     xlim=lim, ylim=lim)
text(fml[,1], fpc[,1], labels=rownames(L), cex=0.8, col="blue", pos=1)
abline(v=0, h=0)
plot(fml[,2], fpc[,2],main="(b) Factor Scores : ml f2 & pc f2",  xlab="f2", ylab="f2",
     xlim=lim, ylim=lim)
text(fml[,2], fpc[,2], labels=rownames(L), cex=0.8, col="blue", pos=1)
abline(v=0, h=0)
plot(fml[,3], fpc[,3],main="(b) Factor Scores : ml f3 & pc f3",  xlab="f3", ylab="f3",
     xlim=lim, ylim=lim)
text(fml[,3], fpc[,3], labels=rownames(L), cex=0.8, col="blue", pos=1)
abline(v=0, h=0)


##### 4-3) CA



#### Hierarchical CA
Z<-scale(X)
rownames(Z)<-seq(1:156)
ds <- dist(Z, method="euclidean")
round(ds, 3)
#단일연결법
single=hclust(ds, method="single")
plot(single, hang=-1, main="(a) Sinle Linkage", cex=0.6)

#완전연결법
complete=hclust(ds, method="complete")
plot(complete,hang=-1, main="(b) Complete Linkage", cex=0.6)

#평균연결법
average=hclust(ds, method="average")
plot(average, hang=-1, main="(c) Average Linkage", cex=0.6)

#와드연결법
ward=hclust(ds, method="ward.D")
plot(ward, hang=-1, main="(d) Ward Linkage", cex=0.6)



#### Non-Hierarchical CA

kmeans <- kmeans(Z, 4) # 4 cluster solution
cluster=data.frame(rownames(X),cluster=kmeans$cluster)
C1=cluster[(cluster[,2]==1),]
C2=cluster[(cluster[,2]==2),]
C3=cluster[(cluster[,2]==3),]
C4=cluster[(cluster[,2]==4),]
C1;C2;C3;C4

aggregate(X, by=list(kmeans$cluster),FUN=mean)




library(cluster)
kmedoids<-pam(Z, 4, metric='euclidean')
cluster<-data.frame(rownames(X), cluster=kmedoids$cluster)
C1=cluster[(cluster[,2]==1),]
C2=cluster[(cluster[,2]==2),]
C3=cluster[(cluster[,2]==3),]
C4=cluster[(cluster[,2]==4),]
C1;C2;C3;C4







