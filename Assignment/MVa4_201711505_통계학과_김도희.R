

#################### MVN HW4
#################### 201711505_통계학과_김도희


#[kellogg.txt]

setwd('C:\\Users\\kheed\\Desktop\\data')
data<-read.table('kellogg.txt', header=T)
rownames<-data[,1]
data<-data[,-1]
rownames(data)<-rownames
colnames(data)<-c('칼로리','단백질','지방','나트륨','다이어트 식이섬유','
                    복합탄수화물','당분','칼륨','비타민과 무기물','유형')
View(data)

head(data)
head(scale(data))
n<-dim(data)[1]; p<-dim(data)[2]

#1)

binary<-matrix(NA,n,p)

for(i in 1:p) {
  for(k in 1:n) {
    if (data[k,i]>=colMeans(data)[i])
      binary[k,i]<-1
    else binary[k,i]<-0
  }
}

binary
colnames(binary)<-colnames(data)
rownames(binary)<-rownames
View(binary)


#2)
X<-binary
J<-matrix(1,n,p)
#유사성행렬 
Cs<-( X%*%t(X) + (J-X)%*%t(J-X) )/p
Cs
length(which(Cs<=0))

D<-matrix(NA,n,n)
for( i in 1:n ) {
  for( k in 1:n ) {
    D[i,k]<-sqrt(Cs[i,i]-2*Cs[i,k]+Cs[k,k])
  }
}
D

# 유클리드 거리
de_bi<-as.matrix(dist(X, method='euclidean'))
de_bi<-as.dist(de_bi)
de_bi



#3)

# 표준화 유클리드 거리
Z<-scale(X)
ds<-dist(Z, method='euclidean')
round(ds,3)
wards=hclust(ds,method='ward.D')
plot(wards,labels=rownames, main='Ward Linkage : Standardized Euclidean Distance')

par(mfrow=c(1,1))
# 단일연결법
single<-hclust(ds, method='single')
plot(single, hang=-1, labels=rownames, main='(a) Single Linkage')

X[c('FroF','AppJ','CorP','Nut&','Froo','Smac'),]
X[c('CorF','Spec','Cris','RiKr'),]
X[c('FruB','RaBr'),]
X[c('MuCB','Crac','NGAR','FrMW','NutW','Rais','AllB','AllF'),]
X[c('Prod'),]
X[c('JRCN','JRFN'),]

# 평균연결법
average<-hclust(ds, method='average')
plot(average, hang=-1, labels=rownames, main='(b) Average Linkage')

X[c('Prod','JRCN','JRFN'),]
X[c('FruB','RaBr','MuCB','Crac','NGAR'),]
X[c('FrMW','NutW','Rais','AllB','AllF'),]


# 와드연결법
ward<-hclust(ds, method='ward.D')
plot(ward, hang=-1,labels=rownames, main='(c) Ward Linkage')


#4)

set.seed(2017)
# K-Means Method
kmeans<-kmeans(Z,5)
cereal.name<-rownames
cluster<-data.frame(cereal.name, cluster=kmeans$cluster)
C1=cluster[(cluster[,2]==1),]
C2=cluster[(cluster[,2]==2),]
C3=cluster[(cluster[,2]==3),]
C4=cluster[(cluster[,2]==4),]
C5=cluster[(cluster[,2]==5),]
C1;C2;C3
C4;C5

# Get cluster means
View(aggregate(data, by=list(kmeans$cluster),FUN=mean))

colMeans(data)

a<-aggregate(data, by=list(kmeans$cluster),FUN=mean)
rbind(a[1,],c(NA,colMeans(data)))
rbind(a[2,],c(NA,colMeans(data)))
rbind(a[3,],c(NA,colMeans(data)))
rbind(a[4,],c(NA,colMeans(data)))
rbind(a[5,],c(NA,colMeans(data)))



# K-Medoids Method
#install.packages('cluster')
library(cluster)
kmedoids<-pam(Z, 5, metric='euclidean')
cluster<-data.frame(cereal.name, cluster=kmedoids$cluster)
C1=cluster[(cluster[,2]==1),]
C2=cluster[(cluster[,2]==2),]
C3=cluster[(cluster[,2]==3),]
C4=cluster[(cluster[,2]==4),]
C5=cluster[(cluster[,2]==5),]
C1;C2;C3
C4;C5

# Get cluster means
View(aggregate(data, by=list(kmedoids$cluster),FUN=mean))


#6)

set.seed(1234)
Z<-scale(data)
de<-as.matrix(dist(Z,method='euclidean'))
de<-as.dist(de)
de
kmean<-kmeans(data,5)
cluster<-data.frame(cereal.name, cluster=kmean$cluster)
C1=cluster[(cluster[,2]==1),]
C2=cluster[(cluster[,2]==2),]
C3=cluster[(cluster[,2]==3),]
C4=cluster[(cluster[,2]==4),]
C5=cluster[(cluster[,2]==5),]
C1;C2;C3
C4;C5


# Get cluster means
View(aggregate(data, by=list(kmean$cluster),FUN=mean))

a<-aggregate(data, by=list(kmean$cluster),FUN=mean)
rbind(a[1,],c(NA,colMeans(data)))
rbind(a[2,],c(NA,colMeans(data)))
rbind(a[3,],c(NA,colMeans(data)))
rbind(a[4,],c(NA,colMeans(data)))
rbind(a[5,],c(NA,colMeans(data)))

