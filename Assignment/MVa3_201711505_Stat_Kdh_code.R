

######HW3_for Multivariate Statistics I
######201711505_통계학과_김도희

######Chapter 3. Factor Analysis
#Consider the  Exercise 3.9 in page 191.

##[자료 2.8.4](protein.txt)은 유럽 25개국의 9가지 단백질 섭취원에 대한 평균 섭취량 자료가 있다.


setwd('C:\\Users\\kheed\\Desktop\\data')
X<-read.table('protein.txt',header=T)
head(X)
X<-X[,-1]
rownames(X)<-X[,1]
X<-X[,-1]
head(X)
p<-ncol(X)
n<-nrow(X)
dim(X)
View(X)

#(1) PCFA를 실시하여 스크리그림을 통하여 인자개수를 정하고 총기여율을 구하라.
Z<-scale(X, scale=T)
Z
R<-round(cor(X),3)
R
eigen.R<-eigen(R)
round(eigen.R$values,2) # Eigenvalues
V<-round(eigen.R$vectors, 2)  # Eigenvectors
V
gof<-eigen.R$values/p*100 # Goodness-of fit
round(gof, 3) # contribution rate
par(mfrow=c(1,1))
plot(eigen.R$values, type='b', main='Scree Graph', xlab='Factor Number', ylab='Eigenvalue')
sum(gof[1:2]); sum(gof[1:3]); sum(gof[1:4])  
#인자개수 4개 채택
V4<-V[,1:4]
L<-V4%*%diag(sqrt(eigen.R$values[1:4])) # Loading matrix
round(L,3)
round(diag(L%*%t(L)), 3)  # Communality
Psi<-diag(R-L%*%t(L)) # Specific Variance
round(Psi, 3)
Rm<-R-(L%*%t(L)+diag(Psi))
round(Rm, 3) # Residual matrix
#install.packages('psych')
library(psych)
pcfa<-principal(R, nfactors=4, rotate='none')
pcfa
round(pcfa$values, 2)
gof<-pcfa$values/p*100 # Goodness-of fit
round(gof,3)
sum(gof[1:4])
#총기여율 약85.82%

#(2) 인자적재값과 인자적재그림을 통하여 인자를 해석하라.
#Factor Loadings Plot : None
pcfa<-principal(Z, nfactors=4, rotate='none') # 원자료의 표준화자료와 공분산행렬 모두 동일한 결과 출력한다.
pcfa
L<-pcfa$loading[,1:4] # 인자적재행렬의 추정
Psi<-pcfa$uniquenesses # 특정분산
diag(Psi) # 특정분산행렬


Rm<-R-(L%*%t(L)+diag(Psi)) # 잔차행렬의 추정
round(Rm,3)

par(mfrow=c(3,2))
lim<-range(pretty(L))
plot(L[,1], L[,2], main='(a) Plot of Factor Loadings : None', xlab='f1', ylab='f2',
     xlim=lim, ylim=lim)
text(L[,1], L[,2], labels=rownames(L), cex=1, col='blue', pos=1)
abline(v=0,h=0)
arrows(0,0,L[,1], L[,2], col=2, code=2, length=0.1)
plot(L[,2], L[,3], main='(b) Plot of Factor Loadings : None', xlab='f2', ylab='f3',
     xlim=lim, ylim=lim)
text(L[,2], L[,3], labels=rownames(L), cex=1, col='blue', pos=1)
abline(v=0,h=0)
arrows(0,0,L[,2], L[,3], col=2, code=2, length=0.1)
plot(L[,3], L[,4], main='(c) Plot of Factor Loadings : None', xlab='f3', ylab='f4',
     xlim=lim, ylim=lim)
text(L[,3], L[,4], labels=rownames(L), cex=1, col='blue', pos=1)
abline(v=0,h=0)
arrows(0,0,L[,3], L[,4], col=2, code=2, length=0.1)
plot(L[,1], L[,3], main='(d) Plot of Factor Loadings : None', xlab='f1', ylab='f3',
     xlim=lim, ylim=lim)
text(L[,1], L[,3], labels=rownames(L), cex=1, col='blue', pos=1)
abline(v=0,h=0)
arrows(0,0,L[,1], L[,3], col=2, code=2, length=0.1)
plot(L[,2], L[,4], main='(e) Plot of Factor Loadings : None', xlab='f2', ylab='f4',
     xlim=lim, ylim=lim)
text(L[,2], L[,4], labels=rownames(L), cex=1, col='blue', pos=1)
abline(v=0,h=0)
arrows(0,0,L[,2], L[,4], col=2, code=2, length=0.1)
plot(L[,1], L[,4], main='(f) Plot of Factor Loadings : None', xlab='f1', ylab='f4',
     xlim=lim, ylim=lim)
text(L[,1], L[,4], labels=rownames(L), cex=1, col='blue', pos=1)
abline(v=0,h=0)
arrows(0,0,L[,1], L[,4], col=2, code=2, length=0.1)


#Factor Loadings Plot : Varimax
pcfa<-principal(Z, nfactors=4, rotate='varimax')
pcfa
L<-pcfa$loading[,1:4]
Psi<-pcfa$uniquenesses
Rm<-R-(L%*%t(L)+diag(Psi))
round(Rm,3)

par(mfrow=c(3,2))
lim<-range(pretty(L))
plot(L[,1], L[,2], main='(a) Plot of Factor Loadings : Varimax', xlab='f1', ylab='f2',
     xlim=lim, ylim=lim)
text(L[,1], L[,2], labels=rownames(L), cex=1, col='blue', pos=1)
abline(v=0,h=0)
arrows(0,0,L[,1], L[,2], col=2, code=2, length=0.1)
plot(L[,2], L[,3], main='(b) Plot of Factor Loadings : Varimax', xlab='f2', ylab='f3',
     xlim=lim, ylim=lim)
text(L[,2], L[,3], labels=rownames(L), cex=1, col='blue', pos=1)
abline(v=0,h=0)
arrows(0,0,L[,2], L[,3], col=2, code=2, length=0.1)
plot(L[,3], L[,4], main='(c) Plot of Factor Loadings : Varimax', xlab='f3', ylab='f4',
     xlim=lim, ylim=lim)
text(L[,3], L[,4], labels=rownames(L), cex=1, col='blue', pos=1)
abline(v=0,h=0)
arrows(0,0,L[,3], L[,4], col=2, code=2, length=0.1)
plot(L[,1], L[,3], main='(d) Plot of Factor Loadings : Varimax', xlab='f1', ylab='f3',
     xlim=lim, ylim=lim)
text(L[,1], L[,3], labels=rownames(L), cex=1, col='blue', pos=1)
abline(v=0,h=0)
arrows(0,0,L[,1], L[,3], col=2, code=2, length=0.1)
plot(L[,2], L[,4], main='(e) Plot of Factor Loadings : Varimax', xlab='f2', ylab='f4',
     xlim=lim, ylim=lim)
text(L[,2], L[,4], labels=rownames(L), cex=1, col='blue', pos=1)
abline(v=0,h=0)
arrows(0,0,L[,2], L[,4], col=2, code=2, length=0.1)
plot(L[,1], L[,4], main='(f) Plot of Factor Loadings : Varimax', xlab='f1', ylab='f4',
     xlim=lim, ylim=lim)
text(L[,1], L[,4], labels=rownames(L), cex=1, col='blue', pos=1)
abline(v=0,h=0)
arrows(0,0,L[,1], L[,4], col=2, code=2, length=0.1)


#(3) 인자점수그림을 통해 유럽 25개국 군집의 형성과 특성을 살펴보라.
fpc<-pcfa$scores
round(fpc,3)
par(mfrow=c(3,2))
par(pty='s')
lim<-range(pretty(fpc))
plot(fpc[,1], fpc[,2], main='(a) Factor Scores : f1 and f2', xlab='f1', ylab='f2',
     xlim=lim, ylim=lim)
text(fpc[,1], fpc[,2], labels=rownames(fpc), cex=1, col='blue', pos=1)
abline(v=0,h=0)
plot(fpc[,2], fpc[,3], main='(b) Factor Scores : f2 and f3', xlab='f2', ylab='f3',
     xlim=lim, ylim=lim)
text(fpc[,2], fpc[,3], labels=rownames(fpc), cex=1, col='blue', pos=1)
abline(v=0,h=0)
plot(fpc[,3], fpc[,4], main='(c) Factor Scores : f3 and f4', xlab='f3', ylab='f4',
     xlim=lim, ylim=lim)
text(fpc[,3], fpc[,4], labels=rownames(fpc), cex=1, col='blue', pos=1)
abline(v=0,h=0)
plot(fpc[,1], fpc[,3], main='(d) Factor Scores : f1 and f3', xlab='f1', ylab='f3',
     xlim=lim, ylim=lim)
text(fpc[,1], fpc[,3], labels=rownames(fpc), cex=1, col='blue', pos=1)
abline(v=0,h=0)
plot(fpc[,2], fpc[,4], main='(e) Factor Scores : f2 and f4', xlab='f2', ylab='f4',
     xlim=lim, ylim=lim)
text(fpc[,2], fpc[,4], labels=rownames(fpc), cex=1, col='blue', pos=1)
abline(v=0,h=0)
plot(fpc[,1], fpc[,4], main='(f) Factor Scores : f1 and f4', xlab='f1', ylab='f4',
     xlim=lim, ylim=lim)
text(fpc[,1], fpc[,4], labels=rownames(fpc), cex=1, col='blue', pos=1)
abline(v=0,h=0)




#(4) (1)의 인자개수에 대해 MLFA를 실시하고 (2)~(3)을 시행한 후에 결과를 서로 비교하라.
###MLFA 시행
mlfa<-factanal(Z, factors=2, rotation='none') # p=0.0709로 공통인자 수 m이 적절하다는 귀무가설 기각
mlfa<-factanal(Z, factors=3, rotation='none')
mlfa<-factanal(Z, factors=4, rotation='none')
mlfa<-factanal(Z, factors=5, rotation='none')
mlfa
##인자적재값과 인자적재그림을 통한 인자 해석
#Factor Loadings Plot : None
mlfa<-factanal(Z, factors=4, rotation='none')
mlfa
L<-mlfa$loading[,1:4] # factor loading
Psi<-mlfa$uniquenesses # specific variance
Rm<-R-(L%*%t(L)+diag(Psi))
round(Rm,3) # Residual matrix
par(mfrow=c(3,2))
lim<-range(pretty(L))
plot(L[,1], L[,2], main='(a) Plot of Factor Loadings : None', xlab='f1', ylab='f2',
     xlim=lim, ylim=lim)
text(L[,1], L[,2], labels=rownames(L), cex=1, col='blue', pos=1)
abline(v=0,h=0)
arrows(0,0,L[,1], L[,2], col=2, code=2, length=0.1)
plot(L[,2], L[,3], main='(b) Plot of Factor Loadings : None', xlab='f2', ylab='f3',
     xlim=lim, ylim=lim)
text(L[,2], L[,3], labels=rownames(L), cex=1, col='blue', pos=1)
abline(v=0,h=0)
arrows(0,0,L[,2], L[,3], col=2, code=2, length=0.1)
plot(L[,3], L[,4], main='(c) Plot of Factor Loadings : None', xlab='f3', ylab='f4',
     xlim=lim, ylim=lim)
text(L[,3], L[,4], labels=rownames(L), cex=1, col='blue', pos=1)
abline(v=0,h=0)
arrows(0,0,L[,3], L[,4], col=2, code=2, length=0.1)
plot(L[,1], L[,3], main='(d) Plot of Factor Loadings : None', xlab='f1', ylab='f3',
     xlim=lim, ylim=lim)
text(L[,1], L[,3], labels=rownames(L), cex=1, col='blue', pos=1)
abline(v=0,h=0)
arrows(0,0,L[,1], L[,3], col=2, code=2, length=0.1)
plot(L[,2], L[,4], main='(e) Plot of Factor Loadings : None', xlab='f2', ylab='f4',
     xlim=lim, ylim=lim)
text(L[,2], L[,4], labels=rownames(L), cex=1, col='blue', pos=1)
abline(v=0,h=0)
arrows(0,0,L[,2], L[,4], col=2, code=2, length=0.1)
plot(L[,1], L[,4], main='(f) Plot of Factor Loadings : None', xlab='f1', ylab='f4',
     xlim=lim, ylim=lim)
text(L[,1], L[,4], labels=rownames(L), cex=1, col='blue', pos=1)
abline(v=0,h=0)
arrows(0,0,L[,1], L[,4], col=2, code=2, length=0.1)



#Factor Loadings Plot : Varimax
mlfa<-factanal(Z, factors=4, rotation='varimax') # 원자료의 표준화자료와 공분산행렬 모두 동일한 결과 출력한다.
mlfa
L<-mlfa$loading[,1:4]
L
Psi<-mlfa$uniquenesses
Rm<-R-(L%*%t(L)+diag(Psi))
round(Rm, 3)
par(mfrow=c(3,2))
lim<-range(pretty(L))
plot(L[,1], L[,2], main='(a) Plot of Factor Loadings : Varimax', xlab='f1', ylab='f2',
     xlim=lim, ylim=lim)
text(L[,1], L[,2], labels=rownames(L), cex=1, col='blue', pos=1)
abline(v=0,h=0)
arrows(0,0,L[,1], L[,2], col=2, code=2, length=0.1)
plot(L[,2], L[,3], main='(b) Plot of Factor Loadings : Varimax', xlab='f2', ylab='f3',
     xlim=lim, ylim=lim)
text(L[,2], L[,3], labels=rownames(L), cex=1, col='blue', pos=1)
abline(v=0,h=0)
arrows(0,0,L[,2], L[,3], col=2, code=2, length=0.1)
plot(L[,3], L[,4], main='(c) Plot of Factor Loadings : Varimax', xlab='f3', ylab='f4',
     xlim=lim, ylim=lim)
text(L[,3], L[,4], labels=rownames(L), cex=1, col='blue', pos=1)
abline(v=0,h=0)
arrows(0,0,L[,3], L[,4], col=2, code=2, length=0.1)
plot(L[,1], L[,3], main='(d) Plot of Factor Loadings : Varimax', xlab='f1', ylab='f3',
     xlim=lim, ylim=lim)
text(L[,1], L[,3], labels=rownames(L), cex=1, col='blue', pos=1)
abline(v=0,h=0)
arrows(0,0,L[,1], L[,3], col=2, code=2, length=0.1)
plot(L[,2], L[,4], main='(e) Plot of Factor Loadings : Varimax', xlab='f2', ylab='f4',
     xlim=lim, ylim=lim)
text(L[,2], L[,4], labels=rownames(L), cex=1, col='blue', pos=1)
abline(v=0,h=0)
arrows(0,0,L[,2], L[,4], col=2, code=2, length=0.1)
plot(L[,1], L[,4], main='(f) Plot of Factor Loadings : Varimax', xlab='f1', ylab='f4',
     xlim=lim, ylim=lim)
text(L[,1], L[,4], labels=rownames(L), cex=1, col='blue', pos=1)
abline(v=0,h=0)
arrows(0,0,L[,1], L[,4], col=2, code=2, length=0.1)


##인자점수그림을 통해 유럽 25개국 군집의 형성과 특성살피기
mlfa<-factanal(Z, factors=4, rotation='varimax', score='regression')
mlfa
fml<-mlfa$scores
round(fml,3)
par(mfrow=c(3,2))
par(pty='s')
lim<-range(pretty(fml))
plot(fml[,1], fml[,2], main='(a) Factor Scores : f1 and f2', xlab='f1', ylab='f2',
     xlim=lim, ylim=lim)
text(fml[,1], fml[,2], labels=rownames(fpc), cex=1, col='blue', pos=1)
abline(v=0,h=0)
plot(fml[,2], fml[,3], main='(b) Factor Scores : f2 and f3', xlab='f2', ylab='f3',
     xlim=lim, ylim=lim)
text(fml[,2], fml[,3], labels=rownames(fpc), cex=1, col='blue', pos=1)
abline(v=0,h=0)
plot(fml[,3], fml[,4], main='(c) Factor Scores : f3 and f4', xlab='f3', ylab='f4',
     xlim=lim, ylim=lim)
text(fml[,3], fml[,4], labels=rownames(fpc), cex=1, col='blue', pos=1)
abline(v=0,h=0)
plot(fml[,1], fml[,3], main='(d) Factor Scores : f1 and f3', xlab='f1', ylab='f3',
     xlim=lim, ylim=lim)
text(fml[,1], fml[,3], labels=rownames(fpc), cex=1, col='blue', pos=1)
abline(v=0,h=0)
plot(fml[,2], fml[,4], main='(e) Factor Scores : f2 and f4', xlab='f2', ylab='f4',
     xlim=lim, ylim=lim)
text(fml[,2], fml[,4], labels=rownames(fpc), cex=1, col='blue', pos=1)
abline(v=0,h=0)
plot(fml[,1], fml[,4], main='(f) Factor Scores : f1 and f4', xlab='f1', ylab='f4',
     xlim=lim, ylim=lim)
text(fml[,1], fml[,4], labels=rownames(fpc), cex=1, col='blue', pos=1)
abline(v=0,h=0)



### PCFA, MLFA 결과 비교
# principal함수와 factanal 함수의 입력자료로 공분산행렬이나 상관행렬자료 대신
# 중심화 내지는 표준화 자료행렬 사용
par(mfrow=c(2,2))
plot(fml[,1], fpc[,1], main='(a) Factor Scores : ml f1 and pc f1', xlab='ml f1', ylab='pc f1',
     xlim=lim, ylim=lim)
text(fml[,1], fpc[,1], labels=rownames(fml), cex=1, col='blue', pos=1)
abline(v=0, h=0)
plot(fml[,2], fpc[,3], main='(b) Factor Scores : ml f2 and pc f2', xlab='ml f2', ylab='pc f2',
     xlim=lim, ylim=lim)
text(fml[,2], fpc[,3], labels=rownames(fml), cex=1, col='blue', pos=1)
abline(v=0, h=0)
plot(fml[,3], fpc[,3], main='(c) Factor Scores : ml f3 and pc f3', xlab='ml f3', ylab='pc f3',
     xlim=lim, ylim=lim)
text(fml[,3], fpc[,3], labels=rownames(fml), cex=1, col='blue', pos=1)
abline(v=0, h=0)
plot(fml[,4], fpc[,4], main='(d) Factor Scores : ml f4 and pc f4', xlab='ml f4', ylab='pc f4',
     xlim=lim, ylim=lim)
text(fml[,4], fpc[,4], labels=rownames(fml), cex=1, col='blue', pos=1)
abline(v=0, h=0)





#(5) 인자행렬도를 통해 단백질 섭취원 인자와 25개국 개체 간의 연관성을 살펴보라.

svd.Z<-svd(Z)
U<-svd.Z$u
V<-svd.Z$v
D<-diag(svd.Z$d)
F<-(sqrt(n-1)*U)[,1:4] # Factor Scores Matrix : F
L<-(sqrt(1/(n-1))*V%*%D)[,1:4] # Factor Loadings Matrix : Lambda
C<- rbind(F, L)
rownames(F)<-rownames(X)
rownames(L)<-colnames(X)
eig<-(svd.Z$d)^2
per<-eig/sum(eig)*100
gof<-sum(per[1:4])
per
gof

par(mfrow=c(1,2))
par(pty='s')
lim1<-range(pretty(L))
lim2<-range(pretty(F))
biplot(F[,1:2], L[,1:2], xlab='f1', ylab='f2', main='(a1) Unrotated Biplot', xlim=lim2, ylim=lim2,
       cex=1, pch=16)
abline(v=0, h=0)

varimax<-varimax(L)
Lt<-varimax$loadings
T<-varimax$rotmat
T
Ft<-F%*%T
biplot(Ft[,1:2], Lt[,1:2], xlab='f1', ylab='f2', main='(b1) Varimax Rotated Biplot',
       xlim=lim2, ylim=lim2, cex=1, pch=16)
abline(v=0, h=0)

biplot(F[,2:3], L[,2:3], xlab='f2', ylab='f3', main='(a2) Unrotated Biplot', xlim=lim2, ylim=lim2,
       cex=1, pch=16)
abline(v=0, h=0)

biplot(Ft[,2:3], Lt[,2:3], xlab='f2', ylab='f3', main='(b2) Varimax Rotated Biplot',
       xlim=lim2, ylim=lim2, cex=1, pch=16)
abline(v=0, h=0)

biplot(F[,3:4], L[,3:4], xlab='f3', ylab='f4', main='(a3) Unrotated Biplot', xlim=lim2, ylim=lim2,
       cex=1, pch=16)
abline(v=0, h=0)

biplot(Ft[,3:4], Lt[,3:4], xlab='f3', ylab='f4', main='(b3) Varimax Rotated Biplot',
       xlim=lim2, ylim=lim2, cex=1, pch=16)
abline(v=0, h=0)

biplot(F[,c(1,3)], L[,c(1,3)], xlab='f1', ylab='f3', main='(a4) Unrotated Biplot', xlim=lim2, ylim=lim2,
       cex=1, pch=16)
abline(v=0, h=0)

biplot(Ft[,c(1,3)], Lt[,c(1,3)], xlab='f1', ylab='f3', main='(b4) Varimax Rotated Biplot',
       xlim=lim2, ylim=lim2, cex=1, pch=16)
abline(v=0, h=0)

biplot(F[,c(2,4)], L[,c(2,4)], xlab='f2', ylab='f4', main='(a5) Unrotated Biplot', xlim=lim2, ylim=lim2,
       cex=1, pch=16)
abline(v=0, h=0)

biplot(Ft[,c(2,4)], Lt[,c(2,4)], xlab='f2', ylab='f4', main='(b5) Varimax Rotated Biplot',
       xlim=lim2, ylim=lim2, cex=1, pch=16)
abline(v=0, h=0)

biplot(F[,c(1,4)], L[,c(1,4)], xlab='f1', ylab='f4', main='(a6) Unrotated Biplot', xlim=lim2, ylim=lim2,
       cex=1, pch=16)
abline(v=0, h=0)

biplot(Ft[,c(1,4)], Lt[,c(1,4)], xlab='f1', ylab='f4', main='(b6) Varimax Rotated Biplot',
       xlim=lim2, ylim=lim2, cex=1, pch=16)
abline(v=0, h=0)


