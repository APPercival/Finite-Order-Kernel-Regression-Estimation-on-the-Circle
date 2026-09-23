#Install if necessary
#install.package("foreach")
#install.package("doParallel")

#Load packages necessary for parallel computation
library(foreach)
library(doParallel)

#Load your formatted data; angles in column 1, response variable in column 2
data=read.csv("sampleData.csv")

#Define the ceiling and floor functions
newCeiling<-function(u){if(u%%1==0){u+1}else{ceiling(u)}}
newFloor<-function(u){if(u%%1==0){u-1}else{floor(u)}}

n<-length(YData)


num_cores <- detectCores() - 1 # Leave one core free so as to not overload machine
cl <- makeCluster(num_cores)
registerDoParallel(cl)

sCandidates = seq(1.5, 5.5, 0.05) #Adjust this range as needed. These are the values of s that will be evaluated in the CV(s) function 

#perform the cross-validation step in parallel
CV <- foreach(k = 1:length(sCandidates), .combine = 'c') %dopar% {
  
  s = sCandidates[k]
  r = newCeiling(s) + 3
  Trunc = newFloor(((1 / (pi * (r - 1))) * (n - 1)^((s + r) / (2 * s + 1)))^(1 / (r - 1))) + 1
  g = function(u) { 1 / (1 + abs(u)^r) }
  h = (n - 1)^(-1 / (2 * s + 1))
  nu_seq = 1:Trunc
  g_terms <- g(h * nu_seq)
  
  CVvec = rep(0, n)
  topCV = rep(0, n)
  bottomCV = rep(0, n)
  
  for(i in 1:n) {
    diff_angles <- angleData[i] - angleData[-i]
    YData_subset <- YData[-i]
    temp_vec1 <- drop(g_terms %*% cos(outer(nu_seq, diff_angles)))
    temp_vec2 <- 1 + 2 * temp_vec1
    bottomCV[i] <- sum(temp_vec2)
    topCV[i] <- sum(temp_vec2 * YData_subset)
    CVvec[i] = (YData[i] - topCV[i] / bottomCV[i])^2
  }
  
  mean(CVvec) 
}

stopCluster(cl)

s=sCandidates[which.min(CV)]#Select the s which minimised the CV function
r=newCeiling(s)+3
Trunc=newFloor(((1/(pi* (r - 1))) *(n - 1)^((s + r)/(2* s + 1)))^(1/(  r - 1)))+1
g=function(u){1/(1+abs(u)^r)}
h=(n-1)^(-1/(2*s+1))
nu_seq=1:Trunc
g_terms<-g(h*nu_seq)


theta_vec<-seq(-pi,pi,.01) #This is the resolution which you want to plot the KRE. Adjust as desired

#Plot the KRE over the range
plot_vec=matrix(rep(0,n*length(theta_vec)),c(n,length(theta_vec)))
for(i in 1:n){
  plot_vec[i,]=angleData[i]-theta_vec
}
plot_mat<-array(sapply(1:Trunc,function(nu){cos(nu*plot_vec)}),c(n,length(theta_vec),Trunc))
temp=sweep(plot_mat,3,g_terms,'*')
mat2=1+2*apply(temp,1:2,sum)
bottom<-apply(mat2,2,sum)
mat3<-sapply(1:n,function(i){YData[i]*mat2[i,]})
top<-apply(mat3,1,sum)
angleData_plot=angleData
for(i in 1:n){if(angleData_plot[i]>pi){angleData_plot[i]<-angleData_plot[i]-2*pi}}

#Adjust the margins
par(mar=rep(5,4))  
plot(angleData_plot,YData,cex=0.5,pch=19,xlab=expression(theta),ylab=expression(hat(f)[s]),main=paste("KRE plot with s=",s), col="gray",xaxt="n",yaxt="n",xlim=c(-pi,pi),xaxs="i",mar=rep(3,4))
axis(1,at=seq(-pi,pi,pi/4),labels=c(paste("-\U03C0"),paste("-3\U03C0/4"),paste("-\U03C0/2"),paste("-\U03C0/4"),paste("0"),paste("\U03C0/4"),paste("\U03C0/2"),paste("3\U03C0/4"),paste("\U03C0")))
axis(2,at=seq(-6,6,2))
lines(theta_vec,top/bottom,col="black",lwd=4)


