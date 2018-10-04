x=seq(0,2*pi,length=100)
plot(x,cos(x),type="l",col="blue",lwd=3)
lines(x,sin(x),col="green",lwd=3)
abline (h=0.0,col="black")
legend("bottomleft",col=c("blue","green"),legend =c("Coseno","Seno"), lwd=3, bty = "n")

areaSeno<-function(x1){
    return(sin(x1))
  
}
areaCoseno<-function(x2){
  return(cos(x2))
  
}

x1=seq(0,pi/4,length=100)
x2=seq(pi/4,pi/2,length=100)

polygon(c(x1,x2),c(areaSeno(x1),areaCoseno(x2)),border=NA,col=col2rgb("yellow",0.5))

contador=0
amplitud=90/10
arr=c()
sum=0
var=1

x=seq(0 , 10, 0.1)
ysen=c()
for(i in 1:length(x)){
  ysen[i]=sin(x[i])
}


seg <- seq.int(0, 10, length.out = 12)
segcos <- seq.int(0, 10, length.out = 12)

fx <- vector(length = length(seg))
for (i in 1:length(seg)) {
  fx[i] <- f3(seg[i])
}
fxcos <- vector(length = length(segcos))

for (i in 1:length(segcos)) {
  fxcos[i] <- f4(segcos[i])
}

for(i in 1:trapecios){
  if(i==trapecios){
    arr[i]=ysen[i*amplitud]
  }else{
    arr[i]=abs(ysen[i*amplitud])+abs(ysen[(i+1)*amplitud])
  }
  contador=((arr[i]/2)*amplitud)*0.07
  
  print (contador)
  sum=(sum+contador)
}
cat("la suma de las areas son ",sum)



