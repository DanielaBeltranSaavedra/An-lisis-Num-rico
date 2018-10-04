N <- 20;
Nrange <- seq(5,N,by=1)
area=rep(0,N)
fun <- -0.2
for(n in Nrange){
  x <- seq (-20,10,len=n)
  dx <- x[2]-x[1]
  f <- (1/(2*pi*n)^(1/2))*exp(((-1/2)*((x-3)^2))/(2*(n^2)))
  
  area[n] <- dx*(sum(f[0:(n-1)])+(f[1]+f[n])/2)
  print(area[n]-2)
}

txt <- sprintf("el resultado de la suma de rima es =%5.3e para n=%d rectangles \n",area[n],n)
cat(txt)
par(mar=c(1,2,2,5))
plot(Nrange, area[Nrange],type ="s",lwd=3,ylim=c(0,5), xlab="n",ylab="area")
grid(lty=1);
abline (h=2.5,col="red")
par(new =T)
plot(Nrange,err[Nrange],col="blue",axes=F,xlab=NA,ylab=NA,ylim=c(0,5))
axis(side=4);
mtext(text="error",sid=4,line=3)