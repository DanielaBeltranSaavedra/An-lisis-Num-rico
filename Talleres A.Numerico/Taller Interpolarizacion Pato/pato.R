
x=c(0   ,2 ,4   ,6 ,7   ,10,12.5,15  ,20,21,22  ,23,23.5,24,              0 ,0.3,0.8,2,3.3,5,7,8.5,9.5                ,9.5,9.3,9.4,9.4,9.3,9.1,8.5,   10.3,11,11.5,11.7 ,11.9,11.9,11.9              ,13,14,15,16,17,18,20             ,21,22,22.5,23,24)
y=c(11.5,11.9,13.5,13.9,12.3,13,13.5  ,13.3,11,10.5,10.3,9.9,9.7,9.5 ,    11.5,10.8,10.3,10,10.3,10.2,10.5,9.5,9        ,7,5,4.3,4,2.7,2.2,1,          1.3,2 ,2.7 ,3.5,4,6,7                         , 7,6.5,6.3,6.3,6.5,7,8            ,7.5,8,8.5,9,9.5)


plot(x,y, pch=19, cex=1, col = "red", asp=1,xlab="X", ylab="Y", main="Diagrama ")
cubicSpline <- function(x,y) {
  a = rep(y)
  n = length(x)
  h <- (c(x,0) - c(0,x))[2:n]
  alph <- (3/c(1,h,1,1)*(c(a,1,1) - c(1,a,1)) - 3/c(1,1,h,1)*(c(1,a,1)-c(1,1,a)))[3:n]
  A <- c(1,rep(0,times=n-1))
  for (i in 1:(n-2)) {
    A <- rbind(A,c( rep(0,times=i-1) , c(h[i],2*(h[i]+h[i+1]),h[i+1]) , rep(0,times=n-i-2) ) )
  }
  A <- rbind(A,c(rep(0,times=n-1),1))
  b <- c(0,alph,0)
  c <- solve(A, b)
  b <- ((c(a,0) - c(0,a))/c(1,h,1) - c(1,h,1)/3*(c(c,0) + 2*c(0,c)))[2:n]
  d <- ((c(c,0) - c(0,c))/(3*c(1,h,1)))[2:n]
  ans = rbind(a[1:n-1],b,c[1:n-1],d)
}
draw <- function(x,y) {
  t = 1:length(x)
  sx = cubicSpline(t,x)
  sy = cubicSpline(t,y)
  for (i in 1:(length(t)-1)) {
    dat<- data.frame(t=seq(t[i],t[i+1], by=0.1) )
    fx <- function(x) (sx[1,i] + sx[2,i]*(x-t[i]) + sx[3,i]*(x-t[i])^2 + sx[4,i]*(x-t[i]))
    fy <- function(x) (sy[1,i] + sy[2,i]*(x-t[i]) + sy[3,i]*(x-t[i])^2 + sy[4,i]*(x-t[i]))
    dat$y=fy(dat$t)
    dat$x=fx(dat$t)
   ## points(dat$x,dat$y,type='l', col='blue')
  }
}
draw(x,y)
