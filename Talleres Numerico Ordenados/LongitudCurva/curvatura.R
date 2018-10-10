
f = function(x){
  return(sin(x))
}

longitudCurva = function(a,b,n){
  suma = 0
  datos = seq(a,b,by=((b-a)/n)) 
  i = 1
  for(i in 1:(length(datos)-1)){
    deltaY = f(datos[i+1])-f(datos[i])
    deltaX = datos[i+1]-datos[i]
    
    calculo = sqrt(1+(deltaY/deltaX)^2)*deltaX
    suma = calculo + suma
  }
  return(c(suma,datos))
}
a = 0
b = 2*pi
n = 1000
m=10
l = longitudCurva(a,b,n)
print("la longitud de curva es ",l)

cat(l[1])
datos = l[2:length(l)]
graficar = function(datos,a,b,paso,raro){
  x = seq(a,b,by=((b-a)/n))
  y = f(x)
  x2 = seq(a,b,by=((b-a)/m))
  y2= f(x2)
 
  plot(x,y,type="l",asp=1)
  plot(x2,y2,type="l",asp=1)
  lines(datos,f(datos),col="red",lwd=5)
  lines(datos,f(datos),col="blue",lwd=1)
 
}

graficar(datos,1,-1,n,m)
fun<- function(x){
  return(sqrt(1+(cos(x))^2))
} 
#datos