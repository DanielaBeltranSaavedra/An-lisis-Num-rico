iteraciones = 2
tamanio = 10
x = c(0, 2, 3, 4, 6, 7, 5, 8, 9, 10, 13, 11, 12, 14, 15, 16, 18)
y = c( 1, 1, 2, 1, 1, 2, 2.5, 2.5, 4, 2.5, 2.5, 2, 1, 1, 2, 1, 1)
plot(x,y,type="l",bty="L",xlab="X",ylab="Y")
points(x,y,type="l",col="blue",lwd=5)
grid(lty=1);
copo=function(iteraciones, tamanio){
  retorno=c(1:iteraciones)
  for (i in retorno) {
    if(i==1){
      retorno[i]=tamanio
    }
    else{
      retorno[i] = (4/3)*retorno[i-1]
    }
  }
  return (retorno)
}


retorno = copo(perimetro, iteraciones)
print(valores[iteraciones])


