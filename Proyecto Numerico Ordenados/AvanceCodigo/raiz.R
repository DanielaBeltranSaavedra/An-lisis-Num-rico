
puntofijo =function(g, x0, maxIter=100){
  print("------------TABLA DE RAICES------------")
  print("X  | VALOR      | ERROR")

  dividir=10
  hasta=1
  while(hasta<=7){
    tol=1/dividir
    dividir=dividir*dividir
    hasta=hasta+1
    k = 1
    repeat{

      x1 = g(x0)
      dx = abs(x1 - x0)
      x0 = x1
      #Imprimir estado
      cat("   x", k, " | ", x1, " |",tol,"\n")
      k = k+1

      #until
      if(dx< tol|| k > maxIter) break;
    }
    # Mensaje de salida
    if( dx > tol ){

      cat("No hubo convergencia ")
      #return(NULL)
    } else{
      cat("x* es aproximadamente ", x1, " con error menor que ", tol,"\n")
    }
  }
}


phi = function(x) 3*x^2/(1+(x/1)^2)# r=3, K=1
puntofijo(phi, 2.5, 15)



Fx <- function(x) exp(-x) + x -2
F1x <- function(x) 1-exp(-x)
# Halla la raiz de Fx
puntofalso <- function(a,b) {
  k=1
  x<-seq(a,b,0.1)
  plot(x,Fx(x),type="l",col="blue")
  abline(h=0,col="blue")
  #x<-b
  #d<-(Fx(b)*a-Fx(a)*b)/(Fx(b)-Fx(a))
  error<-1
  print("------------TABLA DE RAICES------------")
  print("X  | VALOR      | ERROR")
  while (error > 1.e-64) {
    x<-(Fx(b)*a-Fx(a)*b)/(Fx(b)-Fx(a))
    if (Fx(x) == 0) break
    if (Fx(x)*Fx(a) < 0) {b <- x}
    else {a <- x}
    error<-abs(Fx(x)/F1x(x))
    points(rbind(c(x,0)),pch=19,cex=0.7,col="red")
    cat("   x", k, " | ", x, " |",error,"\n")
    k=k+1

  }
}
puntofalso(0,3)


Fx <- function(x) exp(-x) + x -2
F1x <- function(x) 1-exp(-x)
# Metodo de la Secante
# Halla la raiz de Fx

secante <- function(x0,x1) {
  x<-(Fx(x1)*x0-Fx(x0)*x1)/(Fx(x1)-Fx(x0))
  error <-1
  k=1
  print("------------TABLA DE RAICES------------")
  print("X  | VALOR      | ERROR")
  while (error > 1.e-64) {
    x0<-x1
    x1<-x
    x<-(Fx(x1)*x0-Fx(x0)*x1)/(Fx(x1)-Fx(x0))
    if (Fx(x) == 0) break
    error<-abs(Fx(x)/F1x(x))
    cat("   x", k, " | ", x, " |",error,"\n")
    k=k+1
  }
}
secante(0,3)
