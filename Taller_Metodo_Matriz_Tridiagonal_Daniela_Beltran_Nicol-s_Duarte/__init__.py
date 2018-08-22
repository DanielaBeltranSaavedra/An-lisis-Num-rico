from tkinter.tix import COLUMN
print ("Ingrse la cantidad de columnas")
columnas=int(input())
aux=columnas
matriz=[]
for i in range(0,columnas):
    matriz.append([0]*columnas)
print (matriz)
for i in range(0,columnas):
    for j in range(0,columnas):
        print ("Ingrese los numeros de la matriz")
        numero=int(input())
        matriz[i][j]=numero

print(matriz)
a=[]
b=[]
c=[]
d=[]
for i in range(0,columnas):
    b.append([0]*1)
    d.append([0]*1)
for i in range(0,(columnas-1)):
    a.append([0]*1)
    c.append([0]*1)

print (a)

for i in range (0,columnas):
    print ("Ingrese los numeros de la solucion de la matriz")
    solu=int(input())
    d[i]=solu

print (d)

contb=0
contc=0
conta=0
for i in range(0,columnas):
    for j in range (0,columnas):
        if(i==j):
                auxb=matriz[i][j]
                b[contb]=int(auxb)
                contb=contb+1
        else:
            if(i==j-1):
                        auxc=matriz[i][j]
                        c[contc]=auxc
                        contc=contc+1
                        
            else:
                if(i==j+1):
                            auxa=matriz[i][j]
                            a[conta]=auxa
                            conta=conta+1
print ("aca va a")
print (a)

n=len(d)
w=[]
g=[]
x=[]
print ("ca va len")
print (n)
for i in range(0,n):
        w.append([0]*1)
        g.append([0]*1)
        x.append([0]*1)
auxb=b[0]        
w[0]=auxb
auxw=w[0]
auxd=d[0]
g[0]=int(auxd/auxw)
for i in range(1,n-1):
    aub=b[i]
    aua=a[i]
    print(aua)
    auc=c[i-1]
    print(auc)
    auw=w[i-1]
    print (auw)
    aux=aub-(aua*auc)/auw
    w=+aux
    auxds=d[i] 
    auxgs=g[i-1]
    auxito=(auxds-(aua*auxgs))/auw
    g=+auxito
print (g)
