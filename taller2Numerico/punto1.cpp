#include <iostream>
#include <math.h>

using namespace std;
double horner(double p[], int n, double x)
{
    double y=p[0];
    int i,cont=0;
    for(i=1; i<n; i++)
    {
        y=x*y+p[i];
        cont=cont+2;
    }
    cout<<"Cantidad de operaciones Horner: "<<cont<<endl;
    return y;
}
int eval (double p[],int n,int x)
{
    int cont=0;
    int s=0;
    int i;
    for(i=0; i<n; i++)
    {
        s=s+p[i]*pow(x,n-i);
        cont=cont+2;
    }

    cout<<"Cantidad de operaciones evaluado: "<<cont<<endl;
    return s;
}
void derivada (double p[],int n)
{
    double aux[n-1];
    int cont=0,j=0;
    for(int i=0; i<n; i++)
    {
        if(p[i]!=0)
        {

            aux[j]=(p[i]*(n-i));
            cout<<"p"<<j<<"num"<<aux[j]<<endl;
            cont=+2;
            j++;
        }
        /*
        else
        {
            j++;
        }
        */
    }
    for(int i=0; i<(n-1); i++)
    {

        cout<<aux[i]<<"x"<<(n-1)-i<<" ";
    }
    cout<<endl<<"Hizo : "<<cont<<" Operaciones"<<endl;
}

int main()
{

    int n=0,gradoNumero=0;
    double numeros,x;
    cout<<"ingrese el grado del polinomio"<<endl;
    cin>>n;
    double p[n];
    for(int i=0; i<=n; i++)
    {

        cout<<"ingrese los numeros del polinomio con sus signos"<<endl;
        cin>> numeros;
        cout<<"ingrese el grado que acompaña al numero"<<endl;
        cin>>gradoNumero;
        if(gradoNumero==n-i)
        {
            p[i]=numeros;
        }
        else
        {
            p[i]=0;
            p[i+1]=numeros;
            i=i+1;

        }
    }

    cout<<"ingrese el numero a evaluar"<<endl;
    cin>> x;
    cout<<"Funcion evaluo  "<<(eval(p,n,x))-4<<endl;
    cout<<"Funcion Horner "<<((horner(p,n,x))*x)-4<<endl;
    cout<<"Derivada ";
    derivada(p,n);
}
