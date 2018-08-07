#include <iostream>
#include <math.h>
#include <cmath>
#include <cstdlib>
using namespace std;

int eval (double p[],int n,double x)
{
    float s=0;
    int i;
    for(i=0; i<n; i++)
    {
        s=s+p[i]*pow(x,n-i);
    }
    return s;
}
void metodoMuller(int n,double p[],double x0,double x1,double x2,int iteracion)
{
    double  error=1,cuadratica=1,contador=0;
    while(iteracion!=0)
    {
        double h0=0,h1=0,sigma0=0,sigma1=0,a=0,b=0,c=0,Fx1=0,Fx0=0,Fx2=0,x3=0;

        h0=x1-x0;
        h1=x2-x1;

        Fx0=eval(p,n,x0);
        Fx1=eval(p,n,x1);
        Fx2=eval(p,n,x2);

        sigma0=((Fx1-Fx0)/h0);
        sigma1=((Fx2-Fx1)/h1);
        a=((sigma1-sigma0)/(h1+h0));
        b=((a*h1)+sigma1);
        c=Fx2;
        cuadratica=(sqrt(((b+(pow(b,2)))-(4*a*c))));

        x3=(x2+(((-2)*c))/(b+cuadratica));
        error=(((x3-x2)/x3)*100);
        contador++;
        cout<<"iteracion: "<<contador<<endl;
        cout<<"El valor de x0 es: "<<x0<<endl<<"El valor de x1 es: "<<x1<<endl<<"El valor de x2 es: "<<x2<<endl<<"La raiz estimada es: "<<x3<<endl<<"El error porcentual es : "<<(error)*-1<<"%"<<endl;
        cout<<endl;
        x0=x1;
        x1=x2;
        x2=x3;
        iteracion--;
        if(error==-0.00000)
        {
            iteracion=0;
        }
    }
}
int main()
{
    int n=0,gradoNumero=0,iteracion=0;
    double p[n],numeros,x0=0,x1=0,x2=0;
    cout<<"ingrese el GRADO del POLINOMIO"<<endl;
    cin>>n;
    cout<<"ingrese el numero de ITERACIONES que desea"<<endl;
    cin>>iteracion;

    for(int i=0; i<=n; i++)
    {

        cout<<"ingrese el NUMERO del polinomio con sus signos"<<endl;
        cin>> numeros;
        cout<<"ingrese el GRADO que ACOMPAÑA al numero"<<endl;
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
    cout<<"Ingrese Xo (con punto .)"<<endl;
    cin>>x0;
    cout<<"Ingrese X1 (con punto .)"<<endl;
    cin>>x1;
    cout<<"Ingrese X2 (con punto .)"<<endl;
    cin>>x2;
    metodoMuller(n,p,x0,x1,x2,iteracion);

}
