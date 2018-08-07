#include <iostream>
#include <math.h>

using namespace std;

void calcularRaiz(double x,double y,int iteracion)
{

    double funciongx;
    while(iteracion!=0)
    {
        funciongx=((0.2000000)*(pow(2.7182818,-x)+1));
cout<<x<<" ";
        cout<<funciongx<<endl;
         if(funciongx==x+1){
            cout<<"lo cumple"<<endl;
         }
         else{
            cout<<"no lo cumple"<<endl;
         }
        x=funciongx;
            iteracion--;


    }

}


int main()
{
    double x=(0.0000000),y=(0.8000000);
    int iteraciones=5;
    calcularRaiz(x,y,iteraciones);
    cout << "Hello world!" << endl;
    return 0;
}
