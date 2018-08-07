#include <iostream>
#include <stdlib.h>
#include <conio.h>
#include <stdio.h>
#include <string.h>
#include<math.h>
using namespace std;
void redondeoAbsoluto(double n,double m){

double redondeo;
redondeo=1*(pow(10,(n-m)));
cout<<redondeo;

}
void errorRelativo(double m){
double error;
error=10*(pow(10,-m));
cout<<error;

}
int main()
{
    double n=0,m=0;
    string numero;
    cout<<"Ingrese la cantidad de digitos que puede almacenar"<<endl;
    cin>>m;
     cout << "Ingrese el numero" << endl;
    cin>> numero;
   n=(numero.size())-1;
cout<<"El redondeo absoluto esta dado por : ";
 redondeoAbsoluto( n, m);
cout<<endl<<"El error relativo esta dado por :";
errorRelativo(m);
    return 0;
}
