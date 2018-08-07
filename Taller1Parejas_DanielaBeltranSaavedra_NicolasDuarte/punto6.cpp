#include <iostream>

using namespace std;
void DivisionesAritmeticas (int algoritmo){
    int d;
    while(algoritmo > 0){
        d=(algoritmo%2);
        algoritmo=(algoritmo/2);
        cout<<d<<endl;
    }



}
int main()
{
    int algoritmo;
    cout << "Ingrese el tamaño del algoritmo" << endl;
    cin>> algoritmo;
    DivisionesAritmeticas(algoritmo);

}
/* """la cantidad de divisiones aritmeticas T(n) es integer(log2 n)+1 T(n) expresada en O() es de O(log2(n))"""*/
