#include <iostream>
#include <math.h>
//Codido por: Daniela Beltran y Nicolas Duarte
using namespace std;

double calcularFormula( double gravedad, double masa,double velocidad,double tiempo,double coeficienteArrastre)
{
    double Formula=(((gravedad*masa)/coeficienteArrastre)*(1-exp(-(coeficienteArrastre/masa)*tiempo))-velocidad);
    cout<<"el valor correcto de la formula es : "<< Formula<<endl;
    cout<<"el valor correcto del coeficiente de arrastre es "<<coeficienteArrastre;
    return Formula;
}
double triseccion (double formulaCorrecta, int particiones, double limiteInicio, double limiteFin,double coeficienteArrastre)
{
    int iteraciones=0;
 double limiteMedioInicio=-1;
        double limiteMedioFin=-1;
    while (limiteInicio!=coeficienteArrastre&&limiteFin!=coeficienteArrastre&&limiteMedioFin!=coeficienteArrastre&&limiteMedioInicio!=coeficienteArrastre)
    {
        //cantidadEspacioLimites es la resta de los limites de inicio y fin para saber cual es la diferencia entre ellos

        double cantidadEspacioLimites=limiteFin-limiteInicio;
        // despues en distanciaParticion es el valor que incrementa los limites segun las particiones requerisas entonces divivimos la cantidad de estapacios de loss limites que es 4 dividida en particiones que es 3
        double distanciaParticion=cantidadEspacioLimites/particiones;
         limiteMedioInicio=limiteInicio+distanciaParticion;
         limiteMedioFin=limiteMedioInicio+distanciaParticion;
        cout<<endl<<particiones;
        iteraciones++;
        cout<<"iteracion numero : " <<iteraciones<<" Resultados : "<<limiteInicio<<" - "<<limiteMedioInicio<<" - "<<limiteMedioFin<<" - "<<limiteFin<<endl;
        double valorMenorDeLaIteracion=0,valorSegundoMenorIteracion=0;

        if(limiteInicio==coeficienteArrastre||limiteFin==coeficienteArrastre||limiteMedioFin==coeficienteArrastre||limiteMedioInicio==coeficienteArrastre)
        {
            //encontre el valor

            cout<<"el valor es el mismo";
return 0;
        }


        if(limiteInicio<coeficienteArrastre)
        {
            if(limiteMedioInicio<coeficienteArrastre)
            {
                //olvidarLimiteInicio
                if(limiteMedioFin<coeficienteArrastre)
                {
                    limiteInicio=limiteMedioFin;
                    if(limiteFin>coeficienteArrastre)
                    {
                        limiteFin=limiteFin;
                    }
                }
            }

        }
        if(limiteMedioInicio>coeficienteArrastre){
            if(limiteInicio<coeficienteArrastre){
                limiteInicio=limiteInicio;
                limiteFin=limiteMedioFin;
            }

        }
        if(limiteMedioInicio<coeficienteArrastre){
            if(limiteMedioFin>coeficienteArrastre){
                limiteInicio=limiteMedioInicio;
                limiteFin=limiteMedioFin;
            }
        }







    }





}
double cuatriseccion (double formulaCorrecta, int particiones, double limiteInicio, double limiteFin,double coeficienteArrastre)
{
    int iteraciones=0;
 double limiteMedioInicio=-1;
 double limiteIntermedio=-1;
        double limiteMedioFin=-1;
    while (limiteInicio!=coeficienteArrastre)
    {

        double cantidadEspacioLimites=limiteFin-limiteInicio;
        double distanciaParticion=cantidadEspacioLimites/particiones;
         limiteMedioInicio=limiteInicio+distanciaParticion;
         limiteIntermedio=limiteMedioInicio+distanciaParticion;
         limiteMedioFin=limiteIntermedio+distanciaParticion;
        cout<<endl<<particiones;
        iteraciones++;
        cout<<"iteracion numero : " <<iteraciones<<" Resultados : "<<limiteInicio<<" - "<<limiteMedioInicio<<" - "<<limiteIntermedio<<" - "<<limiteMedioFin<<" - "<<limiteFin<<endl;

        if(limiteInicio==coeficienteArrastre)
        {
            //encontre el valor

            cout<<"el valor es el mismo";

        }
        if(limiteFin==coeficienteArrastre)
        {
            //encontre el valor

            cout<<"el valor es el mismo";

        }
        if(limiteMedioFin==coeficienteArrastre)
        {
            //encontre el valor

            cout<<"el valor es el mismo";

        }
        if(limiteMedioInicio==coeficienteArrastre)
        {
            //encontre el valor

            cout<<"el valor es el mismo";

        }
        if(limiteIntermedio==coeficienteArrastre){
             cout<<"el valor es el mismo";
        }

        if(limiteInicio<coeficienteArrastre)
        {
            if(limiteMedioInicio<coeficienteArrastre)
            {
                //olvidarLimiteInicio
                if(limiteIntermedio<coeficienteArrastre){
                if(limiteMedioFin<coeficienteArrastre)
                {
                    limiteInicio=limiteMedioFin;
                    if(limiteFin>=coeficienteArrastre)
                    {
                        limiteFin=limiteFin;
                    }
                }
                }
            }

        }
        if(limiteInicio<coeficienteArrastre){
            if(limiteMedioInicio>coeficienteArrastre){
                limiteInicio=limiteInicio;
                limiteFin=limiteMedioInicio;
            }
        }
    if(limiteMedioInicio<coeficienteArrastre){
        if(limiteFin>coeficienteArrastre){
            limiteInicio=limiteMedioInicio;
            limiteFin=limiteFin;
        }
    }
    if(limiteMedioInicio<coeficienteArrastre){
        if(limiteMedioFin>coeficienteArrastre){
            limiteInicio=limiteMedioInicio;
            limiteFin=limiteMedioFin;
        }
    }

    if(limiteMedioInicio<coeficienteArrastre){
        if(limiteIntermedio>coeficienteArrastre){
            limiteInicio=limiteMedioInicio;
            limiteFin=limiteIntermedio;
        }
    }
    if(limiteIntermedio<coeficienteArrastre){
        if(limiteMedioFin>coeficienteArrastre){
            limiteInicio=limiteIntermedio;
            limiteFin=limiteMedioFin;
        }
    }

    }





}
int main()
{
    double gravedad=9.8, masa=68.1, velocidad=40, tiempo=10,coeficienteArrastre=14.7802;
    int particiones,limiteInicio=12,limiteFin=16;
    double formulaCorrecta;
    cout<<"Ingrese el numero de particiones que desea 3 o 4 "<<endl;
    cin>>particiones;
    formulaCorrecta=calcularFormula(gravedad,masa,velocidad,tiempo,coeficienteArrastre);
    if(particiones==3){
         triseccion ( formulaCorrecta, particiones,  limiteInicio,  limiteFin,coeficienteArrastre);
    }
   else{
        if(particiones==4)
      cuatriseccion(formulaCorrecta, particiones,  limiteInicio,  limiteFin,coeficienteArrastre);
   }


}
