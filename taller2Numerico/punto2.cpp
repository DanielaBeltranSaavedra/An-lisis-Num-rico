#include <iostream>
#include <math.h>
#include <cstdlib>
using namespace std;

float f(float x);
float df(float x);
float bisec(float a, float b);
int main(){
	float x=0,x1=0,e=100;
/** resuelto por newton*/
	while(e>1e-4){
		x=x1;
		x1=x-f(x)/df(x);
		e=fabs(x-x1);
	}
	cout<<x1<<endl;
/* resuleto por biseccion*/
	float a=0,b=12;
	cout<<bisec(a,b)<<endl;

}

float f(float x){
	return (32-2*x)*(24-2*x)*x-1000;
}

float df(float x){
	float h=1e-6;
	return (f(x+h)-f(x))/h;
}


float bisec(float a,float b){
	float x=(a+b)/2,e=fabs((a-b)/2);
	while(e>1e-4){
		if(f(x)==0){return x;}
		if(f(x)*f(a)<0){b=x;}else{a=x;}
		e=fabs((a-b)/2);
		x=(a+b)/2;
	}
	return x;
}
