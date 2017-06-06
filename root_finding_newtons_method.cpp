

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

// y=x^N  ln y=N ln x   y=e^{N ln x}
double realpower(double x,double N) {
    
    return exp(N*log(x));
    
}

double hsample(double x)
{
    return 3.00*x+5.00;
}

double hsampleprime(double x) {
    return 3.00;
}


double hsampleprimeprime(double x) {
    return 0.0;
}

double dist(double x,double x0,double y0) {
    double h=hsample(x);
    return (x-x0)*(x-x0)+(h-y0)*(h-y0);
}

double distprime(double x,double x0,double y0) {
    double h=hsample(x);
    return 2.0*(x-x0)+2.0*(h-y0)*hsampleprime(x);
}

double distprimeprime(double x,double x0,double y0) {
    double h=hsample(x);
    double hp=hsampleprime(x);
    double hpp=hsampleprimeprime(x);
    return 2.0+2.0*hpp*(h-y0)+2.0*hp*hp;
}

double newtons_method(double x_old,double tolerance,
                      double (*f)(double,double,double),
                      double (*fp)(double,double,double),
                      double x0,double y0 ) {
    
    int errormet,numiter;
    double f_old,f_oldprime,x_new,f_new;
    
    numiter=0;
    errormet=0;
    while (errormet==0) {
        
        f_old=f(x_old,x0,y0);
        f_oldprime=fp(x_old,x0,y0);
        x_new=x_old-f_old/f_oldprime;
        
        if (fabs(x_new-x_old)<tolerance)
            errormet=1;
        x_old=x_new;
        numiter++;
        f_new=f(x_new,x0,y0);
        std::cout << "numiter, x_new, d'(x_new) " << numiter << ' ' <<
        x_new << ' ' << f_new << '\n';
        
    }
    f_new=f(x_new,x0,y0);
    std::cout << "numiter, x_new, d'(x_new) " << numiter << ' ' <<
    x_new << ' ' << f_new << '\n';
    return x_new;
}

// (x0,y0)
// d(x)=(x-x0)^2+(h(x)-y0)^2
// d'(x)=2(x-x0)+2(h-y0)h'
// d''(x)=2+2h''(h-y0)+2h'h'
int main() {
    
    double tolerance;
    double x_old,x_new;
    
    tolerance=1.0e-10;
    double mypi=4.0*atan(1.0);
    double x0=4.0;
    double y0=1.0;
    x_old=0.0;
    std::cout << "x_old,x0,y0,dist " <<
    x_old << ' ' << x0 << ' ' << y0 << ' ' << dist(x_old,x0,y0) << '\n';
    x_new=newtons_method(x_old,tolerance,distprime,
                         distprimeprime,x0,y0);
    std::cout << "xclosest,yclosest,dist " << x_new << ' ' <<
    hsample(x_new) << ' ' << dist(x_new,x0,y0) << '\n';
    
}















/*double cube(double a);

int main()
{
    double x = 2.8319, b=3.1681, error;
    //x=b;
    error = 3.000 - x;
    
    for (int i = 0; i < 75; i++)
    {
        cout << "x=" <<x;
        double g = (7.0*x)/8.0 + 81.0/(8*cube(x));
        x = g;
        cout << "\ti=" << i << "\terror=" << 3.00 - x << endl;
    }
    return 0;
}

double cube(double a)
{
    float cu;
    cu=a*a*a;
    return(cu);
}
*/
