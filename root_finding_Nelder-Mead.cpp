#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//AVERAGE VALUE OF SIN .707

using namespace std;

double hX(double x)
{
    double a = 3.0*x +5.0;
    
    return a;
}

double hXprime(double x)
{
    double a = 3.0;
    return a;
}

double hXprimeprime(double x)
{
    double a= 0.0;
    return a;
}

double dist(double x,double x0,double y0)
{
    double h=hX(x);
    return (x-x0)*(x-x0)+(h-y0)*(h-y0);
}

double distprime(double x,double x0,double y0)
{
    double h=hX(x);
    return 2.0*(x-x0)+2.0*(h-y0)*hXprime(x);
}

//step1
void order(double &x1, double &x2, double xNot, double y0)
{
    if (dist(x1, xNot, y0) <= dist(x2,xNot,y0))
    {
        x1 = x1;
        x2= x2;
    }
    else if (dist(x2, xNot, y0) <= dist(x1,xNot,y0))
    {
        double temp = x1;
        x1 = x2;
        x2= temp;
    }
}

//STEP2 centriod
double centriod(double x1, double x0, int n)
{
            x0 = (x1/n);
    return x0;
}

//STEP3 reflect
double reflect( double x2, double x0)
{
    double xR, alpha = 1.0;
    xR = x0 + alpha*(x0 - x2 );
        return xR;
}
//step3 expand
double expand(double xR, double x0)
{
    double xE, gama = 2.0;
    xE = x0 +gama*(xR - x0);
    return xE;
}

double contraction( double x0, double x2)
{
    double xC, rho =.5;
    xC = x0 + rho*(x2 - x0);
    return xC;
}
double shrink(double x1, double x2)
{
    double sigma = .5;
    x2 = x1 + sigma*(x2-x1);
    return x2;
}
            //MAIN//
int main()
{
    int n = 1, iterate=0;
    double x0=0.0, x1 = -1.0, x2 = 9.0, xR, xE,xC, xNot=4.0, y0=1.0, tolerence=10e-4;
    
    while (fabs(x1-x2) > tolerence)
    {
        order( x1, x2,  xNot, y0);
        x0=centriod(x1,x0,n);
        
            xR=reflect(x2,  x0);
            if(dist(xR, xNot, y0) <= dist(x1, xNot, y0) )
            {
                xE = expand(xR, x0);
                if (dist(xE, xNot, y0) < dist(xR, xNot, y0))
                {
                    x2 = xE;
                    xR=reflect(x2,  x0);
                    order( x1, x2,  xNot, y0);
                    x0=centriod(x1,x0,n);
                    if(dist(x1, xNot, y0) <=dist(xR, xNot, y0) )
                        x2 = xR;
                }
                
            }
            else if (dist(xR, xNot, y0) > dist(x1, xNot, y0))
            {
                x2=xR;
                order( x1, x2,  xNot, y0);
                x0=centriod(x1,x0,n);
            }
            xC= contraction( x0,x2);
            if (dist(xC, xNot, y0) < dist(xC, xNot, y0))
                x2=xC;
            else
                x2 =shrink( x1, x2);
        cout << "iterate: " << iterate << "\tx0,x1,x2:" << x0 << "," << x1 << "," << x2 << endl;
        cout << "h(x0),h(x1),h(x2):" << hX(x0) << "," << hX(x1) << "," << hX(x2) <<endl;
        cout << "g(x0),g(x1),g(x2):" << dist(x0, xNot, y0) << "," << dist(x1,xNot,y0) << "," << dist(x2,xNot,y0) <<endl;
        cout << "g'(x0),g'(x1),g'(x2):" << distprime(x0, xNot, y0) << "," << distprime(x1,xNot,y0) << "," << distprime(x2,xNot,y0) <<endl;
        cout << endl;
        iterate++;
    }
 
    return 0;
}


