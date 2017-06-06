#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


using namespace std;


double hX(double x)
{
    double a= 3.0*x + 5.0;
    
    return a;
}

double hXprime(double x)
{
    double a = 3.0;
    return a;
}

double hXprimeprime(double x)
{
    double a=0.0;
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



void bisection(double &a, double &b, double &aOld, double &bOld, double x0, double y0)
{
    double c = 0;
    c = (a + b)/2.00;
    
    double distA=dist(a, x0, y0) , distB= dist(b, x0, y0), distC =dist(c, x0, y0);
    double gPrimeA = distprime(a, x0, y0), gPrimeB = distprime(b, x0, y0),gPrimeC = distprime(c, x0, y0);
    
    if( (gPrimeA*gPrimeC) < 0)
    {
        cout << "a,b,c:"<< a << ","<< b << ","<< c<< endl;
        cout << "h(a),h(b),h(c): "<<hX(a)<<","<<hX(b)<<","<<hX(c)<<endl;
        cout << "g(a),g(b),g(c):"<< distA <<","<<distB <<","<< distC<< endl;
        cout << "g(a)',g(b)',g(c)':"<<gPrimeA<<","<<gPrimeB<<","<<gPrimeC<<endl;
        bOld= b;
        b = c;
        cout << "\n[aNew, bNew]=[" << a<< ","<< b<<"]\n";
        
    }
    else if( (gPrimeB*gPrimeC) < 0)
    {
        cout << "a,b,c:"<< a << ","<< b << ","<< c<< endl;
        cout << "h(a),h(b),h(c): "<<hX(a)<<","<<hX(b)<<","<<hX(c)<<endl;
        cout << "g(a),g(b),g(c):"<< distA <<","<<distB <<","<< distC<< endl;
        cout << "g(a)',g(b)',g(c)':"<<gPrimeA<<","<<gPrimeB<<","<<gPrimeC<<endl;
        
        aOld=a;
        a = c;
        cout << "\n[aNew, bNew]=[" << a<< ","<< b<<"]\n";
        
    }
}


/*
 double distprimeprime(double x,double x0,double y0)
 {
 double h=hX(x);
 double hp=hXprime(x);
 double hpp=hXprimeprime(x);
 return 2.0+2.0*hpp*(h-y0)+2.0*hp*hp;
 }
 */

int main()
{
    double x0 =4.0, y0 = 1.0, a = -1.0, b = 9.0;
    double aOld=0, bOld = b, tolerance=1.0e-4, xClosest=0, yClosest= 0, distanceLeast;
    int i =1;
    cout << "[a, b]=[" << a<< ","<< b<<"]\n";
    
    do
    {
        bisection(a, b, aOld,bOld,x0,y0);
        cout << "interation: " << i<< endl;
        i++;
        
    }while(fabs(dist(a, x0, y0)-dist(b, x0, y0)) > tolerance );
    //while(fabs(a-b)>tolerance || a == b);
    
    double c = 0;
    c = (a + b)/2.00;
    cout << "a,b,c:"<< a << ","<< b << ","<< c<< endl;
    cout << "h(a),h(b),h(c): "<<hX(a)<<","<<hX(b)<<","<<hX(c)<<endl;
    cout << "g(a),g(b),g(c):"<< dist(a, x0, y0) <<","<<dist(b, x0, y0) <<","<< dist(c, x0, y0)<< endl;
    cout << "g(a)',g(b)',g(c)':"<<distprime(a, x0, y0)<<","<<distprime(b, x0, y0)<<",";
    cout << distprime(c, x0, y0)<<endl;
    
    if( dist(a, x0, y0) < dist(b, x0, y0) && dist(a, x0, y0) < dist(c, x0, y0) )
    {
        xClosest = a;
        yClosest = hX(a);
        distanceLeast =dist(a, x0, y0);
    }
    else if( dist(b, x0, y0) < dist(a, x0, y0) && dist(b, x0, y0) < dist(c, x0, y0) )
    {
        xClosest = b;
        yClosest = hX(b);
        distanceLeast =dist(b, x0, y0);
    }
    else if( dist(c, x0, y0) < dist(b, x0, y0) && dist(c, x0, y0) < dist(a, x0, y0) )
    {
        xClosest = c;
        yClosest = hX(c);
        distanceLeast =dist(c, x0, y0);
    }
    else
    {
        xClosest = a;
        yClosest = hX(a);
        distanceLeast =dist(a, x0, y0);
    }
    
    cout << "\nx_closest=" << xClosest << "\ty_closest=" << yClosest << "\tdistance="<< distanceLeast;
    cout << endl;
    return 0;
}


