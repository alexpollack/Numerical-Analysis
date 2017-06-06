#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>

using namespace std;

double fx (double xi);
double gx (double xi, const double pi);

int main()
{
    const double pi = 3.1415;
    int a = 0, N = 30;
    double xi, b = 2*pi;
    
    ofstream fxstream;
    fxstream.open("fxfile.txt");
    ofstream gxstream;
    gxstream.open("gxfile.txt");
    
    for ( int i = 0; i < N; i++)
    {
        xi = a + ((b - a)/N)*i;
        
        fxstream << xi << "\t" << fixed << showpoint;
        fxstream << setprecision(3) << fx(xi) << "\n";
    }
    for ( int i = 0; i < N; i++)
    {
        xi = a + ((b - a)/N)*i;
        
        gxstream << xi << "\t" << gx(xi, pi) << "\n";
    }
    
    fxstream.close();
    gxstream.close();
    cout << "f(x) and g(x) evaluated to N = " << N << endl;
    cout << "f(x) saved to file: fxfile.txt\n";
    cout << "g(x) saved to file: gxfile.txt\n";
}

double fx (double xi)
{
    double fx = cos(xi);
    return fx;
}

double gx (double xi, const double pi)
{
    double gx;
    
    if( xi < pi)
        return gx = -1;
    else
        return gx = 1;
}



