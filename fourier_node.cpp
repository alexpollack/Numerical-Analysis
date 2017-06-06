#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

double fx_sample(double x)
{
    return x*x;
}

double fx_x(double x)
{
    return x;
}


//composite midpoint method
void composite_midpoint(int n,double (*fx)(double),
                        double a,double b,double& sum) {
    
    double xi=0.0,h;
    int i;
    ofstream midpoint;
    midpoint.open("MidpointFile.txt");
    std::cout << "calculating integral with a,b,n = " << a << ' ' <<
    b << ' ' << n << '\n';
    h=(b-a)/n;
    sum=0.0;
    midpoint << xi << "\t" << (h*sum)*3 << endl;
    for (i=0;i<n;i++)
    {
        xi=a+h*(i+0.5);
        sum=sum+fx(xi);
        midpoint << xi << "\t" << (h*sum)*3 << endl;
    }
    sum=h*sum;
    midpoint.close();
}

//composite simplsons method
void composite_simpsons(int n,double (*fx)(double),double a,double b,double& sum)
{
    double h;
    double* X=new double[n+1];
    
    cout << "calculating integral with a,b,n = " << a << ' ' << b << ' ' << n << '\n';
    ofstream simpsons;
    simpsons.open("simpsonsFile.txt");
    h=(b-a)/n;
    
    for (int i=0;i<n+1;i++)
    {
        X[i]=a+h*(i);
    }
    sum=fx(X[0]) + fx(X[n]);
    simpsons << X[0] << "\t" << (h)*sum<< endl;
    for (int j=1; j <n; j++)
    {
        if(j%2==0)
            sum = sum + 2.0*fx(X[j]);
        else if(j%2!=0)
            sum = sum + 4.0*fx(X[j]);
        simpsons << X[j] << "\t" << (h)*sum<< endl;
    }
    simpsons << X[n] << "\t" << (h)*sum+h<< endl;
    sum=(h/3)*sum;
    simpsons.close();
    
}


void Fourier(int n,double (*fx)(double),double a,double b,double& sum)
{
    double h;
    ofstream Fourier;
    Fourier.open("FourierFile.txt");
    cout << "calculating integral with a,b,n = " << a << ' ' << b << ' ' << n << '\n';
    int N=5;
    double* X=new double[n+1];
    double* x=new double[N];
    double* sumation=new double[N];
    double* sine=new double[N];
    
    a=-3.14159265359;
    b=3.14159265359;
    h=(b-a)/n;
    double H=(b-a)/N;
    for (int i=0;i<n+1;i++)
    {
        X[i]=a+h*(i);
    }
    for (int i=0;i<N;i++)
    {
        x[i]=a+H*(i);
    }
    int i=0;
    double sum2 =0.0;
    Fourier << a << "\t" << 0 << endl;
    for(int k=1; k<N; k++)
    {
        sum= fx(a)*sin( (k*a) ) + sin( (b*k) )*fx(b);
        
        for (int j=1; j <n; j++)
        {
            if(j%2==0)
                sum = sum + 2.0*fx(X[j])*sin( (X[j]*k) );
            else if(j%2!=0)
                sum = sum + 4.0*fx(X[j])*sin( (X[j]*k));
        }
        sum= (h)*sum*(1/b);
        sumation[i]=sum;
        i++;
    }
    
    i=0;
    for(int k=1; k<N; k++)
    {
        sine[i]=sin(k*x[i]);
        i++;
    }
    for(int k=1; k<N; k++)
    {
        Fourier << x[k] << "\t" << sumation[k]*sine[k] << endl;
    }
    
    Fourier << b << "\t" << 0 << endl;
    Fourier.close();
}


int main() {
    
    int n;
    double a,b,sum,exact;
    
    a=0.0;
    b=1.0;
    n=50;
    exact=1.0/3.0;
   
    composite_midpoint(n,fx_sample,a,b,sum);
    cout << "a= " << a << '\n';
    cout << "b= " << b << '\n';
    cout << "n= " << n << '\n';
    cout << "midpoint approximation= " << sum  << '\n';
    cout << "exact= " << exact  << '\n';
    cout << "error= " << fabs(sum-exact)  << '\n';
    
    composite_simpsons(n,fx_sample,a,b,sum);
    cout << "a= " << a << '\n';
    cout << "b= " << b << '\n';
    cout << "n= " << n << '\n';
    cout << "simpson's approximation= " << sum  << '\n';
    cout << "exact= " << exact  << '\n';
    cout << "error= " << fabs(sum-exact)  << '\n';
    a=-3.14159265359;
    b=3.14159265359;
      Fourier(n,fx_x,a,b,sum);
    cout << "a= " << a << '\n';
    cout << "b= " << b << '\n';
    cout << "n= " << n << '\n';
    cout << "bK= " << sum  << '\n';
    
    
    
    std::cout << "press <1> then <enter> in order to exit the program\n";
    int dummy_input;
    std::cin >> dummy_input;
    
}
