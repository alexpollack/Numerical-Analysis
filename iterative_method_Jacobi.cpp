#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <iomanip>

using namespace std;

double dotProduct( double x[], double y[], int n);
void addVector( double x[], double y[], int n, double *z);
void ADotX(const double x[], const double d[],const double u[],const double l[], int n, double *z);
void resid(const double b[],const double d[],const double u[],const double l[],int n,const double z[], double *r);
void dInvR( const double d[], const double r[], int n, double *dx);
void InverseD( const double d[], int n, double *InvD);
void jacobiX(double *x,double *jacX, const double d[], const double b[], const double u[], const double l[], int n);
double rnorm( const double r[], int n);


int main()
{
    int n=3, k=0,i;
    double x[n], jacX[n],b[n],d[n],u[n-1],l[n-1], z[n], r[n],InvD[n]/*,dx[n-1]*/;
    double epsilon=1.0, alpha=2.0, beta=0.2,rho[n],h,Ri;
    double error = .0000001, sum;
    
    cout << "Jacobi Method\n";
    h= (1.0/n);
    //rho
    i=1;
    for(int j =0; j<=n; j++)
    {
        Ri= i * h;
        if( i==1 )
            rho[j] = epsilon;
        else if( i > 1 && Ri <= .5 )
            rho[j] = alpha;
        else if( Ri > .5 && i<(n+1) )
            rho[j] = beta;
        else if( i == (n+1) )
            rho[j] = epsilon;
        i++;
    }
    //b
    i=1;
    for(int j=0; j<n; j++)
    {
        Ri= ( (i*1.0)/n );
        b[j] = (h*h) * ( (1.0-Ri)*(1.0-Ri) ) * (Ri*Ri);
        i++;
    }
    //diag
    for(int i =0; i<n;i++)
        d[i] = ( rho[i] + rho[i+1] );
    //intial x=0 guess, k=0
    for(int i=0; i<n;i++)
    {
        x[i]=0;
    }
    cout << "\nInitial Guess X0:\n";
    for(int i=0;i<n;i++)
    {
        cout << "X[" << i << "]: " << x[i] << endl;
    }
    //upper
    for(int i =0; i<n-1;i++)
    {
        u[i] = (-1.0) * rho[i+1];
        l[i] = (-1.0) * rho[i+1];
    }
    cout <<"\nUsing the Tridiagonal Matrix:"<< endl;
    cout << endl;
    for(int i=0;i<n;i++)
        cout << "D[" << i << "]: " << d[i] << endl;
    cout << endl;
    for(int i=0;i<n-1;i++)
    {
        cout << "U[" << i << "]: " << u[i] << endl;
    }
    cout << endl;
    for(int i=0;i<n-1;i++)
    {
        cout << "L[" << i << "]: " << l[i] << endl;
    }
    cout << endl;
    for(int i=0;i<n;i++)
    {
        cout << "B[" << i << "]: " << b[i] << endl;
    }
    //inverse of D
    InverseD(d,n,InvD);
    
    //jacobi method to find best X
    do
    {
        jacobiX(x,jacX,d,b,u,l,n);
        ADotX(x,d,u,l,n,z);
        resid(b,d,u,l,n,z,r);
        k++;
        //rnorm
        sum = rnorm(r, n);
    }while( sum > error);
    
    cout << "\nBest Approximation of X in " << k << " Iterations:\n";
    for(int i=0;i<n;i++)
    {
        cout << "X[" << i << "]: " << x[i] << endl;
    }
    
    
    return 0;
}

double dotProduct( double x[], double y[], int n)
{
    double sum =0.0;
    
    for(int i =0; i<n; i++)
    {
        sum = sum + ( x[i] * y[i] );
    }
    return sum;
}

void addVector( double x[], double y[], int n, double *z )
{
    for(int i =0; i<n; i++)
    {
        z[i] = x[i] + y[i];
    }
}

void ADotX(const double x[], const double d[],const double u[],const double l[], int n, double *z)
{
    z[0] = ( d[0] * x[0] ) + ( u[0] * x[1] );
    for(int i =1; i<n-1; i++)
    {
        z[i] = ( l[i-1] * x[i-1] ) + ( d[i] * x[i] ) + ( u[i] * x[i+1] );
    }
    z[n-1] = ( x[n-2] * l[n-2] ) + ( d[n-1] * x[n-1] );
}

void resid(const double b[],const double d[],const double u[],const double l[],int n,const double z[], double *r)
{
    for(int i=0; i<n; i++)
    {
        r[i] = fabs( b[i] - z[i]);
    }

}

void dInvR( const double d[], const double r[], int n, double *dx)
{
    for(int i=0;i<n;i++)
    {
        dx[i] = r[i] / d[i];
    }
}

void InverseD( const double d[], int n, double *InvD)
{
    for(int i=0;i<n;i++)
    {
        InvD[i] = (1.0 / d[i]);
    }
}

void jacobiX(double *x,double *jacX, const double d[], const double b[], const double u[], const double l[], int n)
{
    jacX[0] = (b[0] - (u[0] * x[1]) ) / d[0];
    for(int i = 1; i<n-1;i++)
    {
        jacX[i] = (b[i] - (u[i] * x[i+1]) - (l[i-1] * x[i-1]) ) / d[i];
    }
    jacX[n-1] = (b[n-1] - (l[n-2] * x[n-2]) ) / d[n-1];
    for(int i=0;i<n;i++)
        {
            x[i]=jacX[i];
        }
}

double rnorm(const double r[], int n)
{
    double sum;
    sum=0.0;
    for(int i =0;i<n;i++)
    {
        sum = sum + (r[i] * r[i]);
    }
    sum = sqrt(sum);
    return sum;
}










