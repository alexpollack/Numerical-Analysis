#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <iomanip>

using namespace std;

double dotProduct( double x[], double y[], int n);
void ADotX(const double x[], const double d[],const double u[],const double l[], int n, double *z);
void resid(const double b[],const double d[],const double u[],const double l[],int n,const double z[], double *r);
double rnorm( const double r[], int n);

int main()
{
    int n=5,i;
    double x[n],b[n],d[n],u[n-1],l[n-1], z[n], r[n];
    double epsilon=1.0, alpha=2.0, beta=3.0,rho[n],h,Ri;
    double error = .000001;
    double alphacg,betacg;
    double cgr[n],cgp[n],cgx[n];
    cout << "Conjugate Gradient Method\n";
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
        x[i]=0.0;
    }
    cout << "\nInitial Guess 'X0'\n";
    for(int i=0;i<n;i++)
    {
        cout << "X0[" << i << "]: " << x[i] << endl;
    }
    //upper
    for(int i =0; i<n-1;i++)
    {
        u[i] = (-1.0) * rho[i+1];
        l[i] = (-1.0) * rho[i+1];
    }
    cout <<"\nUsing the Tridiagonal Matrix:"<< endl;
    for(int i=0;i<n;i++)
        cout << "D[" << i << "]: " << d[i] << endl;
    for(int i=0;i<n-1;i++)
    {
        cout << "U[" << i << "]: " << u[i] << endl;
    }
    for(int i=0;i<n-1;i++)
    {
        cout << "L[" << i << "]: " << l[i] << endl;
    }
    cout << endl;
    for(int i=0;i<n;i++)
    {
        cout << "B[" << i << "]: " << b[i] << endl;
    }
    ADotX(x,d,u,l,n,z);
    resid(b,d,u,l,n,z,r);
    for(int i=0;i<n;i++)
        cgp[i]=r[i];
    int end=0,it=0;
    do
    {
        ADotX(cgp,d,u,l,n,z);
        alphacg = dotProduct(r,r,n)/dotProduct(cgp,z,n);
        for(int i=0;i<n;i++)
            cgx[i]=x[i]+alphacg*cgp[i];
        for(int i=0;i<n;i++)
            cgr[i]=r[i]-alphacg*z[i];
        if(fabs(rnorm(r,n))<error)//.001
            end=1;
        else
        {
            betacg=dotProduct(cgr,cgr,n)/dotProduct(r,r,n);
            for(i=0;i<n;i++)
                cgp[i]=cgr[i]+betacg*cgp[i];
            for(int i=0;i<n;i++)
                r[i]=cgr[i];
            for(int i=0;i<n;i++)
                x[i]=cgx[i];
        }
        
        it++;
    }while(end==0);
    
    cout << "\nBest Approximation of X in " << it << " Iterations:\n";
    for(int i=0;i<n;i++)
    {
        cout << "X[" << i << "]: " << x[i] << endl;
    }
    cout << endl;
    

    
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


