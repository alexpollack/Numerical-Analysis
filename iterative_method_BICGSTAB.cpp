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
double rnorm( const double r[], int n);
void Gx(const double x[], const double z[], const double b[], double *gX, int n);

int main()
{
    int n=5, k=0,i, brk=0;
    double x[n],b[n],d[n],u[n-1],l[n-1], z[n], r[n];
    double epsilon=1.0, alpha=2.0, beta=0.2,rho[n+1],h,Ri,error = .0000001;
    double bir[n],rhoVect[n],alphaBI,omega[n],v[n],p[n], betaBI,H[n],t[n],s[n],check[n],test[n];
    
    cout << "Biconjugate Gradient Stabilized Method\n";
    
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
    cout << endl;
    
    ADotX(x,d,u,l,n,z);
    resid(b,d,u,l,n,z,r);
    for(int i=0;i<n;i++)
        bir[i]=r[i];
    for(int i=0;i<n;i++)
        v[i]=p[i]=0.0;
    rhoVect[0]=omega[0]=alphaBI=1.0;
    i=1;
    do
    {
        rhoVect[i]=dotProduct(r,bir,n);
        betaBI= (rhoVect[i] / rhoVect[i-1]) * (alphaBI/omega[i-1]);
        for(int j=0;j<n;j++)
        {
            p[j]=bir[j] + betaBI*(p[j]-omega[i-1]*v[j]);
        }
        ADotX(p,d,u,l,n,v);
        alphaBI= rhoVect[i] / dotProduct(r,v,n);
        for(int j=0;j<n;j++)
        {
            H[j]=x[j]+alphaBI*p[j];
        }
        ADotX(H,d,u,l,n,check);
        resid(b,d,u,l,n,check,test);
        if(fabs(rnorm(test,n)) < error)
        {
            brk=1;
            for(int j=0;j<n;j++)
                x[j]=H[j];
        }
        else
        {
            for(int j=0;j<n;j++)
                s[j]=bir[j]-alphaBI*v[j];
            ADotX(s,d,u,l,n,t);
            omega[i]=dotProduct(t,s,n)/dotProduct(t,t,n);
            for(int j=0;j<n;j++)
                x[j]=H[j]+omega[i]*s[j];
        }
        ADotX(x,d,u,l,n,check);
        resid(b,d,u,l,n,check,test);
        if(brk !=0 || fabs(rnorm(test,n)) < error)
            brk=1;
        else
        {
            for(int j=0;j<n;j++)
                bir[j]=s[j]-omega[i]*t[j];
        }
        
                //rnorm
        k++;
        i++;
        //sum = rnorm(r, n);
    }while( brk==0);
    
    cout << "Best Approximation of X in " << k << " Iterations:\n";
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

void Gx(const double x[], const double z[], const double b[], double *gX, int n)
{
    for(int i=0; i<n;i++)
    {
        gX[i]= (x[i] * z[i]) - (2.0 * b[i] * x[i]);
    }
}




