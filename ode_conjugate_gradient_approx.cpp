#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <iomanip>

using namespace std;

double integratePx(double x);
double px(double x);
double qx(double x);
double fx(double x);
double dotProduct( double x[], double y[], int n);
void ADotX(const double x[], const double d[],const double u[],const double l[], int n, double *z);
void resid(const double b[],const double d[],const double u[],const double l[],int n,const double z[], double *r);
double rnorm( const double r[], int n);


int main()
{
    int N=5;
    int n=N,i;
    double d[N+1],u[N],l[N],sum=0.0,y0=1.0;
    double xi[N+1],x[N+1],iphi[N+1],jphi[N+1],A[N+1][N+1],b[N+1],H=(1.0/N);
    double z[n+1], r[n+1],epsilon=1.0, alpha=2.0, beta=3.0,rho[n+1],h,Ri,c[N+1];
    double error = .000001;
    double alphacg,betacg;
    double cgr[n+1],cgp[n+1],cgx[n+1];
    
    ofstream fileOut;
    fileOut.open("ODEplotsolution.txt");
    for(int i=0;i<N+1;i++)
        x[i]=0.0;
    for(int i=0;i<N+1;i++)
        xi[i]=i*H;
    
    //ph
    int n1=1,n2=1;
    for(int i=0;i<N+1;i++)
    {
        n2=1;
        for(int j=0;j<N+1;j++)
        {
            sum=0.0;
            for(int k=1;k<N+1;k++)
            {
                if( xi[i-1] < (xi[k]+xi[k-1])/2.0 && (xi[k]+xi[k-1])/2.0 < xi[i])
                    iphi[k]=(1.0/H);
                else if( xi[i] < (xi[k]+xi[k-1])/2.0 && (xi[k]+xi[k-1])/2.0 < xi[i+1])
                    iphi[k]=(-1.0/H);
                else
                    iphi[k]=0.0;
                if(  xi[j-1] < (xi[k]+xi[k-1])/2.0 && (xi[k]+xi[k-1])/2.0 < xi[j])
                    jphi[k]=(1.0/H);
                else if(xi[j] < (xi[k]+xi[k-1])/2.0 && (xi[k]+xi[k-1])/2.0 < xi[j+1])
                    jphi[k]=(-1.0/H);
                else
                    jphi[k]=0.0;
                sum=sum+(px((xi[k]+xi[k-1]/2.0))*iphi[k]*jphi[k]);
            }
            if(i==j)
                d[i]=sum;
            else if(j==i+1)
                u[i]=sum;
            else if(j==i-1)
                l[j]=sum;
            A[i][j]=sum;
            n2++;
        }
        n1++;
    }
    //qh
    n1=1,n2=1;
    for(int i=0;i<N;i++)
    {
        n2=1;
        for(int j=0;j<N;j++)
        {
            sum=0.0;
            for(int k=1;k<N;k++)
            {
                if( xi[i-1] < (xi[k]+xi[k-1])/2.0 && (xi[k]+xi[k-1])/2.0 < xi[i])
                    iphi[k]=(1.0/H)*( (xi[k]+xi[k-1])/2.0 - xi[i]);
                else if( xi[i] < (xi[k]+xi[k-1])/2.0 && (xi[k]+xi[k-1])/2.0 < xi[i+1])
                    iphi[k]=(1.0/H)*( xi[i]- (xi[k]+xi[k-1])/2.0);
                else
                    iphi[k]=0.0;
                if(  xi[j-1] < (xi[k]+xi[k-1])/2.0 && (xi[k]+xi[k-1])/2.0 < xi[j])
                    jphi[k]=(1.0/H)*( (xi[k]+xi[k-1])/2.0 - xi[j]);
                else if(xi[j] < (xi[k]+xi[k-1])/2.0 && (xi[k]+xi[k-1])/2.0 < xi[j+1])
                    jphi[k]=(1.0/H)*( xi[j]- (xi[k]+xi[k-1])/2.0);
                else
                    jphi[k]=0.0;
                
                sum=sum+(qx((xi[k]+xi[k-1]/2.0))*iphi[k]*jphi[k]);
            }
            if(i==j)
                d[i]=d[i]+sum;
            else if(j==i+1)
                u[i]=u[i]+sum;
            else if(j==i-1)
                l[j]=l[j]+sum;
            A[i][j]=A[i][j]+sum;
            n2++;
        }
        n1++;
    }
    sum=0.0;
    for(int k=1;k<N+1;k++)
    {
        if( xi[0] < (xi[k]+xi[k-1])/2.0 && (xi[k]+xi[k-1])/2.0 < xi[1])
            iphi[k]=(1.0/H)*( xi[1] - (xi[k]+xi[k-1])/2.0);
        else
            iphi[k]=0.0;
        sum=sum+(fx((xi[k]+xi[k-1]/2.0))*iphi[k]) - y0*d[0];
    }
    b[0]=sum;
    sum=0.0;
    for(int i=1;i<N+1;i++)
    {
        sum=0.0;
        for(int k=1;k<N+1;k++)
        {
            if( xi[i-1] < (xi[k]+xi[k-1])/2.0 && (xi[k]+xi[k-1])/2.0 < xi[i])
                iphi[k]=(1.0/H)*( (xi[k]+xi[k-1])/2.0 - xi[1]);
            else if( xi[i] < (xi[k]+xi[k-1])/2.0 && (xi[k]+xi[k-1])/2.0 < xi[i+1])
                iphi[k]=(1.0/H)*( xi[1] - (xi[k]+xi[k-1])/2.0);
            else
                iphi[k]=0.0;
            sum=sum+(fx((xi[k]+xi[k-1]/2.0))*iphi[k]);
        }
        b[i]=sum;
    }
    
    cout << "Conjugate Gradient Method\n";
    h= (1.0/n);
    cout << "\nInitial Guess 'X0'\n";
    for(int i=0;i<n+1;i++)
    {
        cout << "X0[" << i << "]: " << x[i] << endl;
    }
    cout <<"\nUsing the Tridiagonal Matrix:"<< endl;
    cout << "A:\n";
    for(int i=0;i<N+1;i++)
    {
        cout <<"|";
        for(int j=0;j<N+1;j++)
        {
            cout << A[i][j];
            if(j<N)
                cout << "\t";
        }
        cout << "|\n";
    }
    cout << endl;
    for(int i=0;i<n+1;i++)
    {
        cout << "B[" << i << "]: " << b[i] << endl;
    }
    ADotX(x,d,u,l,n,z);
    resid(b,d,u,l,n,z,r);
    for(int i=0;i<n+1;i++)
        cgp[i]=r[i];
    int end=0,it=0;
    do
    {
        ADotX(cgp,d,u,l,n,z);
        alphacg = dotProduct(r,r,n)/dotProduct(cgp,z,n);
        for(int i=0;i<n+1;i++)
            cgx[i]=x[i]+alphacg*cgp[i];
        for(int i=0;i<n+1;i++)
            cgr[i]=r[i]-alphacg*z[i];
        if(fabs(rnorm(r,n))<error)//.001
            end=1;
        else
        {
            betacg=dotProduct(cgr,cgr,n)/dotProduct(r,r,n);
            for(i=0;i<n+1;i++)
                cgp[i]=cgr[i]+betacg*cgp[i];
            for(int i=0;i<n+1;i++)
                r[i]=cgr[i];
            for(int i=0;i<n+1;i++)
                x[i]=cgx[i];
        }
        
        it++;
    }while(end==0);
    for(int i=0;i<N+1;i++)
        c[i]=x[i];
    
    cout << "\nBest Approximation of X in " << it << " Iterations:\n";
    for(int i=0;i<n+1;i++)
    {
        cout << "X[" << i << "]: " << x[i] << endl;
    }
    cout << "\nPlotted Pair: (xi,ci)" << it << " Iterations:\n";
    for(int i=0;i<n+1;i++)
    {
        cout << "(" << xi[i] << ", "<< c[i]<< ")\n";
    }

    for(int i=0;i<N+1;i++)
        fileOut << xi[i] << "\t" << c[i] << endl;
    
    fileOut.close();
    
        return 0;
}

double px(double x)
{
    double p=2.0*x*x;
    return p;
}

double integratePx(double x)
{
    return x;
}

double qx(double x)
{
    double q=x;
    return q;
}

double fx(double x)
{
    double f=0.0;
    return f;
}

double dotProduct( double x[], double y[], int n)
{
    double sum =0.0;
    
    for(int i =0; i<n+1; i++)
    {
        sum = sum + ( x[i] * y[i] );
    }
    return sum;
}

void ADotX(const double x[], const double d[],const double u[],const double l[], int n, double *z)
{
    z[0] = ( d[0] * x[0] ) + ( u[0] * x[1] );
    for(int i =1; i<n+1; i++)
    {
        z[i] = ( l[i-1] * x[i-1] ) + ( d[i] * x[i] ) + ( u[i] * x[i+1] );
    }
    z[n+2] = ( x[n-1] * l[n-1] ) + ( d[n] * x[n] );
}

void resid(const double b[],const double d[],const double u[],const double l[],int n,const double z[], double *r)
{
    for(int i=0; i<n+1; i++)
    {
        r[i] = fabs( b[i] - z[i]);
    }
    
}

double rnorm(const double r[], int n)
{
    double sum;
    sum=0.0;
    for(int i =0;i<n+1;i++)
    {
        sum = sum + (r[i] * r[i]);
    }
    sum = sqrt(sum);
    return sum;
}












