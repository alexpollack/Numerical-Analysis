#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

double fx_smooth(double x) {
    return 1.0/(1.0+x*x);
}


double fx_tanh(double x)
{
    double epsilon = 1.0, f_trans;
    double x_eps = (x/epsilon);
    f_trans = tanh(x_eps);
    
    return f_trans;
}

// return y=P_n(x)   P_n interpolates the points in xsten.
void lagrange_interp(double* xsten,double (*fx)(double),
                     int n,double x,double& y) {
    
    double* fsten=new double[n+1];
    int i,j;
    double L;
    
    for (i=0;i<=n;i++)
    {
        fsten[i]=fx(xsten[i]);
    }
    y=0.0;
    for (i=0;i<=n;i++)
    {
        L=1.0;
        for (j=0;j<=n;j++)
        {
            if (i!=j)
            {
                L=L*(x-xsten[j])/(xsten[i]-xsten[j]);
            }
        }
        y=y+fsten[i]*L;
    }
    
    delete[] fsten;
    
}

//NEWTONS DIVIDED DIFFERENCE METHOD
void divided_diff(double* xsten,double (*fx)(double),int n,double x,double& y)
{
    
    double* fsten=new double[n+1];
    int i,N=n+1;
    double Z=1.0;
    y=0.0;
    double F[N][N];
    
    for (i=0;i<=n;i++)
    {
        fsten[i]=fx(xsten[i]);
    }
    for (int i = 0; i<=n;i++)
    {
        F[i][0]=fsten[i];
    }
    
    for (int i = 1; i<=n; i++)
    {
        for (int j=1; j<=i;j++)
        {
            F[i][j]=(F[i][j-1]-F[i-1][j-1])/(xsten[i]-xsten[i-j]);
        }
    }
    for (int i=0;i<=n;i++)
    {
        if (i>0)
            Z=Z*(x-xsten[i-1]);
        y=y+F[i][i]*Z;
    }
    delete[] fsten;
    
}

//ENO
void ENO(double* xsten,double (*fx)(double),int n,double x,double& y)
{
    
    double* fsten=new double[n+1];
    int i,N=n+1,k[N],kL, kH, m=1;
    double Z=1.0,z[N],G[N],a,b;
    
    double F[N][N];
    
    for (i=0;i<=n;i++)
    {
        fsten[i]=fx(xsten[i]);
    }
    for (int i = 0; i<=n;i++)
    {
        F[i][0]=fsten[i];
    }
    
    for (int i = 1; i<=n; i++)
    {
        for (int j=1; j<=i;j++)
        {
            F[i][j]=(F[i][j-1]-F[i-1][j-1])/(xsten[i]-xsten[i-j]);
        }
    }
    
    
    y=0.0;
    for (i=0;i<=n;i++)
    {
        if (xsten[i]<=x && xsten[i+1] >=x)
        {
            k[0]=i;
        }
    }
    k[1] = k[0] + 1;
    kL = k[0]-1;
    kH = k[1] + 1;
    z[0] = xsten[k[0]];
    z[1] = xsten[k[1]];
    
    G[0] = F[k[0]][0];
    G[1] = F[k[1]][1];
    
    do
    {
        m = m + 1;
        a = F[kH-1][m];
        b = F[kH][m];
        if (fabs(a) < fabs(b) )
        {
            k[m] = kL;
            kL = (kL-1);
            G[m] = a;
        }
        else
        {
            k[m] = kH;
            kH = (kH + 1);
            G[m] = b;
        }
        
        z[m] = xsten[k[m]];
        
        
        
    }
    while( kL >= 0 && kH <= n);
    double S=1.0,C,R=1.0;
    Z=0.0;
    for (int i=0;i<=m;i++)
    {
        R=1.0;
        if(i!=0)
        {
            for(int j=0;j<=i-1;j++)
            {
                Z = -z[j];
                S=x+Z;
                R=R*S;
            }
        }
        C=G[i]*R;
        y=y+C;
    }
    
    delete[] fsten;
    
}


int main()
{
    
    int n=7;
    double fTrans;
    
    double xlo=-5.0;
    double xhi=5.0;
    double h=(xhi-xlo)/n;
    double* xsten=new double[n+1];
    for (int i=0;i<=n;i++) {
        xsten[i]=xlo+i*h;
    }
    
    ofstream lagrangefile;
    lagrangefile.open("lagrange.dat");
    ofstream originalfile;
    originalfile.open("original.dat");
    ofstream lagrangeTransfile;
    lagrangeTransfile.open("lagrange_Trans.dat");
    ofstream fTransFile;
    fTransFile.open("f_Transfer.dat");
    ofstream NewtonsTableFile;
    NewtonsTableFile.open("NewtonsTable.dat");
    ofstream NewtonsTanhFile;
    NewtonsTanhFile.open("NewtonsTanh.dat");
    ofstream EnoFile;
    EnoFile.open("ENO.dat");
    ofstream EnoTransFile;
    EnoTransFile.open("EnoTrans.dat");
    cout << "output filenames are lagrange.dat, original.dat, f_Transfer.dat \n";
    
    int nplot=1000;
    double hplot=(xhi-xlo)/nplot;
    
    for (int i=0;i<=nplot;i++) {
        double xplot=xlo+i*hplot;
        double fplot, foriginal;
        lagrange_interp(xsten,fx_smooth,n,xplot,fplot);
        lagrangefile << xplot << '\t' << fplot << '\n';
        divided_diff(xsten,fx_smooth,n,xplot,fplot);
        NewtonsTableFile << xplot << '\t' << fplot << '\n';
        ENO(xsten,fx_smooth,n,xplot,fplot);
        EnoFile << xplot << '\t' << fplot << '\n';
        foriginal=fx_smooth(xplot);
        originalfile << xplot << '\t' << foriginal << '\n';
        lagrange_interp(xsten,fx_tanh,n,xplot,fplot);
        lagrangeTransfile << xplot << '\t' << fplot << '\n';
        divided_diff(xsten,fx_tanh,n,xplot,fplot);
        NewtonsTanhFile << xplot << '\t' << fplot << '\n';
        ENO(xsten,fx_tanh,n,xplot,fplot);
        EnoTransFile << xplot << '\t' << fplot << '\n';
        fTrans = fx_tanh(xplot);
        fTransFile << xplot << '\t' << fTrans << '\n';
        
    }
    
    lagrangefile.close();
    originalfile.close();
    lagrangeTransfile.close();
    fTransFile.close();
    NewtonsTableFile.close();
    NewtonsTanhFile.close();
    EnoFile.close();
    EnoTransFile.close();
    
    delete[] xsten;
    
    
}






