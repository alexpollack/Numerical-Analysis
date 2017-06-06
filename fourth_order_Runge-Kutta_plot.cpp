#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

double Fyt (double Y)
{
    double Fyt;
    Fyt=Y+1;
    return Fyt;
}

double exact (double t)
{
    double exact;
    exact=2*(pow(2.718,t)) - 1;
    return exact;
}


int main()
{
    ofstream exactFile;
    ofstream apprFile;
    
    exactFile.open("exactFile.txt");
    apprFile.open("apprFile.txt");
    
    int k=0,N=5;
    double t=0.0, T=1.0;
    double* Y=new double[N+1];
    
    double delta_t=T/N;
    double tK;
    
    Y[0]=1;
    double Z,k1,k2,k3,k4;
    int j=1;
    apprFile << t <<"\t"<< Y[0]<<endl;
    
    for(int k=0;k<N;k++)
    {
        tK=k*delta_t;
        
        k1=delta_t*Fyt(Y[k]);
        
        k2=delta_t*Fyt((Y[k]+k1*.5));
        
        k3=delta_t*Fyt((Y[k]+k2*.5));
        
        k4=delta_t*Fyt((Y[k]+k3));
        
        Y[j]=Y[k] + (k1 + 2*k2 + 2*k3 + k4)/6.0;
        
        t=t+delta_t;
        apprFile << t <<"\t"<< Y[j]<<endl;
        
        j++;
    }
 
    t=0.0, T=1.0;
    delta_t=T/N;
    int b=1;
    for (int i=0; i<=N;i++)
    {
        Z=exact(t);
        exactFile << t <<"\t"<< Z <<endl;
        
        t=t+delta_t;
    }

    
    exactFile.close();
    apprFile.close();
    
    
    
    return 0;
}