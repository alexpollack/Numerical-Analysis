#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

void nfactorial(int n,double& nfact) {
    
    int i;
    
    nfact=1.0;
    for (i=1;i<=n;i++) {
        nfact=nfact*i;
    }
    
}

double power(double x,int N) {
    
    double pow=1.0;
    for (int i=1;i<=N;i++)
        pow*=x;
    
    return pow;
    
}

double PN(int N,double x) {
    double nfact,term;
    
    double y=1.0;
    for (int i=1;i<=N;i++) {
        nfactorial(i,nfact);
        term=power(x,i)/nfact;
        y=y+term;
    }
    return y;
    
}

int main() {
    
    int N,Nplot;
    double y,exact,x,a,b,h;
    
    N=8;
    x=0.25;
    y=PN(N,x);
    exact=exp(x);
    std::cout << "x= " << x << " and N= " << N << '\n';
    std::cout << "Taylor series approximation is: " << y << '\n';
    std::cout << "exp(x) is: " << exact << '\n';
    std::cout << "error is: " << fabs(exact-y) << '\n';
    ofstream PNfile;
    PNfile.open("PNfile");
    ofstream expfile;
    expfile.open("expfile");
    std::cout << "The file `PNfile' has PN(x) 0<x<1\n";
    std::cout << "The file `expfile' has exp(x) 0<x<1\n";
    a=0.0;
    b=1.0;
    Nplot=100;
    h=(b-a)/Nplot;
    for (int iplot=0;iplot<=Nplot;iplot++) {
        x=a+h*iplot;
        y=PN(N,x);
        PNfile << x << ' ' << y << '\n';
        y=exp(x);
        expfile << x << ' ' << y << '\n';
    }
    PNfile.close();
    expfile.close();
}
