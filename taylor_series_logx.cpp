#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <iomanip>
using namespace std;

int nFact(int n);

void PnXToFile(int N, double x0, double a, double b);

void PnX(int N, double x0, double x, long double exact);

int main()
{
    
    double x0 =.5, a = .25, b = .75, funcXL = 2 * a;
    long double exact = log(funcXL);
    int N = 8;

    
    cout << "Exact f(x) is: " << exact << endl;
    cout << "x0= " << x0 << ", x= " << a << ", N= " << N << '\n';
    cout << "Taylor series approximation is:\nN\tapprox \terror\n";
    PnX(N, x0, a, exact);
    cout << '\n';

    
    PnXToFile(N, x0, a, b);
    
    return 0;
    
}

void PnX(int N, double x0, double x, long double exact)
{
    long double pnxVal = 0, fkx, J = x - x0, error;
    
    long double top1, top2, top3, bot1, bot2;
    
    for (int n = 1; n < N; n++)
    {
        int j = n + 1, k = n - 1;
        
        top1 =pow(-1, j);
        top2=pow(J, n); top3=nFact(k); bot1=pow(x0, n); bot2=nFact(n);
        top3=1.0;
        bot2=n;
        
        fkx =  (top1*top2*top3)/(bot1*bot2);
        //fkx =  (pow(-1, j)*pow(J, n)*nFact(k))/(pow(x0, n)*nFact(n));
        

        pnxVal = pnxVal + fkx;
        error = fabs(exact - pnxVal);
        cout << n << "\t";
        cout << pnxVal << "\t" << error << endl;
        j++, k++;
        error = 0;
    }
    
    
}

int nFact(int n)
{
    int nfact=1;
    for (int i=1 ;i <=  n; i++)
    {
        nfact = nfact*i;
    }
    return nfact;
}

void PnXToFile(int N, double x0, double a, double b)
{
    long double pnxVal = 0, fkx, exact;
    
    long double top1, top2, top3, bot1, bot2, x;
    
    ofstream PNfile;
    PNfile.open("PNfile.txt");
    ofstream exactFile;
    exactFile.open("exactFile.txt");
    
    
    for (int iplot=0;iplot<51;iplot++)
    {
        x=a+(iplot*.01);
        double J = x - x0;
    
        for (int n = 1; n < N; n++)
            {
                int j = n + 1;
        
                top1 =pow(-1, j);
                top2=pow(J, n);
                bot1=pow(x0, n);
                top3=1.0;
                bot2=n;
                fkx =  (top1*top2*top3)/(bot1*bot2);
        
                pnxVal = pnxVal + fkx;
            }
        double X = 2*x;
        exact = log(X);
    
        PNfile << fixed << showpoint << setprecision(2) << x << "\t";
        PNfile << fixed << showpoint << setprecision(9) << pnxVal << endl;
        exactFile << fixed << showpoint << setprecision(2)<< x << "\t";
        exactFile << fixed << showpoint << setprecision(9)<< exact << endl;
        pnxVal = 0;
    }
    PNfile.close();
    exactFile.close();
}

