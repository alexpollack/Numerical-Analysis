#include <iostream>
#include <cmath>
#include <math.h>

using namespace std;

//L2 norm function declaration
double norm(double u[], int n);

int main()
{
    double a1[4]={3.0,3.0,3.0,3.0},a2[4]={2.0,-1.0,2.0,-1.0},a3[4]={12.0,8.0,4.0,0.0};
    double h1[4],h2[4],h3[4],h4[4],H1[4][4],H1A[4][3],H2[4][4],v1[4],v2[3];
    double e1[4]={1.0,0.0,0.0,0.0},e2[3]={1.0,0.0,0.0};
    double OP[4][4],temp[3][3],I4[4][4]={{1.0,0.0,0.0,0.0},{0.0,1.0,0.0,0.0},{0.0,0.0,1.0,0.0},{0.0,0.0,0.0,1.0}};
    double R2[3][3],Q[4][4],R[4][3],sum=0.0,altA2[3],I3[3][3]={{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
    
    //Print A to screen
    cout << "A:\n";
    for(int i=0;i<4;i++)
    {
        cout << "|" << a1[i] << "\t" << a2[i] << "\t" << a3[i] << "|\n";
    }
    //V1
    //First step in the first householder
    for(int i=0;i<4;i++)
        v1[i]= a1[i] + (norm(a1,4)*e1[i]);
    //Outer product
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
            OP[i][j]=(v1[i]*v1[j]);
    }
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
            OP[i][j]= (2.0/pow(norm(v1,4),2)) * OP[i][j];
    }
    //The first householder matrix H1
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
            H1[i][j]=I4[i][j] - OP[i][j];
    }
    //H1dotA
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
            sum=sum+(H1[i][j]*a1[j]);
        if(fabs(sum)<.001)
            H1A[i][0]=0.0;
        else
            H1A[i][0]=sum;
        sum =0.0;
    }
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
            sum=sum+(H1[i][j]*a2[j]);
        if(fabs(sum)<.001)
            H1A[i][1]=0.0;
        else
            H1A[i][1]=sum;
        sum =0.0;
    }
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
            sum=sum+(H1[i][j]*a3[j]);
        if(fabs(sum)<.001)
            H1A[i][2]=0.0;
        else
            H1A[i][2]=sum;
        sum =0.0;
    }
    //First step in finding the second householder matrix
    for(int i=0;i<3;i++)
        altA2[i]=H1A[i+1][1];
    //V2
    for(int i=0;i<4;i++)
        v2[i]= altA2[i] + (norm(altA2,3)*e2[i]);
    //Outer product
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
            temp[i][j]=(v2[i]*v2[j]);
    }
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
            temp[i][j]= (2.0/pow(norm(v2,3),2)) * temp[i][j];
    }
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
            R2[i][j]=I3[i][j] - temp[i][j];
    }
    //H2, setting non diagnal row1 and col1 to zero (c++ saves as 0.00000xxx)
    H2[0][0]=1.0;
    for(int i=1;i<4;i++)
        H2[0][i]=0.0;
    for(int i=1;i<4;i++)
        H2[i][0]=0.0;
    for(int i=1;i<4;i++)
    {
        for(int j=1;j<4;j++)
            H2[i][j]=R2[i-1][j-1];
    }
    //Q=H2H1, Q is the product matrix of all H's
    for(int i=0;i<4;i++)
    {
        h1[i]=H1[i][0];
        h2[i]=H1[i][1];
        h3[i]=H1[i][2];
        h4[i]=H1[i][3];
    }
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
            sum=sum+(H2[i][j]*h1[j]);
        Q[i][0]=sum;
        sum =0.0;
    }
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
            sum=sum+(H2[i][j]*h2[j]);
        Q[i][1]=sum;
        sum =0.0;
    }
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
            sum=sum+(H2[i][j]*h3[j]);
        Q[i][2]=sum;
        sum =0.0;
    }
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
            sum=sum+(H2[i][j]*h4[j]);
        Q[i][3]=sum;
        sum =0.0;
    }
    //Print Q to screen
    cout << "Q:\n";
    for(int i=0;i<4;i++)
    {
        cout << "|";
        for(int j=0;j<4;j++)
        {
            cout << Q[i][j];
            if(j<=2)
                cout<<" ";
        }
        cout <<"|" << endl;
    }
    //R=H2H1A, R is the product matrix of all H and A
    //(Or R is the product of Q and A)
    sum=0.0;
    for(int i=0;i<4;i++)
    {
        a1[i]=H1A[i][0];
        a2[i]=H1A[i][1];
        a3[i]=H1A[i][2];
    }
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
            sum=sum+(H2[i][j]*a1[j]);
        if(fabs(sum)<.001)
            R[i][0]=0.0;
        else
            R[i][0]=sum;
        sum =0.0;
    }
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
            sum=sum+(H2[i][j]*a2[j]);
        if(fabs(sum)<.001)
            R[i][1]=0.0;
        else
            R[i][1]=sum;
        sum =0.0;
    }
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
            sum=sum+(H2[i][j]*a3[j]);
        if(fabs(sum)<.001)
            R[i][2]=0.0;
        else
            R[i][2]=sum;
        sum =0.0;
    }
    //Print R to screen
    cout << "R:\n";
    for(int i=0;i<4;i++)
    {
        cout << "|";
        for(int j=0;j<3;j++)
        {
            cout << R[i][j];
            if(j<=1)
                cout << "\t";
        }
        cout << "|" << endl;
    }
    
    return 0;
}

//Function that will compute&return the L2 norm of a vector
double norm(double u[], int n)
{
    double uNorm,sum=0.0;
    
    for(int i=0;i<n;i++)
        sum= sum + (u[i]*u[i]);
    uNorm =sqrt(sum);
    return uNorm;
}















