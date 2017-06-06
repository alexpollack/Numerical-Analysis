#define TNT_BOUNDS_CHECK 1
// #define TNT_DEBUG 1
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


// The following includes are needed too:
// "tnt_array1d.h", "tnt_i_refvec.h"
#include "tnt_i_refvec.h"
#include "tnt_array1d.h"
#include "tnt_array2d.h"

using namespace TNT;

double my_power(double a,int n) {
    
    double b;
    if (n==0) {
        b=1.0;
    } else if (n==1) {
        b=a;
    } else if (n==-1) {
        b=1.0/a;
    } else if (n>1) {
        b=a;
        for (int i=1;i<=n-1;i++)
            b*=a;
    } else if (n<-1) {
        b=1.0/a;
        for (int i=1;i<=-n-1;i++)
            b/=a;
    }
    return b;
}

void print_2d_array(Array2D<double> A) {
    
    int m=A.dim1();
    int n=A.dim2();
    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            std::cout << A[i][j] << ' ';
        }
        std::cout << '\n';
    }
}

void print_1d_array(Array1D<double> A) {
    
    int m=A.dim1();
    for (int i=0;i<m;i++) {
        std::cout << A[i] << ' ';
    }
    std::cout << '\n';
}

void matrix_vector_mult(Array2D<double> A,Array1D<double> X,
                        Array1D<double>& AX) {
    
    int m=A.dim1();
    int n=A.dim2();
    if (X.dim1()!=n) {
        std::cout << "X has invalid size";
        exit(0);
    }
    if (AX.dim1()!=m) {
        std::cout << "AX has invalid size";
        exit(0);
    }
    
    for (int i=0;i<m;i++) {
        AX[i]=0.0;
        for (int j=0;j<n;j++)
            AX[i]+=A[i][j]*X[j];
    }
}

// Z=alpha X + beta y
void vector_combine(Array1D<double>& Z,Array1D<double> X,
                    Array1D<double> Y,double alpha,double beta) {
    
    int m=Z.dim1();
    if (X.dim1()!=m) {
        std::cout << "X has invalid size";
        exit(0);
    }
    if (Y.dim1()!=m) {
        std::cout << "Y has invalid size";
        exit(0);
    }
    
    for (int i=0;i<m;i++) {
        Z[i]=alpha*X[i]+beta*Y[i];
    }
}


// warning: overwrites A
// status=1 if success
// status=0 if failure
void matrix_solve(Array2D<double> A,Array1D<double> B, Array1D<double>& X,int& status,int n)
{
    
    double alpha,holdvalue;
    int i,j,k,holdj;
    
    status=1;
    for (i=1;i<=n-1;i++) {
        holdj=i;
        holdvalue=fabs(A[i-1][i-1]);
        for (j=i+1;j<=n;j++) {
            if (fabs(A[j-1][i-1])>holdvalue) {
                holdj=j;
                holdvalue=fabs(A[j-1][i-1]);
            }
        }
        if (holdj!=i) {
            for (j=i;j<=n;j++) {
                holdvalue=A[i-1][j-1];
                A[i-1][j-1]=A[holdj-1][j-1];
                A[holdj-1][j-1]=holdvalue;
            }
        }
        holdvalue=B[i-1];
        B[i-1]=B[holdj-1];
        B[holdj-1]=holdvalue;
        if (fabs(A[i-1][i-1])<1.0E-32)
            status=0;
        else {
            for (j=i+1;j<=n;j++) {
                alpha=A[j-1][i-1]/A[i-1][i-1];
                for (k=i;k<=n;k++) {
                    A[j-1][k-1]=A[j-1][k-1]-alpha*A[i-1][k-1];
                }
                B[j-1]=B[j-1]-alpha*B[i-1];
            } // j
        }
    } // i
    
    for (i=n;i>=1;i--) {
        if (status!=0) {
            holdvalue=B[i-1];
            for (j=i+1;j<=n;j++) {
                holdvalue=holdvalue-A[i-1][j-1]*X[j-1];
            }
            if (fabs(A[i-1][i-1])<1.0E-32)
                status=0;
            else
                X[i-1]=holdvalue/A[i-1][i-1];
        }
    } // i
}


//SCALED//
void SCALED_PARTIAL(Array2D<double> A,Array1D<double> B, Array1D<double>& X,int& status,int n)
{
    
    double alpha,holdvalue;
    int i,j,k;
    
    double holdA,holdB, s[n];
    int p;
    
    //scan rows for s
    for(int i=0; i < n; i++)
    {
        s[i]=0.0;
        for(int j=0; j < n; j++)
        {
            if( fabs(A[i][j]) > s[i] )
                s[i]=fabs(A[i][j]);
        }
    }
    
    
    status=1;
    for (i=1;i<=n-1;i++)
    {
        p=0;
        //finding a p value
        for(int i2=0; i2 < n; i2++)
        {
            if( (fabs(A[p][0]))/s[p] < (fabs(A[i2][0]))/s[i2] )
                p=i2;
        }
        //swaping the pth row with the first row
        for(int j=0; j < n; j++)
        {
            holdA=A[i][j];
            A[i][j]=A[p][j];
            A[p][j]=holdA;
        }
        holdB=B[i];
        B[i]=B[p];
        B[p]=holdB;
        
        // swap the scales
        
        if (fabs(A[i-1][i-1])<1.0E-32)
            status=0;
        else
        {
            for (j=i+1;j<=n;j++)
            {
                alpha=A[j-1][i-1]/A[i-1][i-1];
                for (k=i;k<=n;k++)
                {
                    A[j-1][k-1]=A[j-1][k-1]-alpha*A[i-1][k-1];
                }
                B[j-1]=B[j-1]-alpha*B[i-1];
            } // j
        }
        
        
        
    } // i
    
    for (i=n;i>=1;i--)
    {
        if (status!=0)
        {
            holdvalue=B[i-1];
            for (j=i+1;j<=n;j++)
            {
                holdvalue=holdvalue-A[i-1][j-1]*X[j-1];
            }
            if (fabs(A[i-1][i-1])<1.0E-32)
                status=0;
            else
                X[i-1]=holdvalue/A[i-1][i-1];
        }
    } // i
}

//COMPLETE//
void complete(Array2D<double> A,Array1D<double> B, Array1D<double>& X,int& status,int n)
{
    
    double alpha,holdvalue;
    int i,j,k,a=0,b=0;
    double holdA,holdB, s;
    
    status=1;
    for (i=1;i<=n-1;i++)
    {
        s=0.0;
        //scan rows for s
        for(int i1=i-1; i1 < n; i1++)
        {
            for(int j=i-1; j < n; j++)
            {
                if( fabs(A[i1][j]) > s )
                {
                    s=fabs(A[i1][j]);
                    a=i1;
                    b=j;
                }
            }
        }
        //swaping the ath row and bth colum with the first rowand column
        //row
        for(int j=0; j < n; j++)
        {
            holdA=A[i][j];
            A[i][j]=A[a][j];
            A[a][j]=holdA;
        }
        holdB=B[i];
        B[i]=B[a];
        B[a]=holdB;
        //column
        for(int i1=0;i1<n;i1++)
        {
            holdA=A[i1][i];
            A[i1][i]=A[i][b];
            A[i1][b]=holdA;
        }
        
        if (fabs(A[i-1][i-1])<1.0E-32)
            status=0;
        else
        {
            for (j=i+1;j<=n;j++)
            {
                alpha=A[j-1][i-1]/A[i-1][i-1];
                for (k=i;k<=n;k++)
                {
                    A[j-1][k-1]=A[j-1][k-1]-alpha*A[i-1][k-1];
                }
                B[j-1]=B[j-1]-alpha*B[i-1];
            } // j
        }
    } // i
    
    for (i=n;i>=1;i--)
    {
        if (status!=0)
        {
            holdvalue=B[i-1];
            for (j=i+1;j<=n;j++)
            {
                holdvalue=holdvalue-A[i-1][j-1]*X[j-1];
            }
            if (fabs(A[i-1][i-1])<1.0E-32)
                status=0;
            else
                X[i-1]=holdvalue/A[i-1][i-1];
        }
    } // i
}

int main() {
    
    int n=5;
    Array2D<double> A(n,n,0.0);
    for (int i=0;i<n;i++) {
        double alpha=(i+0.5)/n;
        for (int j=0;j<n;j++) {
            // Van Der Monde Matrix
            //   A[i][j]=my_power(alpha,j);
            // Hilbert Matrix
            A[i][j]=1.0/(i+j+1.0);
        }
    }
    Array1D<double> x(n,0.0);
    Array1D<double> B(n,0.0);
    B[n-1]=1.0;
    
    Array1D<double> X(n,0.0);
    int status;
    std::cout << "A is:\n";
    print_2d_array(A);
    std::cout << "B is:\n";
    print_1d_array(B);
    Array2D<double> ASAVE(n,n,0.0);
    Array1D<double> BSAVE(n,0.0);
    BSAVE=B.copy();
    ASAVE=A.copy();
    //matrix_solve(A,B,X,status,n);
    //std::cout << "MATRIX SOLVE\n";
    
    SCALED_PARTIAL(A,B,X,status,n);
    std::cout << "SCALED PARTIAL PIVOT\n";
    
    //complete(A,B,X,status,n);
    //std::cout << "COMPLETE PIVOT\n";
    if (status==1)
        std::cout << "success\n";
    else
        std::cout << "failure\n";
    
    std::cout << "X is: \n";
    print_1d_array(X);
    Array1D<double> AX(n,0.0);
    Array1D<double> RESID(n,0.0);
    matrix_vector_mult(ASAVE,X,AX);
    vector_combine(RESID,BSAVE,AX,1.0,-1.0);
    std::cout << "RESID is:\n";
    print_1d_array(RESID);
    
    return 0;
}






































