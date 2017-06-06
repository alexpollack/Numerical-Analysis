#include <iostream>
#include <cmath>
#include <math.h>
#include <iomanip>

using namespace std;

double norm(double u[], int n);
double normInf(double u[], int n);

int main()
{
    int p,q,n=2,N=4,flag,it=0;
    double d,A[N][N],A1[N][N],a[N],v[N],V[N][N],e[N],I[N][N],H[N][N],R[N][N],Q[N][N],Qtot[N][N],sum=0.0,lower[N-1];
    double sigma[N][N],usvd[N][N],vsvd[N][N],x[N],eV[N],Asave[N][N],alpha,temp[N],beta,hold,vholder[N][N],tmp[N][N];
    double eig[N][N],nrm[N],MPSI[N][N],sigmainv[N][N],zerodiff=0;
    double b[N];
    double Arecons[N][N],AA[N][N],reorder,repeat[N];
    //b
    for(int i=0;i<N;i++)
        b[i]=1.0;
    //Identity matrix
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            if(i==j)
                I[i][j]=1.0;
            else
                I[i][j]=0.0;
        }
    }
    //Q=I
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
            Q[i][j]=I[i][j];
    }
    //usvd=I
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
            usvd[i][j]=I[i][j];
    }
    //usvd=I
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
            vsvd[i][j]=I[i][j];
    }
    //Matrix A
    p=0;
    for(int j=0;j<n;j++)
    {
        for(int i=0;i<n;i++)
        {
            for(q=0;q<N;q++)
                A[p][q]=0.0;
            d=0.0;
            if(i>0)/////
            {
                A[p][p-1]=-1.0;///
                d=d+1.0;
            }
            if(i<n-1)//
            {
                A[p][p+1]=-1.0;
                d=d+1.0;
            }
            if(j>0)////
            {
                A[p][p-n]=-1.0;
                d=d+1.0;
            }
            if(j<n-1)
            {
                A[p][p+n]=-1.0;
                d=d+1.0;
            }
            A[p][p]=d;
            p=p+1;
        }
    }
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
            Asave[i][j]=A[i][j];
    }
    cout << "A:\n";
    for(int i =0;i<N;i++)
    {
        cout << "|";
        for(int j=0;j<N;j++)
        {
            cout << A[i][j];
            if(j<N-1)
                cout <<"\t";
        }
        cout << "|\n";
    }
    cout << "\nb:\n";
    for(int i=0;i<N;i++)
        cout << "|"<<b[i]<<"|\n";
    int c=0;
     // repeat A^k=Q_k R_k (Householder)  A^k+1 = R_k Q_k
     // A^k+1=RQ=Q_k^T A^k Q_k=Q_k^T Q_k-1^T A^k-1 Q_k-1 Q_k  ......  (Q_1 ... Q_k)^T A (Q_1 .... Q_k)
     // D = Q^T A Q     Q=Q_1 Q_2 ... Q_k
    do
    {
        for(int j=0;j<N-1;j++)//j=1///////////////
        {
            if(j==0)
            {
                for(int i=0;i<N;i++)
                {
                    if(i==0)
                        e[i]=1.0;
                    else
                        e[i]=0.0;
                }
                for(int i=0;i<N;i++)
                    a[i]=A[i][0];/////////////
                for(int i=0;i<N;i++)
                    v[i]=a[i]+(norm(a,N)*e[i]);
                for(int i=0;i<N;i++)
                {
                    for(int k=0;k<N;k++)
                        H[i][k]= I[i][k] - ((2.0/(pow(norm(v,N),2)))*(v[i]*v[k]));
                }
                for(int i=0;i<N;i++)
                {
                    for(int k=0;k<N;k++)
                        Q[i][k]=H[i][k];
                }
            }
            //////////////////////
            else
            {
                for(int i=0;i<N;i++)
                {
                    if(i==j)
                        e[i]=1.0;
                    else
                        e[i]=0.0;
                }
                for(int l=0;l<N;l++)
                {
                    for(int i=0;i<N;i++)
                    {
                        for(int k=0;k<N;k++)
                            sum = sum + (Q[i][k]*A[k][l]);
                        if(fabs(sum)<.001)
                            R[i][l]=0.0;
                        else
                            R[i][l]=sum;
                        sum=0.0;
                    }
                }
                for(int i=0;i<j;i++)
                    a[i]=0.0;
                for(int i=j;i<N;i++)
                    a[i]=R[i][j];
                for(int i=0;i<N;i++)
                    v[i]=a[i]+(norm(a,N)*e[i]);
                for(int i=0;i<N;i++)
                {
                    for(int k=0;k<N;k++)
                        V[i][k]=v[i]*v[k];
                }
                for(int i=0;i<N;i++)
                {
                    for(int k=0;k<N;k++)
                        H[i][k]= I[i][k] - ((2.0/(pow(norm(v,N),2)))*(V[i][k]));
                }
                for(int l=0;l<N;l++)
                {
                    for(int i=0;i<N;i++)
                    {
                        for(int k=0;k<N;k++)
                            sum = sum + (H[i][k]*Q[k][l]);
                        Qtot[i][l] = sum;
                        sum = 0.0;
                    }
                }
                for(int i=0;i<N;i++)
                {
                    for(int k=0;k<N;k++)
                        Q[i][k]=Qtot[i][k];
                }
            }///////////////////////////////////////////////////////////////////////////
        }
        c++;
        
        sum=0.0;
        //FINDING R
        for(int l=0;l<N;l++)
        {
            for(int i=0;i<N;i++)
            {
                for(int j=0;j<N;j++)
                    sum = sum + (Q[i][j]*A[j][l]);
                if(fabs(sum)<.001)
                    R[i][l]=0.0;
                else
                    R[i][l]=sum;
                sum=0.0;
            }
        }
        for(int i=0;i<N;i++)
        {
            for(int k=0;k<N;k++)
                Q[i][k]=Qtot[k][i];
        }
        sum=0.0;
        for(int m=0;m<N;m++)
        {
            for(int o=0;o<N;o++)
            {
                for(int z=0;z<N;z++)
                    sum=sum+ (R[o][z]*Q[z][m]);
                if(fabs(sum)<=.001)
                    A1[o][m]=0.0;
                else
                    A1[o][m]=sum;
                sum=0.0;
            }
        }
        
        //for continuing to find similar matrix, singular values i.e. eigenvalues
        for(int m=0;m<N;m++)
        {
            for(int o=0;o<N;o++)
                A[m][o]=A1[m][o];
        }
        for(int m=0;m<N-1;m++)
        {
            lower[m]=A[m+1][m];
        }
    }while(fabs(norm(lower,(N-1))) >.001);
    
    ///
    
    //SIGMA
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            if(i==j)
                eV[i]=A[i][j];
        }
    }
    for(int i=0;i<N;i++)
    {
        for(int j=i+1;j<N-1;j++)
        {
            if(eV[j]>eV[i])
            {
                reorder=eV[j];
                eV[j]=eV[i];
                eV[i]=reorder;
            }
        }
    }
    for(int i=0;i<N;i++)
            repeat[i]=0;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            if(fabs(eV[i] - eV[j]) < 0.001 && i!=j )
                repeat[i]=j;
        }
    }
    //SIGMA
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            if(i==j)
            {
                sigma[i][j]=eV[i];
            }
            else
                sigma[i][j]=0.0;
        }
    }
    cout << "\nSigma:\n";
    for(int i =0;i<N;i++)
    {
        cout << "|";
        for(int j=0;j<N;j++)
        {
            cout << sigma[i][j];
            if(j<N-1)
                cout <<"\t";
        }
        cout << "|\n";
    }
    //SIGMA INVERSE
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            if(i==j && sigma[i][j]!=0.0)
                sigmainv[i][j]=1.0/sigma[i][j];
            else
                sigmainv[i][j]=sigma[i][j];
        }
    }
    //AA
    sum=0.0;
    for(int l=0;l<N;l++)
    {
        for(int i=0;i<N;i++)
        {
            for(int k=0;k<N;k++)
                sum = sum + (Asave[i][k]*Asave[k][l]);
            if(fabs(sum)<.001)
                AA[i][l]=0.0;
            else
                AA[i][l]=sum;
            sum=0.0;
        }
        sum=0.0;
    }
    //U
    alpha=0.0;
    for(int l=0;l<N;l++)
    {
        zerodiff=0;
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
                eig[i][j]=Asave[i][j] - (eV[l]*I[i][j]);//Asave
        }
        sum=0.0;
        for(int i=0;i<N;i++)
        {
            if(fabs(eig[i][i])<.001)
                zerodiff=1.0;
        }
        
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                if(fabs(eig[i][j])<0.001)
                    eig[i][j]=0.0;
            }
        }
        sum=0.0;
        int stop=0;
        for(int i=1;i<N;i++)
        {
            if(fabs(eig[i][l])==1.0)
            {
                for(int j=0;j<N;j++)
                {
                    temp[j]=eig[l][j];
                    eig[l][j]=eig[i][j];
                    eig[i][j]=temp[j];
                }
                stop=1;
            }
            if(stop==1)
                break;
        }
        stop=0;
        if(eig[l][l]!=0.0)
            alpha=eig[l][l];
        for(int i=0;i<N;i++)
        {
            eig[l][i]=eig[l][i]/alpha;
            if(fabs(eig[l][i])==0.0)
                eig[l][i]=0.0;
        }
        for(int i=0;i<N;i++)
        {
            for(int k=i;k<N;k++)
            {
                if(fabs(eig[i][i])!=0.0)
                    alpha= -1.0*(eig[k+1][i]/eig[i][i]);
                else if(fabs(eig[i][i])==0.0)
                {
                    for(int m=i+1;m<N;m++)
                    {
                        if(fabs(eig[m][i])!=0.0)//==1.0
                        {
                            for(int mm=0;mm<N;mm++)
                            {
                                temp[mm]=eig[i][mm];
                                eig[i][mm]=eig[m][mm];
                                eig[m][mm]=temp[mm];
                            }
                            stop=1;
                        }
                    }
                    if(fabs(eig[i][i])!=0.0)
                    {
                        hold=eig[i][i];
                        for(int m=0;m<N;m++)
                        {
                            eig[i][m]=eig[i][m]/hold;
                            if(fabs(eig[i][m])==0.0)
                                eig[i][m]=0.0;
                        }
                        alpha= -1.0*(eig[k+1][i]/eig[i][i]);
                    }
                }
                
                for(int j=i;j<N;j++)
                {
                    if(k<N-1)
                        eig[k+1][j]=eig[k+1][j] + alpha*eig[i][j];
                    if(fabs(eig[k+1][j])<.001)
                        eig[k+1][j]=0.0;
                }
            }
        }
        /*for(int i=0;i<N;i++)
        {
            for(int j=i;j<N;j++)
            {
                if(i==j && eig[i][j]!=0.0)
                    alpha=eig[i][j];
                eig[i][j]=eig[i][j]/alpha;
                if(fabs(eig[i][j])==0.0)
                    eig[i][j]=0.0;
            }
        }*/
        sum=0.0;
        for(int i=0;i<N;i++)
            sum=sum+eig[i][i];
        if(fabs(sum)==N)
        {
            sum=0.0;
            for(int i=N-1;i>=0;i--)
            {
                if(i==N-1)
                {
                    usvd[i][l]=1.0;
                }
                else
                {
                    for(int j=N-1;j>i;j--)
                        sum=sum+(-1.0)*eig[i][j]*usvd[j][l];
                    usvd[i][l]=sum;
                    sum=0.0;
                }
                if(fabs(usvd[i][l])==0.0)
                    usvd[i][l]=0.0;
                else if(fabs(usvd[i][l])<.001)
                    usvd[i][l]=0.0;
            }
        }
        /////////
        else
        {
            //BACKSUB//
            stop=0;
            for(int i=0;i<N;i++)
            {
                if(fabs(eig[i][i])==0.0 && fabs(eig[i+1][i])!=0.0)
                {
                    for(int j=0;j<N;j++)
                    {
                        temp[j]=eig[i+1][j];
                        eig[i+1][j]=eig[i][j];
                        eig[i][j]=temp[j];
                    }
                    stop=1;
                }
                if(stop==1)
                    break;
            }
            
            for(int i=0;i<N;i++)
            {
                for(int j=0;j<N;j++)
                {
                    if(fabs(eig[i][j])<0.01)
                        eig[i][j]=0.0;
                }
            }
            for(int i=0;i<N;i++)
            {
                if(fabs(eig[i][i])!=1.0 && fabs(eig[i][i])!=0.0)
                {
                    alpha=eig[i][i];
                    for(int j=i;j<N;j++)
                    {
                        eig[i][j]=eig[i][j]/alpha;
                        if(fabs(eig[i][j])<0.001)
                            eig[i][j]=0.0;
                    }
                }
            }
            
            for(int i=0;i<N;i++)
            {
                if(eig[i][i]<0.0)
                {
                    for(int j=0;j<N;j++)
                    {
                        if(fabs(eig[i][j])!=0.0)
                            eig[i][j]=(-1.0)*eig[i][j];
                    }
                }
            }
            
            stop=0;
            beta=0.0;
            for(int i=N-1;i>=0;i--)
            {
                for(int k=i;k>=0;k--)
                {
                    if(eig[i][i]!=0.0)//DIV0!
                    {
                        beta= -1.0*(eig[k-1][i]/eig[i][i]);
                        for(int j=N-1;j>=0;j--)
                        {
                            if(k>0)
                                eig[k-1][j]=eig[k-1][j] + beta*eig[i][j];
                            if(fabs(eig[k-1][j])<.001)
                                eig[k-1][j]=0.0;
                        }
                    }

                }
            }
            
            for(int i=0;i<N;i++)
            {
                for(int j=0;j<N;j++)
                {
                    if(fabs(eig[i][j])<0.01)
                        eig[i][j]=0.0;
                }
            }
            sum=0.0;
            for(int i=N-1;i>=0;i--)
            {
                sum=0.0;
                if(fabs(eig[i][i])==0.0)
                {
                    if(fabs(eV[l+1] - eV[l]) <0.001)
                    {
                        if(zerodiff==1.0 || zerodiff==2.0)
                        {
                            if(zerodiff==1.0)
                                usvd[i][l]=1.0;
                            else
                                usvd[i][l]=0.0;
                            zerodiff=2.0;
                            
                        }
                        else
                            usvd[i][l]=1.0;
                    }
                    else if(fabs(eV[l-1] - eV[l]) <0.001)/////////////////////////////
                    {
                        if(zerodiff==1.0 || zerodiff==2.0)
                        {
                            if(zerodiff==1.0)
                                usvd[i][l]=0.0;
                            else
                                usvd[i][l]=1.0;
                            zerodiff=2.0;
                            
                        }
                        else
                            usvd[i][l]=1.0;
                    }
                    else
                    {
                        if(zerodiff==1.0 || zerodiff==2.0)
                        {
                            if(zerodiff==1.0)
                                usvd[i][l]=1.0;
                            else
                                usvd[i][l]=0.0;
                            zerodiff=2.0;
                            
                        }
                        else
                            usvd[i][l]=1.0;
                    }
                }
                else
                {
                    
                    for(int j=N-1;j>i;j--)
                        sum=sum+(usvd[j][l]*eig[i][j]);
                    if(fabs(sum)<0.001)
                        usvd[i][l]=0.0;
                    else
                        usvd[i][l]=-1.0*sum;
                }
                sum=0.0;
            }
            
        }
    }
    
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
            nrm[j]=usvd[j][i];
        for(int j=0;j<N;j++)
            usvd[j][i]=usvd[j][i]/norm(nrm,N);
    }
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            //if(j!=1 && j!=2 && j!=3 && j!=5 && j!=6 && j!=8)
            // usvd[i][j]=(-1.0*usvd[i][j]);
            if(fabs(usvd[i][j])==0.0)
                usvd[i][j]=0.0;
            else if(fabs(usvd[i][j])<.001)
                usvd[i][j]=0.0;
        }
    }
    //V!
    for(int i =0;i<N;i++)
    {
        for(int j=0;j<N;j++)
            vsvd[i][j]=usvd[j][i];
    }
    sum=0.0;
    /*for(int l=0;l<N;l++)
    {
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
                sum = sum + (sigma[i][j]*vsvd[j][l]);
            if(fabs(sum)<.001)
                vholder[i][l]=0.0;
            else
                vholder[i][l]=sum;
            sum=0.0;
        }
    }*/
    sum=0.0;
    //MOORE PENROSE INVERSE
    sum=0.0;
    for(int l=0;l<N;l++)
    {
        for(int i=0;i<N;i++)
        {
            for(int k=0;k<N;k++)
                sum = sum + (sigmainv[i][k]*vsvd[k][l]);
            if(fabs(sum)<.001)
                tmp[i][l]=0.0;
            else
                tmp[i][l]=sum;
            sum=0.0;
        }
    }
    sum=0.0;
    for(int l=0;l<N;l++)
    {
        for(int i=0;i<N;i++)
        {
            for(int k=0;k<N;k++)
                sum = sum + (usvd[i][k]*tmp[k][l]);
            if(fabs(sum)<.001)
                MPSI[i][l]=0.0;
            else
                MPSI[i][l]=sum;
            sum=0.0;
        }
    }
    //Xhat
    sum=0.0;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
            sum= sum + MPSI[i][j]*b[j];
        if(fabs(sum)<.001)
            x[i]=0.0;
        else
            x[i]=sum;
        sum=0.0;
    }
    //A=USIGMAVT
    sum=0.0;
    for(int l=0;l<N;l++)
    {
        for(int i=0;i<N;i++)
        {
            for(int k=0;k<N;k++)
                sum = sum + (sigma[i][k]*vsvd[k][l]);
            if(fabs(sum)<.001)
                tmp[i][l]=0.0;
            else
                tmp[i][l]=sum;
            sum=0.0;
        }
        sum=0.0;
    }
    sum=0.0;
    for(int l=0;l<N;l++)
    {
        for(int i=0;i<N;i++)
        {
            for(int k=0;k<N;k++)
                sum = sum + (usvd[i][k]*tmp[k][l]);
            if(fabs(sum)<.001)
                Arecons[i][l]=0.0;
            else
                Arecons[i][l]=sum;
            sum=0.0;
        }
        sum=0.0;
    }
    cout << "\nSingular Values:\n";
    for(int i=0;i<N;i++)
        cout << "|" << eV[i] <<"|\n";
    cout <<"\nU:\n";
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
            cout << "|"<< usvd[i][j] << "|\t";
        cout << endl;
    }
    cout << "\nSigma:\n";
    for(int i =0;i<N;i++)
    {
        cout << "|";
        for(int j=0;j<N;j++)
        {
            cout << sigma[i][j];
            if(j<N-1)
                cout <<"\t";
        }
        cout << "|\n";
    }
    cout << "\nV TRANSPOSE:\n";
    for(int i =0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            cout << "|"<< vsvd[i][j] <<"|\t";
        }
        cout << endl;
    }
    cout <<"\nMOORE PENROSE INVERSE:\n";
    for(int i =0;i<N;i++)
    {
        cout << "|";
        for(int j=0;j<N;j++)
        {
            cout << MPSI[i][j];
            if(j<N-1)
                cout <<"\t";
        }
        cout << "|\n";
    }
    cout << "\nx:\n";
    for(int i=0;i<N;i++)
        cout << "|"<<x[i]<<"|\n";
    
    cout << "\nU*Sigma*Vt:\n";
    for(int i =0;i<N;i++)
    {
        cout << "|";
        for(int j=0;j<N;j++)
        {
            cout << Arecons[i][j];
            if(j<N-1)
                cout <<"\t";
        }
        cout << "|\n";
    }
    
    return 0;
}

double norm(double u[], int n)
{
    double uNorm,sum=0.0;
    
    for(int i=0;i<n;i++)
        sum= sum + (u[i]*u[i]);
    uNorm =sqrt(sum);
    return uNorm;
}

double normInf(double u[], int n)
{
    double max=0.0;
    for(int i=0;i<n;i++)
    {
        if(fabs(u[i])>max)
            max=fabs(u[i]);
    }
    return max;
}


































