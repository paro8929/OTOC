#include <string>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <math.h> 

using namespace std;

int bigN=2;
int force=0;//flag to force overwrite


const double M13=1./3.;
const double M23=2./3.;
const double M43=4./3.;



double pochhammer(double x,double n)
{
  return tgamma(x+n)/tgamma(x);
}

int si(int i)
{
  if ((i%2)==0)
    return 1;
  else
    return -1;
}

double Hyp1(int n, int m)
{
  double sum=0;
  if (m>n)
    {
       int z=m;
       m=n;
       n=z; 
    }
//now n>=m
  if (n<=50)
    {
      for(int i=0;i<=m;i++)
	{
	  sum+=pochhammer(n-i+1,i)*pochhammer(m-i+1,i)*pochhammer(M13,i)/pochhammer(M23-n,i)/pochhammer(M23-m,i)/tgamma(i+1);	  
	}
    }
  else 
    {      
      for(int i=0;i<=m;i++)
	sum+=exp(lgamma(n+1)-lgamma(n-i+1)+lgamma(m+1)-lgamma(m-i+1)+lgamma(n+1-i-M23)-lgamma(1+n-M23)+lgamma(m+1-i-M23)-lgamma(1+m-M23)+lgamma(M13+i)-lgamma(i+1));
      sum/=tgamma(M13);
    }
  return sum;
}


double Hyp2(int n, int m)
{
  double sum=0;
  if (m>n)
    {
       int z=m;
       m=n;
       n=z;
    }
//now n>=m
  if (n<=50)
    {
      for(int i=0;i<=m;i++)
	{
	  sum+=pochhammer(n-i+1,i)*pochhammer(m-i+1,i)/pochhammer(M23-n,i)/pochhammer(M23-m,i);
	}
    }
  else
    {
      for(int i=0;i<=m;i++)
	sum+=exp(lgamma(n+1)-lgamma(n-i+1)+lgamma(m+1)-lgamma(m-i+1)+lgamma(n+1-i-M23)-lgamma(1+n-M23)+lgamma(m+1-i-M23)-lgamma(1+m-M23));
    }
  return sum;
}


double Aeven(int n,int m)
{
   double res=0;
   if ((n<=45)&&(m<=45))
     res=pochhammer(M13,n)*pochhammer(M13,m)*Hyp1(n,m)/sqrt(9*(M13+n)*(M13+m)*tgamma(n+1)*tgamma(m+1)*tgamma(M23+n)*tgamma(M23+m));
   else
     res=exp(lgamma(M13+n)+lgamma(M13+m)-0.5*lgamma(n+1)-0.5*lgamma(m+1)-0.5*lgamma(M23+n)-0.5*lgamma(M23+m))*Hyp1(n,m)/(3*tgamma(M13)*tgamma(M13)*sqrt((M13+n)*(M13+m)));
   return res;
}



double Aodd(int n,int m)
{
   double res=0;
   if ((n<=45)&&(m<=45))
     res=pochhammer(M13,n)*pochhammer(M13,m)*Hyp2(n,m)/sqrt(9*(M23+n)*(M23+m)*tgamma(n+1)*tgamma(m+1)*tgamma(M43+n)*tgamma(M43+m));
   else
     {
       res=exp(lgamma(M13+n)+lgamma(M13+m)-0.5*lgamma(n+1)-0.5*lgamma(m+1)-0.5*lgamma(M43+n)-0.5*lgamma(M43+m))*Hyp2(n,m)/(3*tgamma(M13)*tgamma(M13)*sqrt((M23+n)*(M23+m)));
     }
   return res;
}


int main(int argc, char *argv[])
{
 

  if (argc>1)
    bigN=stoi(argv[1]);

  if (argc>2)
    force=1;

  //printf("check %.16g %.16g\n",Aeven(50,64),Hyp1(50,64));
  
  //printf("Ghy %f\n",GHyp1(256,171));
  //for (int i=0;i<256;i++)
  //  printf("%i check %.16g %.16g \n",i,Aeven(i,64),Aodd(i,64));	

  //printf("check %.16g %.16g\n",Aeven(67,99),Aeven(67,100));

  printf("INFO: Allocating memory...\n");
  
  //double mat[bigN+1][bigN+1];

  double **mat;
  mat=new double*[bigN+1];
  for (int i=0;i<bigN+1;i++)
    mat[i]=new double[bigN+1];
  
  printf("\t...done.\n");

  string name="matrices/even"+to_string(bigN)+".dat";
  ifstream infile;
  infile.open(name);
  if ((infile.fail())||(force))
  {

    printf("INFO: Generating even matrix N=%i...\n",bigN);
    infile.close();
    FILE *outfile;
    outfile = fopen(name.c_str(),"w");
    for (int i=0;i<bigN+1;i++)
      for (int j=i;j<bigN+1;j++)
        {
          mat[i][j]=Aeven(i,j);
          mat[j][i]=mat[i][j];
	  if (mat[i][j]<1.e-10)
	    printf("i=%i j=%i %.16g\n",i,j,mat[i][j]);
        }

    for (int i=0;i<bigN+1;i++)
      {
        for (int j=0;j<bigN+1;j++)
          fprintf(outfile,"%.16g\t",mat[i][j]);
        fprintf(outfile,"\n"); 
      } 
    fclose(outfile);
  }
  else printf("ERROR: file %s already exists, nothing to do \n",name.c_str());
  infile.close();

  printf("\t...done.\n");

  
  name="matrices/odd"+to_string(bigN)+".dat";
  infile.open(name);
  if ((infile.fail())||(force))
  {

    printf("INFO: Generating odd matrix N=%i...\n",bigN);
    infile.close();
    FILE *outfile;
    outfile = fopen(name.c_str(),"w");
    for (int i=0;i<bigN+1;i++)
      for (int j=i;j<bigN+1;j++)
        {
          mat[i][j]=Aodd(i,j);
          mat[j][i]=mat[i][j];
        }

    for (int i=0;i<bigN+1;i++)
      {
        for (int j=0;j<bigN+1;j++)
          fprintf(outfile,"%.16g\t",mat[i][j]);
        fprintf(outfile,"\n");
      }
    fclose(outfile);
  }
  else printf("ERROR: file %s already exists, nothing to do \n",name.c_str());
  infile.close();

  printf("\t...done.\n");
  
   printf("INFO: Freeing memory...\n");


  for (int i=0;i<bigN+1;i++)
    delete [] mat[i];

  delete [] mat;
  
  printf("\t...done.\n");
}
