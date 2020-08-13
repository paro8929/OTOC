#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>


using namespace std;

int bigN=2;
int ABORT=0;
int MAX=2;
const int cutM=32;

const double M13=1./3.;
const double M23=2./3.;
const double M43=4./3.;

double * en;
double xn[cutM][cutM];

double Rbeven(int n,int m,double t)
{
  double res=0;
  for (int i=0;i<MAX;i++)
    {
      res+=xn[n][i]*xn[m][i]*((en[2*i+1]-en[2*m])*cos((en[2*i+1]-en[2*n])*t)+(en[2*i+1]-en[2*n])*cos((en[2*i+1]-en[2*m])*t));
      //printf("adding %f where %f %f\n",xn[n][i]*xn[m][i]*((en[2*i+1]-en[2*m])*cos((en[2*i+1]-en[2*n])*t)+(en[2*i+1]-en[2*m])*cos((en[2*i+1]-en[2*m])*t)),xn[n][i]*xn[m][i],(en[2*i+1]-en[2*m]),);
    }
  res*=0.5;
  return res;
}

double Ibeven(int n,int m,double t)
{
  double res=0;
  for (int i=0;i<MAX;i++)
    {
      res+=xn[n][i]*xn[m][i]*(-(en[2*i+1]-en[2*m])*sin((en[2*i+1]-en[2*n])*t)+(en[2*i+1]-en[2*n])*sin((en[2*i+1]-en[2*m])*t));
    }
  res*=0.5;
  return res;
}

double Rbodd(int n,int m,double t)
{
  double res=0;
  for (int i=0;i<MAX;i++)
    {
      res+=xn[i][n]*xn[i][m]*((en[2*i]-en[2*m+1])*cos((en[2*n+1]-en[2*i])*t)+(en[2*i]-en[2*n+1])*cos((en[2*m+1]-en[2*i])*t));
      //printf("adding %f where %f %f\n",xn[i][n]*xn[i][m]*((en[2*i]-en[2*m+1])*cos((en[2*n+1]-en[2*i])*t)+(en[2*i]-en[2*n+1])*cos((en[2*m+1]-en[2*i])*t)),xn[i][n]*xn[i][m],(en[2*i]-en[2*m+1]));
    }
  res*=0.5;
  return res;
}



double Ibodd(int n,int m,double t)
{
  double res=0;
  for (int i=0;i<MAX;i++)
    {
      res+=xn[i][n]*xn[i][m]*((en[2*i]-en[2*m+1])*sin((en[2*n+1]-en[2*i])*t)-(en[2*i]-en[2*n+1])*sin((en[2*m+1]-en[2*i])*t));
      //printf("adding %f where %f %f\n",xn[n][i]*xn[m][i]*((en[2*i+1]-en[2*m])*cos((en[2*i+1]-en[2*n])*t)+(en[2*i+1]-en[2*m])*cos((en[2*i+1]-en[2*m])*t)),xn[n][i]*xn[m][i],(en[2*i+1]-en[2*m]),);
    }
  res*=0.5;
  return res;
}


double ceven(int n,double t)
{
  double res=0;
  for (int i=0;i<MAX;i++)
    {
      double r=Rbeven(n,i,t);
      double ii=Ibeven(n,i,t);
      res+=r*r+ii*ii;
    }
  return res;
}

double codd(int n,double t)
{
  double res=0;
  for (int i=0;i<MAX;i++)
    {
      double r=Rbodd(n,i,t);
      double ii=Ibodd(n,i,t);
      res+=r*r+ii*ii;
    }
  return res;
}

double Z(double T)
{
  double res=0;
  for (int i=0;i<2*MAX+2;i++)
    {
      res+=exp(-en[i]/T);
    }
  return res;
}

double Tinc(int u)
{
  return 0.25*pow(2,u);
}
  

int main(int argc, char *argv[])
{

  printf("INFO:\t This is the OTOC solver for quartic quantum oscillator, v1.0\n");
  printf("\n \t Syntax: ./otoc [K] [MAX] \n");
  printf("\n \t\t [K] ... dimensions of eigenmatrix to be solved\n");
  printf("\n \t\t [MAX] ... max dimension of otoc, e.g. c_{MAX}\n\n");
  
  if (argc>1)
    bigN=stoi(argv[1]);

  if (argc>2)
    MAX=stoi(argv[2]);

  
  en=new double[2*bigN+2];

  //load data
  string name="results/En"+to_string(bigN)+".dat";
  
  ifstream infile;
  infile.open(name.c_str());
  if(infile.fail())
      {
	printf("ERROR: file %s does not exist\n",name.c_str());
	ABORT=1;
      }
  else
    {
      double dummy;
      int a;
      for (int i=0;i<2*bigN+2;i++)
	{
	  infile >> a;
	  infile >> dummy;
	  en[i]=dummy;
	}

      //for (int i=0;i<2*bigN+2;i++)
      //	printf("checking En(%i)=%f\n",i,en[i]);
    }
  infile.close();
  

  
  if(!ABORT)
    {
      name="results/xn"+to_string(bigN)+".dat";
      infile.open(name.c_str());
      if(infile.fail())
	{
	  printf("ERROR: file %s does not exist\n",name.c_str());
	  ABORT=1;
	}
      else
	{
	  double dummy;
	  int a;
	  for (int i=0;i<min(bigN+1,cutM);i++)
	    {
	      infile >> a;
	      for (int j=0;j<min(bigN+1,cutM);j++)
		{
		  infile >> dummy;
		  xn[j][i]=dummy;
		  //if (i==0)
		  //printf("j=%i %f\n",j,xn[j][i]);
		}
	    }
	}
      infile.close();
    }

  if(!ABORT)
    {

      FILE *outfile,*otocfile;
      FILE *Cinf,*SSF;
      
      name="results/cn"+to_string(bigN)+"-MAX"+to_string(MAX)+".dat";
      outfile = fopen(name.c_str(),"w");

      name="results/otoc"+to_string(bigN)+"-MAX"+to_string(MAX)+".dat";
      otocfile = fopen(name.c_str(),"w");

      name="results/oinf"+to_string(bigN)+"-MAX"+to_string(MAX)+".dat";
      Cinf = fopen(name.c_str(),"w");

      name="results/SSF"+to_string(bigN)+"-MAX"+to_string(MAX)+".dat";
      //name="results/SSF"+to_string(bigN)+".dat";
      SSF = fopen(name.c_str(),"w");


      
      fprintf(outfile,"#t\t\t");
      for (int u=0;u<MAX;u++)
	{
	  fprintf(outfile,"c_%i\t\tc_%i\t\t",2*u,2*u+1);
	}
      fprintf(outfile,"\n");

      fprintf(otocfile,"#t\t\t");
      for (int u=0;u<MAX;u++)
	{
	  fprintf(otocfile,"T=%f\t\t",Tinc(u));
	}
      fprintf(otocfile,"\n");


      fprintf(Cinf,"#T\t\tC_T(t=inf)\n");
      double avC[MAX];
      for (int v=0;v<MAX;v++)
	avC[v]=0.;

      
      //printf("check %f %f\n",xn[0][0],Rbeven(0,0,0.5));
      for (int i=0;i<501;i++)
	{
	  double t=i*0.1;
	  fprintf(outfile,"%f\t",t);
	  fprintf(otocfile,"%f\t",t);
	  double otoc[MAX];
	  for (int v=0;v<MAX;v++)
	    otoc[v]=0.;

	  
	  for (int u=0;u<MAX;u++)
	    {
	      double ce=ceven(u,t);
	      double co=codd(u,t);
	      fprintf(outfile,"%f\t%f\t",ce,co);
	      for (int v=0;v<MAX;v++)		
		otoc[v]+=exp(-en[2*u]/Tinc(v))*ce+exp(-en[2*u+1]/Tinc(v))*co;	      
	    }	  
	  fprintf(outfile,"\n");

	  for (int v=0;v<MAX;v++)
	    {
	      otoc[v]/=Z(Tinc(v));
	      fprintf(otocfile,"%f\t\t",otoc[v]);
	    }
	  fprintf(otocfile,"\n");

	  if((i>=200)&&(i<250))
	    {
	      for (int v=0;v<MAX;v++)
		{
		  avC[v]+=otoc[v];
		  //if (v==2)
		  //printf("%i %f\n",i,otoc[v]);
		}
	    }
	  
	    
	  
	  //printf("Z(0.61)=%f\n",Z(0.61));
	  //printf("%f %f %f\n",t,Rbeven(0,2,t),Ibeven(0,2,t));
	  //printf("%f %f %f\n",t,Rbodd(0,2,t),Ibodd(0,2,t));
	  //printf("%f %f %f\n",t,ceven(0,t,2),Rbeven(0,0,t)*Rbeven(0,0,t)+Ibeven(0,0,t)*Ibeven(0,0,t));
	  //printf("%f %f\n",t,Rbeven(0,1,t)*Rbeven(0,1,t)+Ibeven(0,1,t)*Ibeven(0,1,t));
	  //printf("%f %f\n",t,ceven(0,t));
	}

      
      fclose(outfile);
      fclose(otocfile);
      
      for (int v=0;v<MAX;v++)
	avC[v]/=50.;

      for (int v=0;v<MAX;v++)
	fprintf(Cinf,"%f\t %f\n",Tinc(v),avC[v]);
      
      fclose(Cinf);
      
      
      
      name="results/Tcomp"+to_string(bigN)+"-MAX"+to_string(MAX)+".dat";
      outfile = fopen(name.c_str(),"w");
      fprintf(outfile,"#T\t\tZ(T)\t\t<p^2>_T\t\t<x^2>_T\n");

      int xABORT=0;
      double *x2vals;
      x2vals=new double[2*MAX];

      
      name="results/Ax2-"+to_string(bigN)+".dat";
      infile.open(name.c_str());
      if(infile.fail())
	{
	  printf("ERROR: file %s does not exist, will skip <x^2> column...\n",name.c_str());
	  xABORT=1;
	}
      else
	{
	  int a;
	  double dummy;
	  for (int i=0;i<MAX;i++)
	    {
	      infile >> a;
	      infile >> dummy;
	      x2vals[2*i]=dummy;
	      infile >> dummy;
	      x2vals[2*i+1]=dummy;
	      //printf("%i %f %f\n",i,x2vals[2*i],x2vals[2*i+1]);
	    }
	}
      

      
      for (int v=0;v<MAX;v++)
	{
	  double z=Z(Tinc(v));
	  fprintf(outfile,"%f\t%f\t",Tinc(v),z);
	  double dummy=0;
	  for (int u=0;u<MAX;u++)
	    {
	      dummy+=exp(-en[2*u]/Tinc(v))*en[2*u];
	      dummy+=exp(-en[2*u+1]/Tinc(v))*en[2*u+1];
	    }
	  fprintf(outfile,"%f\t",dummy/z);

	  dummy=0;
	  for (int u=0;u<MAX;u++)
	    {
	      dummy+=exp(-en[2*u]/Tinc(v))*x2vals[2*u];
	      dummy+=exp(-en[2*u+1]/Tinc(v))*x2vals[2*u+1];
	    }
	  fprintf(outfile,"%f\t",dummy/z);
	  
	  fprintf(outfile,"\n");
	}

      fclose(outfile);


      fprintf(SSF,"#t\t\t");
      for (int v=0;v<MAX;v++)
	fprintf(SSF,"%f\t\t",Tinc(v));
      fprintf(SSF,"\n");
      
      for (int i=0;i<1001;i++)
	{
	  double t=i*0.01;
	  fprintf(SSF,"%f\t",t);
	  for (int v=0;v<MAX;v++)
	    {
	      double beta=1/Tinc(v);
	      double z=Z(Tinc(v));
	      double tbI=0;
	      double tbR=0;
	      //double t0I=0;
	      //double t0R=0;
	      for (int u=0;u<bigN+1;u++)
		{
		  tbR+=exp(-beta*en[2*u])*cos(en[2*u]*t);
		  tbI+=exp(-beta*en[2*u])*sin(en[2*u]*t);
		  //t0R+=cos(en[2*u]*t);
		  //t0I+=sin(en[2*u]*t);
		  tbR+=exp(-beta*en[2*u+1])*cos(en[2*u+1]*t);
		  tbI+=exp(-beta*en[2*u+1])*sin(en[2*u+1]*t);
		  //t0R+=cos(en[2*u+1]*t);
		  //t0I+=sin(en[2*u+1]*t);
		}
	      //double ssf=(tbR*tbR+tbI*tbI)/(t0R*t0R+t0I*t0I);
	      double ssf=(tbR*tbR+tbI*tbI)/(z*z);
	      fprintf(SSF,"%f\t",ssf);
	    }
	  fprintf(SSF,"\n");
	}
      

      
      fclose(SSF);
      
      delete [] x2vals;
    }

  
  delete [] en;

}
