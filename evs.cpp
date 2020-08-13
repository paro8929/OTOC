#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
 
using namespace Eigen;
using namespace std;

int bigN=2;
const double M13=1./3.;
const double M23=2./3.;
const double M43=4./3.;

int cutM=32;
int generate_x2_flag=0; //if set, this flag will start generating the calculation of <n|x^2|n>, which is numerically costly


MatrixXd test(const char *filename)
{
  MatrixXd result(bigN+1,bigN+1);

  
  ifstream infile;
    infile.open(filename);
    if(infile.fail())
      {
	printf("ERROR: file %s does not exist\n",filename);
      }
    else
      {
	
	double dummy;
	for (int row=0;row<bigN+1;row++)
	  {
	    for(int col=0;col<bigN+1;col++)
	      {
		infile >> dummy;
		result(row,col)=dummy;		
		//printf("%.16g\n",dummy);
	      }
	    //printf("Row %i done\n",row);
	  }	  	
      }
    return result;
}


double pochhammer(double x,double n)
{
  return tgamma(x+n)/tgamma(x);
}


//assume n is 2n and m is 2m+1 
double sumxnm(int n,int m,MatrixXd m_even,MatrixXd m_odd)
{
  MatrixXd V_even = MatrixXd::Random(bigN+1,1);
  MatrixXd V_odd = MatrixXd::Random(bigN+1,1);
  
  V_even=m_even.col(n);
  V_odd=m_odd.col(m);

  //std::cout << "Here is normalized V_even:\n" << V_even << std::endl;

  
  double norm=0;
  for (int i=0;i<bigN+1;i++)
    for (int j=0;j<bigN+1;j++)
      {
	if ((i<40)&&(j<40))
	  norm+=2*V_even(i,0)*V_odd(j,0)*(1.-3*i)*(1.+3*j)/(1.-3*i+3*j)*pochhammer(-1./3.,i)*pochhammer(1./3.,j)/sqrt((1+3*i)*(2+3*j)*tgamma(i+2./3.)*tgamma(j+4./3.)*tgamma(i+1)*tgamma(j+1));
	else
	  norm+=2*V_even(i,0)*V_odd(j,0)*(1.-3*i)*(1.+3*j)/(1.-3*i+3*j)*exp(lgamma(i-M13+1)+lgamma(M13+j)-0.5*lgamma(i+M23)-0.5*lgamma(j+M43)-0.5*lgamma(i+1)-0.5*lgamma(j+1))/((i-M13)*tgamma(M13)*tgamma(-M13)*sqrt((1+3*i)*(2+3*j)));	
      }

  return norm;
}


double Hyp3(int n, int m)
{
  double sum=0;
  if (m>n)
    {
       int z=m;
       m=n;
       n=z; 
    }
//now n>=m
  for(int i=0;i<=m;i++)
    if ((i<=40))
      sum+=pochhammer(n-i+1,i)*pochhammer(m-i+1,i)/pochhammer(M43-n,i)/pochhammer(M43-m,i);
    else
      sum+=exp(lgamma(n+1)+lgamma(m+1)-lgamma(n-i+1)-lgamma(m+1-i)-lgamma(1+n-M13)+lgamma(n-i+M23)-lgamma(1+m-M13)+lgamma(m-i+M23))*(M13-n)*(M13-m)/(M13-n+i)/(M13-m+i);
  return sum;
}

/*
double AHyp3(int n, int m)
{
  double sum=0;
  if (m>n)
    {
       int z=m;
       m=n;
       n=z; 
    }
//now n>=m
  for(int i=0;i<=m;i++)    
    sum+=exp(lgamma(n+1)+lgamma(m+1)-lgamma(n-i+1)-lgamma(m+1-i)-lgamma(1+n-M13)+lgamma(n-i+M23)-lgamma(1+m-M13)+lgamma(m-i+M23))*(M13-n)*(M13-m)/(M13-n+i)/(M13-m+i);
  return sum;
  }*/



double Hyp4(int n, int m)
{
  double sum=0;
  if (m>n)
    {
       int z=m;
       m=n;
       n=z; 
    }
//now n>=m
 
      for(int i=0;i<=m;i++)
	{
	  if (i<40)
	    sum+=pochhammer(n-i+1,i)*pochhammer(m-i+1,i)*pochhammer(1+M23,i)/pochhammer(M43-n,i)/pochhammer(M43-m,i)/tgamma(i+1);
	  else
	    {
	      sum+=exp(+lgamma(n+1)-lgamma(n-i+1)+lgamma(m+1)-lgamma(m-i+1)+lgamma(1+M23+i)-lgamma(i+1)-lgamma(M23+n)+lgamma(M23+n-i)-lgamma(M23+m)+lgamma(M23+m-i))*(M13-n)*(M13-m)/(M13-n+i)/(M13-m+i)/tgamma(1+M23);
	    }
	}
	 
  return sum;
}

/*
double AHyp4(int n, int m)
{
  double sum=0;
  if (m>n)
    {
       int z=m;
       m=n;
       n=z; 
    }
//now n>=m
  for(int i=0;i<=m;i++)
    //sum+=exp(+lgamma(n+1)-lgamma(n-i+1)+lgamma(m+1)-lgamma(m-i+1)+lgamma(1+M23+i)-lgamma(i+1))*(M13-n)*(M13-m)/(M13-n+i)/(M13-m+i)*tgamma(M13-n)/tgamma(M13-n+i)*tgamma(M13-m)/tgamma(M13-m+i);
    // sum+=exp(+lgamma(n+1)-lgamma(n-i+1)+lgamma(m+1)-lgamma(m-i+1)+lgamma(1+M23+i)-lgamma(i+1)-lgamma(M23+n)+lgamma(M23+n-i))*(M13-n)*(M13-m)/(M13-n+i)/(M13-m+i)*sin(M_PI*(M13-n+i))/sin(M_PI*(M13-n))*tgamma(M13-m)/tgamma(M13-m+i);
    sum+=exp(+lgamma(n+1)-lgamma(n-i+1)+lgamma(m+1)-lgamma(m-i+1)+lgamma(1+M23+i)-lgamma(i+1)-lgamma(M23+n)+lgamma(M23+n-i)-lgamma(M23+m)+lgamma(M23+m-i))*(M13-n)*(M13-m)/(M13-n+i)/(M13-m+i);
    
  sum/=tgamma(1+M23);
  return sum;
}
*/

double sumx2even(int n,MatrixXd m_even)
{
  MatrixXd V_even = MatrixXd::Random(bigN+1,1);
  
  V_even=m_even.col(n);

  //std::cout << "Here is normalized V_even:\n" << V_even << std::endl;


  double norm=0;
  for (int i=0;i<bigN+1;i++)
    {
      for (int j=0;j<bigN+1;j++)
	{
	  if ((i<=40)&&(j<=40))
	    norm+=2*V_even(i,0)*V_even(j,0)*pochhammer(-1./3.,i)*pochhammer(-1./3.,j)/sqrt((1+3*i)*(1+3*j)*tgamma(i+M23)*tgamma(j+M23)*tgamma(i+1)*tgamma(j+1))*Hyp3(i,j);
	  else
	    norm+=2*V_even(i,0)*V_even(j,0)*exp(0.5*lgamma(M23+i)+0.5*lgamma(M23+j)-0.5*lgamma(i+1)-0.5*lgamma(j+1))/tgamma(-M13)/tgamma(-M13)/((i-M13)*(j-M13)*sqrt((1+3*i)*(1+3*j)))*Hyp3(i,j);
	}
    }

  return norm;
}

double sumx2odd(int n,MatrixXd m_odd)
{
  MatrixXd V_odd = MatrixXd::Random(bigN+1,1);
  
  V_odd=m_odd.col(n);

  double norm=0;
  for (int i=0;i<bigN+1;i++)
    for (int j=0;j<bigN+1;j++)
      {
	if ((i<=40)&&(j<=40))
	  norm+=2*V_odd(i,0)*V_odd(j,0)*pochhammer(-1./3.,i)*pochhammer(-1./3.,j)/sqrt((2+3*i)*(2+3*j)*tgamma(i+M43)*tgamma(j+M43)*tgamma(i+1)*tgamma(j+1))*Hyp4(i,j);
	else
	  {
	    norm+=2*V_odd(i,0)*V_odd(j,0)*exp(lgamma(i+M23)+lgamma(j+M23)-0.5*lgamma(i+M43)-0.5*lgamma(j+M43)-0.5*lgamma(i+1)-0.5*lgamma(j+1))/((i-M13)*(j-M13)*tgamma(-M13)*tgamma(-M13)*sqrt((2+3*i)*(2+3*j)))*Hyp4(i,j);
	  }	  
      }
  norm*=tgamma(5./3.)*pow(3,M23);

  return norm;
}





int main(int argc, char *argv[])
{

  printf("INFO:\t This is the eigensystem solver for quartic quantum oscillator, v1.0\n");
  printf("\n \t Syntax: ./evs [K] [y/n] \n");
  printf("\n \t\t [K] ... dimensions of eigenmatrix to be solved\n");
  printf("\n \t\t [y/n] ... generate <x^2> values (warning: may take a long time!)\n\n");
  
  //printf("checkme %f %f\n",Aeven(3,4),Aodd(5,6));  

  if (argc>1)
    bigN=stoi(argv[1]);
  if (argc>2)
    generate_x2_flag=1;

  
  MatrixXd m_even = MatrixXd::Random(bigN+1,bigN+1);
  MatrixXd m_odd = MatrixXd::Random(bigN+1,bigN+1);

  printf("INFO: Reading in Matrices...\n");
  
  //MatrixXd evenV = MatrixXd::Random(bigN+1,bigN+1);

  //string vav=std::to_string(bigN);
  string name="matrices/even"+to_string(bigN)+".dat";
  //cout << name << endl;
  //test(name.c_str());
  
  //m_even=readMatrix(name.c_str());
   m_even=test(name.c_str());
  name="matrices/odd"+to_string(bigN)+".dat";
  //cout << name << endl;
  m_odd=test(name.c_str());

  printf("\t...done.\n");

  //std::cout << "Here is the matrix m:\n" << m_even << std::endl;

  //EigenSolver<MatrixXd> es(m, false);

  //cout << "The eigenvalues of the 3x3 matrix of ones are:" 
  //  << endl << es.eigenvalues() << endl;

  SelfAdjointEigenSolver<MatrixXd> es_even,es_odd;
  printf("INFO: Solving parity-even eigensystem...\n");
  es_even.compute(m_even);
  printf("\t...done.\n");
  printf("INFO: Solving parity-odd eigensystem...\n");
  es_odd.compute(m_odd);
  printf("\t...done.\n");

  m_even=es_even.eigenvectors();
  m_odd=es_odd.eigenvectors();

  MatrixXd E_even = MatrixXd::Random(bigN+1,1);
  MatrixXd E_odd = MatrixXd::Random(bigN+1,1);

  
  E_even=es_even.eigenvalues();
  E_odd=es_odd.eigenvalues();

  FILE *outfile;
  name="results/En"+to_string(bigN)+".dat";
  outfile = fopen(name.c_str(),"w");
  for (int i=0;i<bigN+1;i++)
    {
      E_even(i,0)=pow(3,M13)/tgamma(M13)/E_even(i,0);
      E_odd(i,0)=pow(3,M13)/E_odd(i,0);      
      //cout << E_even(i,0) << "\t" << E_odd(i,0) << endl;
    }
  for (int i=0;i<bigN+1;i++)
    {
      fprintf(outfile,"%i\t%.16g\n",2*i,E_even(bigN-i,0));
      fprintf(outfile,"%i\t%.16g\n",2*i+1,E_odd(bigN-i,0));
    }
  
  fclose(outfile);

  // std::cout << "Here is V_even:\n" << m_even.col(0) << std::endl;
  
  //properly normalize wave-functions
  
  for (int i=0;i<bigN+1;i++)
    {
      m_even.col(i)*=sqrt(pow(3,M13)*E_even(i)*0.5);
      m_odd.col(i)*=sqrt(E_odd(i)*0.5/pow(3,M13));
    }

  name="results/Vn"+to_string(bigN)+".dat";
  outfile = fopen(name.c_str(),"w");
  fprintf(outfile,"#n\t");
  for (int i=0;i<2*bigN+2;i++)
    fprintf(outfile,"c^{(%i)}\t\t\t",i);
  fprintf(outfile,"\n");
  
  for (int i=0;i<bigN+1;i++)
    {
      fprintf(outfile,"%i\t",i);
      for (int j=0;j<bigN+1;j++)
	{
	  fprintf(outfile,"%.16g\t",m_even(i,bigN-j));
	  fprintf(outfile,"%.16g\t",m_odd(i,bigN-j));
	}
      fprintf(outfile,"\n");
    }
  
  fclose(outfile);

  

  printf("INFO: Generating matrix elements x_{nm}...\n");
  double outmat[cutM][cutM];
  
  
  //printf("%f vs %f \n",Hyp3(3,4),Hyp4(3,4));
  
  
#pragma omp parallel shared(outmat)
  {
#pragma omp for schedule(dynamic,1)
    for (int j=0;j<min(bigN+1,cutM);j++)
      {
	//fprintf(outfile,"%i\t",j);
	for(int i=0;i<min(bigN+1,cutM);i++)
	  {
	    double mm=sumxnm(bigN-i,bigN-j,m_even,m_odd);
	    outmat[j][i]=mm;
	    //fprintf(outfile,"%.16g\t",mm);
	    //printf("Matrix element %i %i =%f\n",i,j,sumxnm(i,j,m_even,m_odd));
	  }
	//fprintf(outfile,"\n");
      }
    //fclose(outfile);
  }

  name="results/xn"+to_string(bigN)+".dat";
  outfile = fopen(name.c_str(),"w");
  for (int j=0;j<min(bigN+1,cutM);j++)
      {
	fprintf(outfile,"%i\t",j);
	for(int i=0;i<min(bigN+1,cutM);i++)
	  fprintf(outfile,"%.16g\t",outmat[j][i]);
	fprintf(outfile,"\n");
      }
  fclose(outfile);
  
    
  printf("\t...done.\n");

  if (generate_x2_flag)
    {
      
      printf("INFO: Generating <n|x^2|n>...\n");
#pragma omp parallel shared(outmat)
      {
#pragma omp for schedule(dynamic,1)
	for (int j=0;j<min(bigN+1,cutM);j++)
	  {
	    double me=sumx2even(bigN-j,m_even);
	    outmat[j][0]=me;
	    me=sumx2odd(bigN-j,m_odd);
	    outmat[j][1]=me;
	  }
      }
      
      name="results/Ax2-"+to_string(bigN)+".dat";
      outfile = fopen(name.c_str(),"w");
      for (int j=0;j<min(bigN+1,cutM);j++)
	fprintf(outfile,"%i\t %.16g\t %.16g\n",j,outmat[j][0],outmat[j][1]);
      fclose(outfile);
      
      //printf("%f vs %f \n",sumx2even(bigN,bigN,m_even),Hyp4(3,4));
      printf("\t...done.\n");
    }

  //printf("%f %f\n",Hyp3(100,10),Hyp4(100,10));
  //printf("%f vs %f\n",Hyp4(50,50),AHyp4(50,50));
  
  // for (int i=0;i<bigN+1;i++)
  //  for (int j=0;j<bigN+1;j++)
  //    printf("Matrix element %i %i =%f\n",i,j,sumxnm(i,j,m_even,m_odd));
  
  
  //evenV=es.eigenvectors();
  //cout << "The eigenvalues are: " << es.eigenvalues().transpose() << endl;
  //cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors() << endl << endl;
}
