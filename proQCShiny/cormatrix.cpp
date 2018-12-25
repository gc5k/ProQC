#include <Rcpp.h>
using namespace Rcpp;
#include <string>
#include <vector>
using namespace std;

// [[Rcpp::export]] 
NumericMatrix CorMatrix(NumericMatrix Mat)
{
  int i=0,j=0,k=0;
  int num;
  double sum,temp;
  int c=Mat.ncol();
  int r=Mat.nrow();
  double* m=(double*) malloc (sizeof(double) * c*r);
  memset(m, 0, sizeof(double)* c*r);
  for(i=0; i<r;i++)
  {
    for(j=0; j<c; j++)
    {
      m[i*c+j]=Mat(i,j);
    }
  }
  NumericMatrix exCor(r,r);
  for(i=0; i<r;i++)
  {
    for(j=0; j<=i; j++)
    {
      sum=0.0;
      num=0;
      for(k=0; k<c; k++)
      {
        //temp = Mat(i,k)*Mat(j,k);
        temp = m[i*c+k]*m[j*c+k];
        if(abs(temp)>1.0e-8) {num++;sum += temp;}
      }
      exCor(j,i) = exCor(i,j) = sum/(double)num;
    }
  }
  return wrap(exCor);
}


/***
 hello("World");
 */

