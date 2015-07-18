#include <stdio.h>
#include "lapacke.h"

int main (int argc, const char * argv[])
{
   double a[3][3] = {1,1,1,2,3,4,3,5,2};
   lapack_int info,m,n,lda;
   lapack_int ipiv[3];
   int i,j;

   m = 3;
   n = 3;
   lda = 3;

   info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,m,n, *a, lda, ipiv);
	info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, *a, lda, ipiv);

   for(i=0;i<n;i++)
   {
      for(j=0;j<n;j++)
      {
         printf("%lf ",a[i][j]);
      }
      printf("\n");
   }
   return(info);
}

