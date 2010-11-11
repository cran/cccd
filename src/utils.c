#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* sorts smallest to biggest */
int cmp1(const void *x, const void *y)
{
   double *xx,*yy;

   xx = (double *)x;
   yy = (double *)y;

   if((*xx)<(*yy)) return(-1);
   if((*xx)>(*yy)) return(1);
   return(0);
}

