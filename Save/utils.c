#include <stdio.h>
#include <math.h>

#define MIN(x,y) (((x)<(y))?(x):(y))

/* sorts biggest to smallest */
int cmp(const void *x, const void *y)
{
   double **xx,**yy;

   xx = (double **)x;
   yy = (double **)y;

   if((*xx)[0]>(*yy)[0]) return(-1);
   if((*xx)[0]<(*yy)[0]) return(1);
   return(0);
}

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

/* sorts smallest to biggest */
int cmp2(const void *x, const void *y)
{
   double **xx,**yy;

   xx = (double **)x;
   yy = (double **)y;

   if((*xx)[0]<(*yy)[0]) return(-1);
   if((*xx)[0]>(*yy)[0]) return(1);
   return(0);
}

