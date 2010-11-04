#include <stdio.h>
#include <math.h>

#define PRECISION 1E-8

void classify_andromeda_class(
               double *D,     /* distances between observations and witness elts: nXnw       */
			   int *N,        /* number of x observations                                    */
			   int *CLASS,    /* class assigned to the observations: vector of length nw     */
			   double *R,     /* radii: vector of length nw                                  */
			   int *NW,       /* number of witness elements                                  */
			   int *out,      /* output class assignments (-1 means no call)                 */
			   int *doubles,  /* output flag for observations in more than one class' balls  */
			   double *MAXR)  /* if d/r>maxr no call is made                                 */
{
   int i,j;
   int n=*N;                    /* number of observations */
   int nw=*NW;                  /* number of witnesses    */
   double **dists;              /* D as an nxn matrix     */
   double *d;                   /* temporary pointer      */
   double maxr=*MAXR;
   int ind;
   double minv;


   if(doubles) memset(doubles,0,n*sizeof(int));
   d = D;
   for(i=0;i<n;i++){
	  minv = d[0]/R[0];
	  ind = 0;
	  for(j=1;j<nw;j++){
	     if(minv>d[j]/R[j]){
		    minv = d[j]/R[j];
			ind = j;
		 }
	  }
	  if(minv<maxr){
		 out[i] = CLASS[ind];
	  }
	  else{
	     out[i] = -1;
		 /*
		 fprintf(stderr,
			"\noutside: i=%d, ind=%d, class[ind]=%d, minv=%g",i,ind,CLASS[ind],minv);
	     */
	  }
	  if(doubles){
	     for(j=0;j<nw;j++){
			if(CLASS[j] != CLASS[ind]){
			   if(d[j]<R[j]){
				  if((R[j]-d[j])>PRECISION){ /* kludge for round-off errors */
					 doubles[i]++;
					 /*
					 fprintf(stderr,
					"\nj=%d, ind=%d, class[j]=%d, class[ind]=%d, d=%4.16g, r=%4.16g, d/r=%g, d<r=%d",
						j,ind,CLASS[j],CLASS[ind],d[j],R[j],d[j]/R[j],(int)(d[j]<R[j]));
				     */
				  }
			   }
			}
		 }
	  }
	  d += nw;
   }
}
#undef PRECISION
