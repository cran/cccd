#include <stdio.h>
#include <math.h>

extern void classify_andromeda_class(double *D, int *N, int *CLASS, double *R, int *NW,
                                     int *out, int *doubles, double *MAXR);
extern void andromedaCabP(double *D,double *P, int *N, int *NC, double *ALPHA, double *BETA,
                          int *BYCENTER, int *A, double *R, int *CLASS, int *C, int *NW, 
						  int *VERBOSE);

/* 
 * D is the training distance matrix, ND is the test distance matrix.
 */
void complexity_andromeda_class(
               double *D,      /* distances between observations: nXnw */
               double *ND,     /* distances between observations: mXnw */
			   int *N,         /* number of x observations */
			   int *M,         /* number of y observations */
			   double *P,      /* Probabilities (markings) for the observations */
			   int *CLASS,     /* Classes associated with the observations */
			   int *NC,        /* number of classes */
			   int *BYCENTER,  /* if true, use center of ball for class assignment */
			   double *MAXR,   /* if d/r>maxr no call is made */
			   double *ALPHAS, /* grid of alpha values: length na */
			   int *NA,        /* number of alphas */
			   double *BETAS,  /* grid of beta values: length nb */
			   int *NB,        /* number of betas */
			   int  *ERRS,     /* matrix of errors (returned): naXnb */
			   int  *SIZES)    /* matrix of number of witness elts: naXnb */
{
   int i,j,k;
   int n=*N;                    /* number of training observations */
   int nc=*NC;                  /* number of classes */
   int m=*M;                    /* number of testing observations */
   int na=*NA;
   int nb=*NB;
   double maxr=*MAXR;
   int *A;
   double *R;
   int *OCLASS;
   int nw;
   int *out;
   int temp;
   int verbose=0;

   A = (int *)calloc(n,sizeof(int));
   R = (double *)calloc(n,sizeof(double));
   OCLASS = (int *)calloc(n,sizeof(int));
   out = (int *)calloc(n,sizeof(int));

   for(i=0;i<na;i++){
      for(j=0;j<nb;j++){
	     andromedaCabP(D,P,N,NC,&ALPHAS[i],&BETAS[j],BYCENTER,A,R,OCLASS,
		               (int *)NULL,&nw,&verbose);
		 SIZES[i*nb+j] = nw;
		 classify_andromeda_class(ND,M,OCLASS,R,&nw,out,(int *)NULL,MAXR);
		 temp = 0;
		 for(k=0;k<n;k++){
			if(CLASS[k]!=out[k]){
			   temp ++;
			}
		 }
		 ERRS[i*nb+j] = temp;
	  }
   }

   free(A);
   free(R);
   free(OCLASS);
}
