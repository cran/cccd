#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>

void pdist(double *x,int *N, int *D, double *P, double *w, double *dist)
{
   int i,j,k;
   int n=*N,d=*D;
   double p=*P;

   if(w != NULL){
	  for(i=0;i<n;i++){
		 for(j=i;j<n;j++){
			dist[i*n+j] = 0.0;
			for(k=0;k<d;k++){
			   dist[i*n+j] += pow(w[k]*fabs(x[i*d+k]-x[j*d+k]),p);
			}
			dist[i*n+j] = pow(dist[i*n+j],1./p);
			dist[j*n+i] = dist[i*n+j];
		 }
	  }
   }
   else {
	  for(i=0;i<n;i++){
		 for(j=i;j<n;j++){
			dist[i*n+j] = 0.0;
			for(k=0;k<d;k++){
			   dist[i*n+j] += pow(fabs(x[i*d+k]-x[j*d+k]),p);
			}
			dist[i*n+j] = pow(dist[i*n+j],1./p);
			dist[j*n+i] = dist[i*n+j];
		 }
	  }
   }
}

void pdistcos(double *x,int *N, int *D, double *dist)
{
   int i,j,k;
   int n=*N,d=*D;
	double l1,l2;

  for(i=0;i<n;i++){
	 for(j=i;j<n;j++){
		dist[i*n+j] = 0.0;
		l1 = l2 = 0.0;
		for(k=0;k<d;k++){
			dist[i*n+j] += x[i*d+k]*x[j*d+k];
			l1 += x[i*d+k]*x[i*d+k];
			l2 += x[j*d+k]*x[j*d+k];
		}
		if(l1>0 && l2>0){
			dist[i*n+j] /= sqrt(l1*l2);
		}
		dist[i*n+j] = 1.0-dist[i*n+j];
		dist[j*n+i] = dist[i*n+j];
	 }
  }
}

void pdistxy(double *x,double *y,int *Nx, int *Ny, int *D, double *P, 
             double *w, double *dist)
{
   int i,j,k;
   int nx=*Nx,ny=*Ny,d=*D;
   double p=*P;

   if(w != NULL){
	  for(i=0;i<nx;i++){
		 for(j=0;j<ny;j++){
			dist[i*ny+j] = 0.0;
			for(k=0;k<d;k++){
			   dist[i*ny+j] += pow(w[k]*fabs(x[i*d+k]-y[j*d+k]),p);
			}
			dist[i*ny+j] = pow(dist[i*ny+j],1./p);
		 }
	  }
   }
   else {
	  for(i=0;i<nx;i++){
		 for(j=0;j<ny;j++){
			dist[i*ny+j] = 0.0;
			for(k=0;k<d;k++){
			   dist[i*ny+j] += pow(fabs(x[i*d+k]-y[j*d+k]),p);
			}
			dist[i*ny+j] = pow(dist[i*ny+j],1./p);
		 }
	  }
   }
}

#define MAX(x,y) ((x)>(y)?(x):(y))

void pdistinf(double *x,int *N, int *D, double *w, double *dist)
{
   int i,j,k;
   int n=*N,d=*D;

   if(w != NULL){
	  for(i=0;i<n;i++){
		 for(j=i;j<n;j++){
			dist[i*n+j] = 0.0;
			for(k=0;k<d;k++){
			   dist[i*n+j] = MAX(dist[i*n+j],(w[k]*fabs(x[i*d+k]-x[j*d+k])));
			}
			dist[j*n+i] = dist[i*n+j];
		 }
	  }
   }
   else {
	  for(i=0;i<n;i++){
		 for(j=i;j<n;j++){
			dist[i*n+j] = 0.0;
			for(k=0;k<d;k++){
			   dist[i*n+j] = MAX(dist[i*n+j],fabs(x[i*d+k]-x[j*d+k]));
			}
			dist[j*n+i] = dist[i*n+j];
		 }
	  }
   }
}

void pdistxyinf(double *x,double *y,int *Nx, int *Ny, int *D, 
             double *w, double *dist)
{
   int i,j,k;
   int nx=*Nx,ny=*Ny,d=*D;

   if(w != NULL){
	  for(i=0;i<nx;i++){
		 for(j=0;j<ny;j++){
			dist[i*ny+j] = 0.0;
			for(k=0;k<d;k++){
			   dist[i*ny+j] = MAX(dist[i*ny+j],w[k]*fabs(x[i*d+k]-y[j*d+k]));
			}
		 }
	  }
   }
   else {
	  for(i=0;i<nx;i++){
		 for(j=0;j<ny;j++){
			dist[i*ny+j] = 0.0;
			for(k=0;k<d;k++){
			   dist[i*ny+j] = MAX(dist[i*ny+j],fabs(x[i*d+k]-y[j*d+k]));
			}
		 }
	  }
   }
}
