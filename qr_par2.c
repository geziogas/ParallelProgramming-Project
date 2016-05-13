#include <omp.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int timer(void);

double *A,*Q,*R;
double *c,*t;
int n;

int main(int argc, char *argv[]) {
	int i,j,k,m,thrid,num_threads;
	int qr_eval=0,eval=0; //Debuging FLAGS
	clock_t start, diff;
	double error,d,time,s_time;
	omp_lock_t *lock;
	

	//Error message in entering arguments, check
	if (argc!=3){
		printf("Use: ./par [size] [num_threads]\n");
		return 0;
	}
	n = atoi(argv[1]);
	num_threads = atoi(argv[2]);
	omp_set_num_threads(num_threads);
	//printf("Test: Number of threads: %d\n", omp_get_num_threads());
  
   	//Allocate the matrices   
   	A = (double *)malloc(n*n*sizeof(double));
   	t = (double *)malloc(n*n*sizeof(double));
   	Q = (double *)malloc(n*n*sizeof(double));
   	R = (double *)malloc(n*n*sizeof(double));	   	
   
   	//Define an array of n locks
   	//and initialize them.
	lock = (omp_lock_t *)malloc(n*sizeof(omp_lock_t));
	for(i=0;i<n;i++) 
		omp_init_lock(&lock[i]);


	time=omp_get_wtime();	//Start counting exec time
	double r_sum = 0;
	#pragma omp parallel private(i,j,k,r_sum,thrid)
	{
		thrid=omp_get_thread_num();

	   	//Matrix A (initial), matrix 't' a temporary dublicate matrix
	   	//to keep the initial values of A    
	  	#pragma omp for schedule(static,1)
	  	for (i = 0; i<n; i++)
	  	{
		  	for(j=0;j<n;j++)
		  	{
			  	A[i*n+j] = (rand() % 40);
			  	t[i*n+j] = A[i*n+j];			  	
		  	}
		  	omp_set_lock(&lock[i]);		  	
		}
	
		//Parallel version of QR-factorization using 
		//modified Gram-Schmidt algorithm with locks
	
		// First column of ( Q[][0] )
		//Thread 0 calculates the 1st column
		//and unsets the 1st lock.
		r_sum=0;
		if(thrid==0)
		{
			// Calculation of ||A||
			for (i=0; i<n; i++){			
				r_sum = r_sum + A[0*n+i] * A[0*n+i]; 
			}
			R[0*n+0] = sqrt(r_sum);  			
			for (i=0; i<n; i++) {
				Q[0*n+i] = A[0*n+i]/R[0*n+0];							
	      	}
	      	omp_unset_lock(&lock[0]);
		}

	    for (k=1; k<n; k++)
	    {

	    	//Check if Q[][i-1] (the previous column) is computed.   	
	    	omp_set_lock(&lock[k-1]);
	    	omp_unset_lock(&lock[k-1]);
			
	  		#pragma omp for schedule(static,1) nowait
		    for(j=0; j<n; j++)
		    {	
		    	if(j>=k)
		    	{	    	
			        R[(k-1)*n+j]=0;	
			        for(i=0; i<n; i++) {			        	
			        	R[j*n+(k-1)] += Q[(k-1)*n+i] * A[j*n+i];			        	
			        } 
			        for (i=0; i<n; i++) {				        	
			        	A[j*n+i] = A[j*n+i] - R[j*n+(k-1)]*Q[(k-1)*n+i];
			        }			        
			       
			       	//Only one thread calculates the norm ||A||
			       	//and unsets the lock for the next column.
					if(j==k)
					{
						thrid=omp_get_thread_num();
						r_sum=0;
						for (i=0; i<n; i++){			
							r_sum = r_sum + A[k*n+i] * A[k*n+i]; 
						}
						R[k*n+k] = sqrt(r_sum);  
						
						//#pragma omp for schedule(static,1) nowait
						for (i=0; i<n; i++) {
							Q[k*n+i] = A[k*n+i]/R[k*n+k];			
				      	}
				      	//printf("I am thread %d and I unset lock %d.\n", thrid,k);
				      	omp_unset_lock(&lock[k]);
			      	}
		      	}		        
			}			
			r_sum=0;	      	
		}
	}
	time=omp_get_wtime()-time;
   	printf("(Parallel execution)Elapsed time: %f \n",time);
   	
   	//*****************************************************
   	//Code for evaluating the results for small matrices.
   	// (after changing the FLAGS)
   	c = (double *)malloc(n*n*sizeof(double));   	

   	if(eval==1) {
   		printf("\nPrinting A...\n");
		for (i = 0; i<n; i++){
		  	for(j=0;j<n;j++){			  	
			  	printf("%.3f ", t[j*n+i]);
		  	}
		  	printf("\n");
		}

		printf("\nPrinting Q...\n");
		for (i = 0; i<n; i++){
		  	for(j=0;j<n;j++){			  	
			  	printf("%.3f ", Q[j*n+i]);
		  	}
		  	printf("\n");
		}

		printf("\nPrinting R...\n");
		for (i = 0; i<n; i++){
		  	for(j=0;j<n;j++){			  	
			  	printf("%.3f ", R[j*n+i]);
		  	}
		  	printf("\n");
		}			
	}	
	
	//Evaluating the decomposition
	double sum; 
	if(qr_eval==1){
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				sum = 0;
				for (k = 0; k < n; k++) {
					sum = sum + Q[k*n+i] * R[j*n+k];
				}
				c[j*n+i] = sum;
			}
		}
		if(eval==1){
		printf("\nQ*R (Init A matrix) : \n");
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					printf("%.3f ", c[j*n+i]);
				}
				printf("\n");
			}
		}
		error=0;
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {				
				error+=fabs(c[j*n+i]-t[j*n+i]);
			}			
		}
		printf("\nError: %e\n",error);
	}	
	free(A);free(Q);free(R);free(t);free(c);free(lock);
	return 0;
}