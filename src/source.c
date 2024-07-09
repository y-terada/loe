#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

SEXP LOEobjt(SEXP DM, SEXP A, SEXP c)
	{
		int i, j, k, n;
		n = INTEGER(GET_DIM(A))[0];
		double stress=0, C=REAL( c )[0], *D=REAL(DM), *ADM=REAL(A);
		SEXP str;
		PROTECT(str=allocVector(REALSXP, 1));
		
		for(i = 0; i <n; i++){
			for(j = 0; j < n; j++){
				if(ADM[i +j*n]==1){
					for(k = 0; k < n; k++){
						if(i != k){
							if( ADM[i +k*n]==0 && D[i + j*n]+C>D[i + k*n]){
								stress += pow(D[i + n*j]+C-D[i + n*k],2);
							}
						}
					}
				}
			}
		}
		REAL(str)[0]=stress;
		UNPROTECT(1);
		return(str);
	}

void LOEgrad(double *grad, double *X,double *D, int *ADM, int *N, int *P, double *c)
	{
		int i, j, k, n, p, s;
		double C,tmpdij,tmpdik;
		n = *N;
		C = *c;
		p = *P;
		//
		for(i = 0; i <n; i++){// start i
			for(s=0; s<p;s++){//start for s
				for(j = 0; j < n; j++){//start for j
				/*~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~*/
					if(ADM[i +j*n]==1 ){//start if ADM[i +j*n]==1 
						for(k = 0; k < n; k++){//start for k
							if(ADM[i +k*n]==0){//start if ADM[i +k*n]==0
								if(D[i + j*n]+C>D[i + k*n]){
									if(D[i + j*n]<0.00001){
										tmpdij=0.00001;
									}else{
										tmpdij = D[i + j*n];
									}
									if(D[i + k*n]<0.00001){
										tmpdik=0.00001;
									}else{
										tmpdik=D[i + k*n];
									}
									grad[i+n*s] += 2*( D[i + j*n]+C-D[i + k*n] )*(  ( (X[i+n*s]-X[j+n*s] )/tmpdij ) - ( (X[i+n*s]-X[k+n*s] )/tmpdik)  );
									grad[k+n*s] += -2*(D[i + j*n]+C-D[i + k*n])*( (X[k+n*s]-X[i+n*s])/tmpdik  );
									grad[j+n*s]+= 2*(D[i + j*n]+C-D[i + k*n])*( (X[j+n*s]-X[i+n*s])/tmpdij  );
								}
							}// if ADM[i +k*n]==0
						}//end for k
					}//end if ADM[i +j*n]==1
				/*~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~*/
				}//end for j
			}//end for s
		}//end for i
	}
	
void fmajo(double *X, double *str, int *ADM, int *N, int *P, double *c, double *eps, int *maxit, int *report)
	{//start void
		int i, j, k, r, n, s, p, rp, MAXIT;
		double C, tmp=0, muij=0, muji=0, etaij=0, etaji=0, tmpd=0,tmpstr =0,tmpM=0, tmpH=0, EPS=0;
		double *alp,*alpa,*bet,*beta, *M, *H, *tmpX, *D;
		n = *N;
		C = *c;
		p = *P;
		MAXIT = *maxit+1;
		rp = *report;
		EPS = *eps;
		D=(double *)calloc(n*n,sizeof(double));
		M=(double *)calloc(n*n,sizeof(double));
		H=(double *)calloc(n*n,sizeof(double));
		tmpX=(double *)calloc(n*p,sizeof(double));
		/*start: copy data*/
		for(i=0;i<n;i++){
			for(s=0;s<p;s++){
				tmpX[i+s*n] = X[i+s*n];
			}
		}
		/*end: copy data*/
		/*-----------------------------------------------------*/
		/*start: compute distance matrix*/
		for(i=0; i<(n-1); i++){
			for(j =i+1; j <n; j++){
				tmpd =0;
				for(s=0; s<p; s++){
					tmpd+= pow(X[i+s*n]-X[j+s*n],2);
				}
				D[i+j*n] = sqrt(tmpd);
				D[j+i*n] = D[i+j*n]; 
			}
		}
		/*end: compute distance matrix*/
		/*start: compute initial stress*/
		tmpstr = 0;
		for(i = 0; i <n; i++){
			for(j = 0; j < n; j++){
				if(ADM[i +j*n]==1){
					for(k = 0; k < n; k++){
						if( ADM[i +k*n]==0 && D[i + j*n]+C>D[i + k*n]){
							tmpstr += pow(D[i + n*j]+C-D[i + n*k],2);
						}
					}
				}
			}
		}
		/*-----------------------------------------------------*/
		str[0] = tmpstr;
		/*end: compute  initial stress*/
		Rprintf("initial stress= %f \n", str[0]);
		/*=======================================*/
		for(r=1; r < MAXIT; r++){/*start for r*/
			/*start: copy data*/
			for(i=0;i<n;i++){
				for(s=0;s<p;s++){
					tmpX[i+s*n] = X[i+s*n];
				}
			}
			alp=(double *)calloc(n*n,sizeof(double));
			alpa=(double *)calloc(n*n,sizeof(double));
			bet=(double *)calloc(n*n,sizeof(double));
			beta=(double *)calloc(n*n,sizeof(double));
			/*end: copy data*/
			/*-----------------------------------------------------*/
			/*start compute coefficient*/
			tmp = 0.01; //tau>0, small positive instead of d_ij=0
			//Rprintf("%f \n", tmp);
			for(i=0; i <n; i++){
				for(j=0;j <n; j++){
					if(ADM[i+j*n]==1){
						if(D[i+j*n]>tmp){//if start D[i+j*n]>0
							for(k=0; k<n; k++){
								if(ADM[i+k*n]==0){
									alpa[i+k*n] += 2;
									if(D[i+j*n]+D[i+k*n]>=C){
										alp[i+j*n] += 2;
										if(D[i+j*n]+C>=D[i+k*n]){
											bet[i+j*n] += (D[i+j*n]+D[i+k*n]-C)/D[i+j*n];
										}
									}else{
										alp[i+j*n] += (D[i+j*n]+C-D[i+k*n])/D[i+j*n];
									}
									if(D[i+k*n]>0){
										if(D[i+j*n]+C>=D[i+k*n]){
											beta[i+k*n] += (D[i+j*n]+D[i+k*n]+C)/D[i+k*n];
										}
									}
									if(D[i+j*n]+C<D[i+k*n]){
										bet[i+j*n] += 2;
										if(D[i+k*n]>0){
											beta[i+k*n] += 2;
										}
									}
								}
							}
						}else if((D[i+j*n]>0) & (D[i+j*n]<tmp)){
							for(k=0; k<n; k++){
								if(ADM[i+k*n]==0){
									alpa[i+k*n] += 2;
									if(D[i+j*n]+D[i+k*n]>=C){
										alp[i+j*n] += 2;
										if(D[i+j*n]+C>=D[i+k*n]){
											bet[i+j*n] += (tmp+D[i+k*n]-C)/tmp;
										}
									}else{
										alp[i+j*n] += (tmp+C-D[i+k*n])/tmp;
									}
									if(D[i+k*n]>0){
										if(D[i+j*n]+C>=D[i+k*n]){
											beta[i+k*n] += (D[i+j*n]+D[i+k*n]+C)/D[i+k*n];
										}
									}
									if(D[i+j*n]+C<D[i+k*n]){
										bet[i+j*n] += 2;
										if(D[i+k*n]>0){
											beta[i+k*n] += 2;
										}
									}
								}
							}
						}else{//if else D[i+j*n]>0
							for(k =0;k<n;k++){
								if(ADM[i+k*n]==0){
									alpa[i+k*n] += 2;
									if(D[i+j*n]+D[i+k*n]>=C){
										alp[i+j*n] += 2;
									}else{
										alp[i+j*n] += (tmp+C-D[i+k*n])/tmp;
									}
									if(D[i+k*n]>0){
										if(D[i+j*n]+C>=D[i+k*n]){
											beta[i+k*n] += (tmp+D[i+k*n]+C)/D[i+k*n];
										}
									}
									if(D[i+j*n]+C<D[i+k*n]){
										bet[i+j*n] += 2;
										if(D[i+k*n]>0){
											beta[i+k*n] += 2;
										}
									}
								}
							}
						}//if end D[i+j*n]>0
					}//if end ADM[i+j*n]>0
				}//j
			}//i
			/*end compute coefficient*/
			/*-----------------------------------------------------*/
			/*-----------------------------------------------------*/
			/*start compute matrices M and H*/
			for(i=0;i <(n-1); i++){//for start i
				for(j =i+1; j <n; j++){//for start j
					muij = alp[i+j*n]+alpa[i+j*n];
					muji = alp[j+i*n]+alpa[j+i*n];
					M[i+j*n] = -2*muij;
					M[j+i*n] = -2*muji;
					etaij = bet[i+j*n]+beta[i+j*n];
					etaji = bet[j+i*n]+beta[j+i*n];
					H[i+j*n] = -(etaij + etaji) ;
					H[j+i*n] = -(etaij + etaji) ;
				}//for end j
			}//for end i
			free(alp);
			free(alpa);
			free(bet);
			free(beta);
			for(i=0;i <n; i++){//for start i
						tmpM=0;
						M[i+i*n] =0;
						H[i+i*n] =0;
						tmpH=0;
				for(j =0; j <n; j++){//for start j
						tmpM += (-M[i+j*n]/2-M[j+i*n]/2);
						tmpH += (-H[i+j*n]);
				}//for end j
						M[i+i*n] = tmpM;
						H[i+i*n] = tmpH;
			}//for end i
			/*end compute matrices M and H*/
			/*-----------------------------------------------------*/
			/*-----------------------------------------------------*/
			/*start renew*/
			for(i =0 ; i<n; i++){
				for(s=0; s<p;s++){
					tmp =0;
					for(j =0; j < n; j++){
						if(i ==j){
							tmp += 2*H[i+j*n]*tmpX[j+s*n];
						}else{
							tmp += 2*H[i+j*n]*tmpX[j+s*n] - (M[i+j*n]+M[j+i*n])*X[j+s*n];
						}
					}
					X[i+s*n] = tmp/(2*M[i+i*n]);
				}
			}
			/*end renew*/
			/*-----------------------------------------------------*/
			str[r]=0;
			for(i=0; i<(n-1); i++){
				for(j =i+1; j <n; j++){
					tmpd =0;
					for(s=0; s<p; s++){
						tmpd+= pow(X[i+s*n]-X[j+s*n],2);
					}
					D[i+j*n] = sqrt(tmpd);
					D[j+i*n] = D[i+j*n];
				}
			}
			/*start: compute stress*/
			tmpstr = 0;
			for(i = 0; i <n; i++){
				for(j = 0; j < n; j++){
					if(ADM[i +j*n]==1){
						for(k = 0; k < n; k++){
							if( ADM[i +k*n]==0 && D[i + j*n]+C>D[i + k*n]){
								tmpstr += pow(D[i + n*j]+C-D[i + n*k],2);
							}
						}
					}
				}
			}
			str[r] = tmpstr;
			if(r % rp==0){
				Rprintf("iter %d stress = %f \n",r, str[r]);
			}
			if(str[r-1]-str[r]<EPS){
				Rprintf("iter (final) %d stress = %f \n",r, str[r]);
				Rprintf("converged\n");
				break;
			}
		}/*end for r*/
		/*=======================================*/
		if(r==MAXIT){
			Rprintf("final (iter %d) stress = %f \n",r-1, str[r-1]);
			Rprintf("stopped after %d iterations \n",r-1);
		}
		free(tmpX);
		free(M);
		free(H);
		free(D);
	}//end void
	
void getorder(int *AM, double *D, int *N)
	{
		int cnt=0, i, j, k,l,n;
		n = *N;
		
		for(i = 0; i <(n-1); i++){
			for(j = i+1; j<n; j ++){
				for(k = 0; k <(n-1); k++){
					for(l = k+1; l<n; l ++){
						if(D[i+j*n]<D[k+l*n]){
							AM[4*cnt+0] = i+1;
							AM[4*cnt+1] = j+1;
							AM[4*cnt+2] = k+1;
							AM[4*cnt+3] = l+1;
							cnt += 1;
						}
					}
				}
			}
		}
	}
	
void SOEgrad(double *grad, double *X, double *D, int *AM, double *c, int *N, int *P, int *NC)
	{
		int a, i, j, k,l, nc,n,p,s;
		double C=*c, tmpdij, tmpdkl;
		nc = *NC;
		n = *N;
		p = *P;
		
		for(a = 0; a <nc; a++){
			i = AM[a]-1;
			j = AM[a+nc]-1;
			k = AM[a+2*nc]-1;
			l = AM[a+3*nc]-1;
			//to be safe
			if(D[i + j*n]<0.00001){
				tmpdij = 0.00001;
			}else{
				tmpdij = D[i + j*n];
			}
			if(D[k + l*n]<0.00001){
				tmpdkl = 0.00001;
			}else{
				tmpdkl = D[k + l*n];
			}
			//main computation
			if( D[i + j*n]+C>D[k + l*n]){
				if(i != k && i !=l){
					for(s = 0; s <p;s++){
						grad[i+n*s] += 2*((X[i+s*n]-X[j+s*n])/tmpdij)*(D[i + j*n]+C-D[k + l*n]);
						if(j!=k && j!=l){
							grad[j+n*s] += 2*((X[j+s*n]-X[i+s*n])/tmpdij)*(D[i + j*n]+C-D[k + l*n]);
							grad[k+n*s] += 2*((X[l+s*n]-X[k+s*n])/tmpdkl)*(D[i + j*n]+C-D[k + l*n]);
							grad[l+n*s] += 2*((X[k+s*n]-X[l+s*n])/tmpdkl)*(D[i + j*n]+C-D[k + l*n]);
						}else if(j==k){
							grad[j+n*s] += 2*( (X[j+s*n]-X[i+s*n])/tmpdij - (X[j+s*n]-X[l+s*n])/tmpdkl )*(D[i + j*n]+C-D[k + l*n]);
							grad[l+n*s] += 2*((X[k+s*n]-X[l+s*n])/tmpdkl)*(D[i + j*n]+C-D[k + l*n]);
						}else if(j==l){
							grad[j+n*s] += 2*( (X[j+s*n]-X[i+s*n])/tmpdij - (X[j+s*n]-X[k+s*n])/tmpdkl )*(D[i + j*n]+C-D[k + l*n]);
							grad[k+n*s] += 2*((X[l+s*n]-X[k+s*n])/tmpdkl)*(D[i + j*n]+C-D[k + l*n]);
						}
					}
				}else if(i==k){
					for(s = 0; s <p;s++){
						grad[i+n*s] += 2*( (X[i+s*n]-X[j+s*n])/tmpdij - (X[i+s*n]-X[l+s*n])/tmpdkl )*(D[i + j*n]+C-D[k + l*n]);
						grad[j+n*s] += 2*((X[j+s*n]-X[i+s*n])/tmpdij)*(D[i + j*n]+C-D[k + l*n]);
						grad[l+n*s] += 2*((X[k+s*n]-X[l+s*n])/tmpdkl)*(D[i + j*n]+C-D[k + l*n]);
					}
				}else if(i==l){
					for(s = 0; s <p;s++){
						grad[i+n*s] += 2*( (X[i+s*n]-X[j+s*n])/tmpdij - (X[i+s*n]-X[k+s*n])/tmpdkl )*(D[i + j*n]+C-D[k + l*n]);
						grad[j+n*s] += 2*((X[j+s*n]-X[i+s*n])/tmpdij)*(D[i + j*n]+C-D[k + l*n]);
						grad[k+n*s] += 2*((X[l+s*n]-X[k+s*n])/tmpdkl)*(D[i + j*n]+C-D[k + l*n]);
					}
				}
			}
		}
	}
	
void SOEobjt(double *D, int *AM, double *c, int *N, int *NC, double *str)
	{
		int a, i, j, k,l, nc,n;
		double C=*c, stress=0;
		nc = *NC;
		n = *N;
		
		for(a = 0; a <nc; a++){
			i = AM[a]-1;
			j = AM[a+nc]-1;
			k = AM[a+2*nc]-1;
			l = AM[a+3*nc]-1;
			if( D[i + j*n]+C>D[k + l*n]){
				stress = stress + pow(D[i + n*j]+C-D[k + n*l],2);
			}
		}
		str[0] = stress;
	}
	
	
void SDMLOE(double *X, double *str,int *ADM, int *N, int *P, double *c,int *maxiter, double *eps, double *del, double *h, int *report)
	{
		int i, j, k, r, n, p, s, M, one[]={1}, rp;
		double C,*tmpX,tmpd=0,tmpdij,tmpdik,tmpstr =0,*D,*grad,EPS,DEL,H,tmpH, gradnorm=0;
		n = *N;
		C = *c;
		p = *P;
		int np[] = {n*p};
		M = *maxiter+1;
		EPS = *eps;
		H = *h;
		DEL= *del;
		rp = *report;
		D=(double *)calloc(n*n,sizeof(double));
		grad=(double *)calloc(n*p,sizeof(double));
		tmpX = (double *)calloc(n*p,sizeof(double)); 
		/*start: copy data*/
		for(i=0;i<n;i++){
			for(s=0;s<p;s++){
				tmpX[i+s*n] = X[i+s*n];
			}
		}
		/*end: copy data*/
		/*start: compute distance matrix*/
		for(i=0; i<(n-1); i++){
			for(j =i+1; j <n; j++){
				tmpd =0;
				for(s=0; s<p; s++){
					tmpd+= pow(X[i+s*n]-X[j+s*n],2);
				}
				D[i+j*n] = sqrt(tmpd);
				D[j+i*n] = D[i+j*n]; 
			}
		}
		/*end: compute distance matrix*/
		/*start: compute initial stress*/
		tmpstr = 0;
		for(i = 0; i <n; i++){
			for(j = 0; j < n; j++){
				if(ADM[i +j*n]==1){
					for(k = 0; k < n; k++){
						if( ADM[i +k*n]==0 && D[i + j*n]+C>D[i + k*n]){
							tmpstr += pow(D[i + n*j]+C-D[i + n*k],2);
						}
					}
				}
			}
		}
		str[0] = tmpstr;
		/*end: compute  initial stress*/
		Rprintf("initial stress= %f \n", str[0]);
		for(r = 1; r<M; r++){// start i
			for(i=0;i<n;i++){
				for(s=0;s<p;s++){
					tmpX[i+s*n] = X[i+s*n];
					grad[i+s*n] = 0;
				}
			}
			/*start: compute grad*/
			for(i = 0; i <n; i++){// start i
				for(s=0; s<p;s++){//start for s
					for(j = 0; j < n; j++){//start for j
					/*~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~*/
						if(ADM[i +j*n]==1 ){//start if ADM[i +j*n]==1 
							for(k = 0; k < n; k++){//start for k
								if(ADM[i +k*n]==0){//start if ADM[i +k*n]==0
									if(D[i + j*n]+C>D[i + k*n]){
										if(D[i + j*n]<0.001){
											tmpdij=0.001;
										}else{
											tmpdij = D[i + j*n];
										}
										if(D[i + k*n]<0.001){
											tmpdik=0.001;
										}else{
											tmpdik=D[i + k*n];
										}//
										grad[i+n*s] += 2*( D[i + j*n]+C-D[i + k*n] )*(  ( (X[i+n*s]-X[j+n*s] )/tmpdij ) - ( (X[i+n*s]-X[k+n*s] )/tmpdik)  );
										grad[k+n*s] += -2*(D[i + j*n]+C-D[i + k*n])*( (X[k+n*s]-X[i+n*s])/tmpdik  );
										grad[j+n*s]+= 2*(D[i + j*n]+C-D[i + k*n])*( (X[j+n*s]-X[i+n*s])/tmpdij  );
									}
								}// if ADM[i +k*n]==0
							}//end for k
						}//end if ADM[i +j*n]==1
					/*~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~j~*/
					}//end for j
				}//end for s
			}/*end for i*/
			/*end: compute grad*/
			gradnorm = F77_CALL(ddot)(np, grad, one, grad, one);
			str[r] = str[r-1] +1;
			/*Line search (step method)*/
			tmpH = H;
			while(str[r-1]-DEL*tmpH*gradnorm<str[r]){
				for(i =0; i<n;i++){
					for(s=0;s<p;s++){
						//Rprintf("%f ", DEL*tmpH );
						X[i+s*n] = tmpX[i+s*n] - DEL*tmpH*grad[i+s*n];
					}
				}

				str[r]=0;
				for(i=0; i<(n-1); i++){
					for(j =i+1; j <n; j++){
						tmpd =0;
						for(s=0; s<p; s++){
							tmpd+= pow(X[i+s*n]-X[j+s*n],2);
						}
						D[i+j*n] = sqrt(tmpd);
						D[j+i*n] = D[i+j*n]; 
					}
				}
				/*start: compute stress*/
				tmpstr = 0;
				for(i = 0; i <n; i++){
					for(j = 0; j < n; j++){
						if(ADM[i +j*n]==1){
							for(k = 0; k < n; k++){
								if( ADM[i +k*n]==0 && D[i + j*n]+C>D[i + k*n]){
									tmpstr += pow(D[i + n*j]+C-D[i + n*k],2);
								}
							}
						}
					}
				}
				tmpH = tmpH*H;
				str[r] = tmpstr;
				/*end: compute stress*/
			}
			/*Line search (step method)*/
			if(r % rp==0){
				Rprintf("iter %d stress = %f \n",r, str[r]);
			}
			if(str[r-1]-str[r]<EPS){
				Rprintf("iter (final) %d stress = %f \n",r, str[r]);
				Rprintf("converged\n");
				break;
			}
		}/*end for r*/
		if(r==M){
			Rprintf("final (iter %d) stress = %f \n",r-1, str[r-1]);
			Rprintf("stopped after %d iterations \n",r-1);
		}
		free(tmpX);
		free(grad);
		free(D);
	}

	
