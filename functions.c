#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <openssl/evp.h>
#include <openssl/hmac.h>
#include <openssl/rand.h>
#include <openssl/sha.h>
#include <math.h>
#include <time.h>



/*DONE*************************************************************
 * Name: mpf_log
 * 
 * Description: Computes the natural logarithm of a given mpf_t
 * 
 * Arguments:	mpf_t result: result
 * 				mpf_t element: element to compute the logarithm
 * **************************************************************/
 void mpf_log(mpf_t result, mpf_t element) {
	 // Put element=d*2^exp
	 signed long int exp;
	 double d = mpf_get_d_2exp(&exp,element);
	 
	 // Compute log(element)=log(d)+exp*log(2)
	 mpf_set_d(result,log(d)+exp*log(2.0));
 }

/*DONE***********************************************************
*Name: fmpz_p_mod
*Description: returns the polynomial modulo q
*
*Arguments:   fmpz_poly_t a: polynomial
*             int l: number of coefficients in a
*             fmpz_poly_t b: result
*             mpz_t q: modulo
****************************************************************/
void fmpz_p_mod(fmpz_poly_t a, int l, fmpz_poly_t b, mpz_t q){
	fmpz_t qaux;
	fmpz_t q2;
	fmpz_t aux;
	fmpz_init(qaux);
	fmpz_init(q2);
	fmpz_set_mpz(qaux,q);
	fmpz_tdiv_q_ui(q2,qaux,2);
	for(int i = 0; i < l; i++){
		fmpz_init(aux);
		fmpz_poly_get_coeff_fmpz(aux,b,i);
		fmpz_mod(aux,aux,qaux);
		if(fmpz_cmp(q2,aux)<0) fmpz_sub(aux,aux,qaux);
		fmpz_poly_set_coeff_fmpz(a,i,aux);
		fmpz_clear(aux);
	}
	fmpz_clear(qaux);
	fmpz_clear(q2);
}

/*DONE****************************************************************
*Name: fmpz_poly_mul_Rq
*
*Description: Multiplies polynomials in Rq=Zq[x]/(x^n+1) c=a*b (mod Rq)
*
*Arguments: 	fmpz_poly_t c: result
*							fmpz_poly_t a: first polynomial
*							fmpz_poly_t b: second polynomial
*							int n: security parameter
*							mpz_t q: modulo
******************************************************************/
void fmpz_poly_mul_Rq(fmpz_poly_t c, fmpz_poly_t a, fmpz_poly_t b, int n,
											mpz_t q) {
	fmpz_poly_t aux;
	fmpz_poly_init(aux);
	fmpz_poly_mul_karatsuba(aux,a,b);
	fmpz_t aux2;
	fmpz_t aux3;
	for(int i=0; i<n; ++i) {
		fmpz_init(aux2);
		fmpz_init(aux3);
		fmpz_poly_get_coeff_fmpz(aux2,aux,i);
		fmpz_poly_get_coeff_fmpz(aux3,aux,i+n);
		fmpz_sub(aux2,aux2,aux3);
		fmpz_poly_set_coeff_fmpz(c,i,aux2);
		fmpz_clear(aux2);
		fmpz_clear(aux3);
	}
	fmpz_poly_clear(aux);
	

	fmpz_p_mod(c,n,c,q);
}

/*DONE**********************************************************
*Name: pop_back
*
*Description: takes an array, copies it without the last
*             element and deletes the previous one
*
*Arguments:   int vector[]: old vector
*             int new_vector[]: new vector
***********************************************************/
void pop_back(int vector[], int new_vector[], int length) {
  for(int i=0; i<length-1; ++i) new_vector[i]=vector[i];
  //free(vector);
  return;
}

/*DONE**********************************************************
*Name: push_back
*
*Decription:  takes an array and a value and copies it with
*             the value added at the end. Then deletes the
*             previous one.
*
*Arguments:   int vector[]:     old vector
*             int new_vector[]: new vector
*             int value:        value to add
***********************************************************/
void push_back(int vector[], int new_vector[], int value, int length) {
  for(int i=0; i<length; ++i) new_vector[i]=vector[i];
  new_vector[length]=value;
  //free(vector);
  return;
}

/*DONE**********************************************************
*Name: cmpfunc
*
*Description: Compare function for quick sort
*
*Arguments:   const void * a: Pointer to the first element
*                             to compare
*             const void * b: Pointer to the seconfd element
*                             to compare
***********************************************************/
int cmpfunc (const void * a, const void * b) {
   return ( *(int*)b - *(int*)a );
}

/*DONE**********************************************************
*Name: equal_vec
*
*Description: returns whether two vectors (of the same size)
*             are equal
*Arguments:   int* vec1:  first vector
*             int* vec2:  second vector
*             int size:
***********************************************************/
int equal_vec(int* vec1, int* vec2, int size) {
  for(int i=0; i<size; ++i) {
    if(vec1[i]!=vec2[i]) return 0;
  }
  return 1;
}

/*DONE**********************************************************
*Name: smaller_vec
*
*Description: returns whether the first vector is smaller
*             than the second
*
*Arguments:   int* vec1: first vector
*             int* vec2: second vector
*             int size:
***********************************************************/
int smaller_vec(int* vec1, int* vec2, int size) {
  for(int i=0; i<size; ++i) {
    if(vec1[i]<vec2[i]) return 1;
    if(vec1[i]>vec2[i]) return 0;
  }
  return 0;
}

/*DONE**********************************************************
*Name: find_keys
*
*Description: finds key position through binary search
*
*Arguments:   int** key_map1: map of the keys
*             int* key: key we want to find in decreasing order
*             int start: where the search start is
*             int end: where the search end is
***********************************************************/
int find_keys(int** key_map1, int* key, int size, int start, int end) {
  //First we set the ending condition
  if(start==end) return start;

  //Else we find the element in the middle
  int middle= (start+end)/2;
  int* vec_mid = (int*)malloc(size * sizeof(int));
  for(int i=0; i<size; ++i) vec_mid[i]=key_map1[middle][i];
  if (equal_vec(key, vec_mid, size)){
    free(vec_mid);
    return middle;
  }
  if (smaller_vec(key, vec_mid, size)) {
    free(vec_mid);
    return find_keys(key_map1, key, size, start, middle-1);
  }
  free(vec_mid);
  return find_keys(key_map1, key, size, middle+1, end);
}

/*DONE**********************************************************
*Name: create_keys
*
*Description: Function creating the key "map"
*
*Arguments:   int state[]: Vector of integers, signalling the "forbidden" players
*                          already chosen
*             int indexes[]: Vector of integers, signalling
*                            players allowed to be chosen
*             int size: Integer signalling how many
*                       players are left to pick
*             int** M: matrix
*             int length_state: Length of state vector
*             int length_indexes: Length of indexes vector
*             int keys: rows of M
*             int t: columns of M
***********************************************************/
void create_keys(int state[], int indexes[], int size, int** M,
              int length_state, int length_indexes, int keys, int t) {
  // First we give the condition to end the recursion
  if(size==0){
    // We need to sort state since the keys in the dictionary are sorted
    int row =0;
    for(int i=0; i<keys; ++i) {
      if(M[i][0]!=0) row=i+1;
    }
    for(int loop=0; loop<length_state;++loop) {
      M[row][loop]=state[loop];
    }
    return;
  }

  // Then the general case
  else{
    // We need to find all the "forbidden" subsets of players containing state,
    // so we compute all the "forbidden" subsets of players containing state and
    // the last element of indexes and all that do not contain it
    int a=indexes[length_indexes-1];
    int* indexes_new= malloc((length_indexes-1)*sizeof(int));
    pop_back(indexes, indexes_new, length_indexes);
    if(length_indexes>size) {
      create_keys(state,indexes_new,size,M,length_state,length_indexes-1,keys,t);
    }
    int* state_new=malloc((length_state+1)*sizeof(int));
    push_back(state, state_new, a, length_state);
    create_keys(state_new,indexes_new,size-1,M,length_state+1,length_indexes-1,keys,t);
    free(indexes_new);
    free(state_new);
  }
  return;
}

/*DONE**********************************************************
*Name: binom
*
*Description: Computes the binomial coefficient n over t
*
*Arguments:   int n:
*             int t:
***********************************************************/
int binomial(int n, int k) {
   if (k == 0 || k == n) return 1;
   return binomial(n - 1, k - 1) + binomial(n - 1, k);
}

/*DONE**********************************************************
*Name: subsets
*
*Description: Function aiding matrix
*
*Arguments:   int state[]:        Vector of integers, signalling
*                                 the "forbidden" players
*                                 already chosen
*             int indexes[]:      Vector of integers, signalling
*                                 players allowed to be chosen
*             int size:           Integer signalling how many
*                                 players are left to pick
*             int** M:            matrix
*             int length_state:   Length of state vector
*             int length_indexes: Length of indexes vector
*             int u:              rows of M
*             int binom:          columns of M
***********************************************************/
void subsets(int state[], int indexes[], int size, int** M,
              int length_state, int length_indexes, int u, int binom) {
  // First we give the condition to end the recursion
  if(size==0){
    int column =0;
    for(int i=0; i<u; ++i) {
      for(int j=0; j<binom; ++j) {
        if(M[i][j]==1) column=j+1;
      }
    }
    for(int loop=0; loop<length_state;++loop) {
      M[state[loop]][column]=1;
    }
    return;
  }

  // Then the general case
  else{
    // We need to find all the "forbidden" subsets of players containing state,
    // so we compute all the "forbidden" subsets of players containing state and
    // the last element of indexes and all that do not contain it
    int a=indexes[length_indexes-1];
    int* indexes_new= malloc((length_indexes-1)*sizeof(int));
    pop_back(indexes, indexes_new, length_indexes);
    if(length_indexes>size) {
      subsets(state,indexes_new,size,M,length_state,length_indexes-1,u,binom);
    }
    int* state_new=malloc((length_state+1)*sizeof(int));
    push_back(state, state_new, a, length_state);
    subsets(state_new,indexes_new,size-1,M,length_state+1,length_indexes-1,u,binom);
    free(indexes_new);
    free(state_new);
  }
  return;
}

/*DONE**********************************************************
*Name: matrix
*
*Description: creates the matrix of subsets of players
*
*Arguments:   int u:    total of players
              int t:    size of the subsets
              int** M:  matrix of subsets
***********************************************************/
void matrix(int u, int t, int binom, int** M) {
  int* state=malloc(0*sizeof(int));
  int* indexes=malloc(u*sizeof(int));
  for(int i=0; i<u; ++i) {
    indexes[i]=i;
  }
  subsets(state,indexes,t,M,0,u,u,binom);
  free(state);
  free(indexes);
}

/*DONE************************************************************
*Name: rand_gen
*
*Description: returns a random number uniformly distributed
*             between 0 and 1
*
*Arguments:
*************************************************************/
double rand_gen() {
	union {
        uint64_t i;
        unsigned char c[sizeof(uint64_t)];
    } u;

    if (!RAND_bytes(u.c, sizeof(u.c))) {
        fprintf(stderr, "Can't get random bytes!\n");
        exit(1);
    }
    /* 53 bits / 2**53 */
    return (u.i >> 11) * (1.0/9007199254740992.0);
}

/*DONE************************************************************
*Name: normal_sampler
*
*Description: returns a sample of a normal distribution with
*             mean mu and standard deviation sigma
*
*Arguments:   mpf_t sigma:   standard deviation
*							mpf_t sample:
*************************************************************/
void normal_sampler(mpf_t sigma, mpf_t sample) {
	 double u1=rand_gen();
   double u2=rand_gen();

	 //Compute sigma*cos(2*M_PI*u2)*sqrt(-2.*log(u1));
	 u1=sqrt(-2.*log(u1));
	 u2=cos(2*M_PI*u2);
	 mpf_t U1;
	 mpf_init(U1);
	 mpf_set_d(U1,u1);
	 mpf_t U2;
	 mpf_init(U2);
	 mpf_set_d(U2,u2);
	 mpf_mul(sample,U2,U1);
	 mpf_mul(sample,sample,sigma);
   mpf_clear(U1);
	 mpf_clear(U2);
}

/*DONE************************************************************
*Name: round_mpf
*
*Description: Rounds a multiple precission float to its nearest
*             multiple precission integer
*
*Arguments:   mpf_t x:    the value we want to round
*             mpz_t near: where we return the result
*************************************************************/
void round_mpf(mpz_t near, mpf_t x) {
  // First we round x up and down
  mpf_t f;
  mpf_init(f);
  mpf_floor(f,x);
  mpf_t c;
  mpf_init(c);
  mpf_ceil(c,x);

  // Then we compare the differences and take the smallest
  mpf_t subs1;
  mpf_init(subs1);
  mpf_t subs2;
  mpf_init(subs2);
  mpf_sub(subs1,x,f);
  mpf_sub(subs2,c,x);
  int comp = mpf_cmp(subs1,subs2);
  if(comp>0) {
    mpz_set_f(near,c);
    mpf_clear(f);
    mpf_clear(c);
    mpf_clear(subs1);
    mpf_clear(subs2);
    return;
  }
  if(comp==0) {
    mpz_set_f(near,c);
    mpf_clear(f);
    mpf_clear(c);
    mpf_clear(subs1);
    mpf_clear(subs2);
    return;
  }
  else {
    mpz_set_f(near,f);
    mpf_clear(f);
    mpf_clear(c);
    mpf_clear(subs1);
    mpf_clear(subs2);
    return;
  }
}

/*DONE************************************************************
*Name: disc_gauss_Zq
*
*Description: samples a discrete gaussian with standard
*             deviation sigma in Zq
*
*Arguments:   mpf_t sigma: standard deviation
*             mpz_t q: modulo
*             mpz_t sample: where we return the value
*************************************************************/
void disc_gauss_Zq(mpf_t sigma, mpz_t q, fmpz_t sample) {
  // To sample the discrete gaussian we first sample a normal
  // distribution of mean 0 and standard deviation sigma/2pi
	mpf_t x;
	mpf_init(x);
	mpf_t sigmaq;
	mpf_init(sigmaq);
	mpf_t qfloat;
	mpf_init(qfloat);
	mpf_set_z(qfloat,q);
	mpf_div(sigmaq,sigma,qfloat);
	normal_sampler(sigmaq,x);
	//printf("Normal sample: %f \n", x);
  // Then we take only the decimal part of it
  if(mpf_cmp_si(x,0)>=0) {
		mpf_t f;
		mpf_init(f);
		mpf_floor(f,x);
    mpf_sub(x,x,f);
		mpf_clear(f);
  }
  else {
		mpf_t c;
		mpf_init(c);
		mpf_ceil(c,x);
		mpf_sub(x,x,c);
		mpf_clear(c);
  }
	//printf("Normal tor sample: %f\n", x);

  // Now we multiply it by q. To do so we must convert q and
  // x to a multiple precission float
  mpf_t k;
  mpf_init(k);
  mpf_mul(k,x,qfloat);
	//gmp_printf("Tor sample mult q: %.*Ff \n", 6,k);
  mpf_clear(qfloat);
  mpf_clear(x);

  // Now we round it to the nearest integer in mpz_t
	mpz_t aux;
	mpz_init(aux);
  round_mpf(aux,k);
	fmpz_set_mpz(sample,aux);
	mpf_clear(k);
	mpz_clear(aux);
	//gmp_printf("Arrodonit: %Zd\n", sample);
}

/*DONE***********************************************************
*Name: disc_gauss_Rq
*
*Description: samples a discrete gaussian with standard
*             deviation sigma in Rq
*Arguments:   mpf_t sigma:   standard deviation
*             mpz_t* sample:  where we will store the sample
*             mpz_t q:        modulo
************************************************************/
void disc_gauss_Rq(mpf_t sigma, fmpz_poly_t sample, int n, mpz_t q) {
	fmpz_t x;
	for(int i =0; i<n; ++i) {
		fmpz_init(x);
		disc_gauss_Zq(sigma,q,x);
		fmpz_poly_set_coeff_fmpz(sample,i,x);
		fmpz_clear(x);
	}
	fmpz_p_mod(sample,n,sample,q);
	//mpz_p_repr(sample,n,sample,q);
}

/*DONE********************************************************************
*Name: rand_Rq
*
*Description: generates a random element in Rq
*
*Arguments:		fmpz_poly_t c: random element
*							int n: security parameter
*							mpz_t q: modulo
*************************************************************************/
void rand_Rq(fmpz_poly_t c, int n, mpz_t q) {
	double rando;
	mpf_t r;
	mpf_t q2;
	mpf_t aux;
	mpz_t aux2;
	fmpz_t aux3;
	mpf_init(q2);
	mpf_set_z(q2,q);
	for(int i=0; i<n; ++i) {
		mpf_init(r);
		mpf_init(aux);
		mpz_init(aux2);
		fmpz_init(aux3);
		rando = rand_gen();
		mpf_set_d(r,rando);
		mpf_mul(aux,r,q2);
		round_mpf(aux2,aux);
		fmpz_set_mpz(aux3,aux2);
		fmpz_poly_set_coeff_fmpz(c,i,aux3);
		mpf_clear(r);
		mpf_clear(aux);
		mpz_clear(aux2);
		fmpz_clear(aux3);
	}
	fmpz_p_mod(c,n,c,q);
	mpf_clear(q2);
}

/*DONE********************************************************************
*Name: rand_Zq
*
*Description: generates a random element in Zq
*
*Arguments:		fmpz_t c: random element
*							int n: security parameter
*							mpz_t q: modulo
*************************************************************************/
void rand_Zq(fmpz_t c, int n, mpz_t q) {
	double rando;
	mpf_t r;
	mpf_t q2;
	mpf_t aux;
	mpz_t aux2;
	fmpz_t qaux;
	fmpz_t qhalf;
	fmpz_init(qhalf);
	fmpz_init(qaux);
	fmpz_set_mpz(qaux,q);
	fmpz_tdiv_q_ui(qhalf,qaux,2);
	mpf_init(q2);
	mpf_set_z(q2,q);
    mpf_init(r);
    mpf_init(aux);
	mpz_init(aux2);
    rando = rand_gen();
    mpf_set_d(r,rando);
    mpf_mul(aux,r,q2);
  	round_mpf(aux2,aux);
	fmpz_set_mpz(c,aux2);
	if(fmpz_cmp(qhalf,c)<=0) fmpz_sub(c,c,qaux);
	mpf_clear(r);
	mpf_clear(q2);
	mpf_clear(aux);
	mpz_clear(aux2);
	fmpz_clear(qaux);
	fmpz_clear(qhalf);
}

/*DONE***************************************************************
 * Name: inverse_Zq
 * 
 * Definition: Computes the inverse of a number modulo an integer q
 * 			   with representative between -q/2 and q/2
 * 
 * Arguments: mpz_t inverse: where the result will be stored (must be
 * 			  initialized beforehand)
 * 			  mpz_t a: value for which we want to find the inverse
 * 			  mpz_t q: modulo
 * *****************************************************************/
 void inverse_Zq(mpz_t inverse, mpz_t a, mpz_t q) {
	 mpz_t t;
	 mpz_t newt;
	 mpz_t r;
	 mpz_t newr;
	 mpz_t aux;
	 mpz_t quotient;
	 mpz_init(t);
	 mpz_init(newt);
	 mpz_init(r);
	 mpz_init(newr);
	 mpz_init(aux);
	 mpz_init(quotient);
	 mpz_set_ui(t, 0);
	 mpz_set_ui(newt, 1);
	 mpz_set(r, q);
	 mpz_set(newr, a);
	 
	 while(mpz_cmp_ui(newr,0)!=0) {
		 mpz_tdiv_q(quotient,r,newr);
		 
		 mpz_set(aux,t);
		 mpz_set(t,newt);
		 mpz_mul(newt,newt,quotient);
		 mpz_sub(newt,aux,newt);
		 
		 mpz_set(aux,r);
		 mpz_set(r,newr);
		 mpz_mul(newr,newr,quotient);
		 mpz_sub(newr,aux,newr);
	 }
	 
	 if(mpz_cmp_ui(r,1)>0) {
		 // a is not invertible modulo q
		 gmp_printf("%Zd is NOT invertible modulo %Zd\n",a,q);
		 return;
	 }
	 
	 mpz_mod(inverse, t, q);
	 mpz_t q2;
	 mpz_init(q2);
	 mpz_tdiv_q_ui(q2,q,2);
	 if(mpz_cmp(inverse,q2)>0) mpz_sub(inverse,inverse,q);
	 
	 mpz_clear(t);
	 mpz_clear(newt);
	 mpz_clear(r);
	 mpz_clear(newr);
	 mpz_clear(aux);
	 mpz_clear(quotient);
	 mpz_clear(q2);
 }

/*DONE********************************************************************
*Name: gen_shamir_Rq
*
*Description:	Given an element in Rq returns u shares of this element with
*							threshold t
*
*Arguments:		fmpz_poly_t shares[]: array of shares (size t+1)
*							fmpz_poly_t element: element to share
*							int n: security parameter
*							int u: number of players
*							int t: threshold
**********************************************************************/
void gen_shamir_Rq(fmpz_poly_t shares[], fmpz_poly_t element, int n, int u, int t,
							mpz_t q) {
	// First we set all shares to the element
	for(int i=0; i<u; ++i) {
		fmpz_poly_zero(shares[i]);
		fmpz_poly_add(shares[i],shares[i],element);
	}

	// For every coefficient until t we generate a random polynomial in Rq
	// and for every player we multiply it by j^i to add it to its share
	fmpz_poly_t aux;
	fmpz_poly_t aux2;
	fmpz_t aux3;
	for(int i=1; i<t+1; ++i) {
		fmpz_poly_init(aux);
		rand_Rq(aux,n,q);
		// For every player j we add aux*j^i to the share
		for(int j=1; j<u+1; ++j) {
			fmpz_poly_init(aux2);
			fmpz_init(aux3);
			int d=pow(j,i);
			fmpz_set_si(aux3,d);
			fmpz_poly_scalar_mul_fmpz(aux2,aux,aux3);
			fmpz_poly_add(shares[j-1],shares[j-1],aux2);
			fmpz_poly_clear(aux2);
			fmpz_clear(aux3);
		}
		fmpz_poly_clear(aux);
	}
}

/*DONE********************************************************************
*Name: gen_shamir_Zq
*
*Description:	Given an element in Zq returns u shares of this element with
*							threshold t
*
*Arguments:		fmpz_t shares[]: array of shares (size t+1)
*							fmpz_t element: element to share
*							int n: security parameter
*							int u: number of players
*							int t: threshold
**********************************************************************/
void gen_shamir_Zq(fmpz_t shares[], fmpz_t element, int n, int u, int t,
							mpz_t q) {
	// First we set all shares to the element
	for(int i=0; i<u; ++i) fmpz_init_set(shares[i],element);

	// For every coefficient until t we generate a random element in Zq
	// and for every player we multiply it by j^i to add it to its share
	fmpz_t aux;
	fmpz_t aux2;
	fmpz_t aux3;
	for(int i=1; i<t+1; ++i) {
		fmpz_init(aux);
		rand_Zq(aux,n,q);
		// For every player j we add aux*j^i to the share
		for(int j=1; j<u+1; ++j) {
			fmpz_init(aux2);
			fmpz_init(aux3);
			int d=pow(j,i);
			fmpz_set_si(aux3,d);
			fmpz_mul(aux2,aux,aux3);
			fmpz_add(shares[j-1],shares[j-1],aux2);
			fmpz_clear(aux2);
			fmpz_clear(aux3);
		}
		fmpz_clear(aux);
	}
}

/*DONE********************************************************************
*Name: compute_shamir_Rq
*
*Description: Given exactly t+1 shares of the secret it reconstructs it
*
*Arguments:	fmpz_poly_t element: where the result will be stored
* 			fmpz_poly_t shares[]: shamir shares to reconstruct
* 			int players[]: vector of players participating in the 
* 						   reconstruction
* 			int t: threshold
* 			int u: number of players
* 			int n: dimension
* 			mpz_t q: modulo
****************************************************************************/
void compute_shamir_Rq(fmpz_poly_t element, fmpz_poly_t shares[], int players[],
												 int t, int u, int n, mpz_t q) {
	mpz_t inverse;
	mpz_init(inverse);
	mpz_t coeff;
	mpz_init(coeff);
	mpz_t aux_inv;
	mpz_init(aux_inv);
	fmpz_poly_t aux_share;
	fmpz_poly_init(aux_share);
	fmpz_poly_zero(element);
	
	
	// Adding for all players
	for(int i=0; i<=t; ++i) {
		// Compute the Lagrange coefficient modulo q
		mpz_set_si(coeff,1);
		for(int j=0; j<=t; ++j) {
			if(j!=i) {
				mpz_set_si(aux_inv,players[j]-players[i]);
				if(mpz_cmp_ui(aux_inv,0)<0) mpz_add(aux_inv,aux_inv,q);
				mpz_set_ui(inverse,0);
				inverse_Zq(inverse,aux_inv,q);
				mpz_mul(coeff,coeff,inverse);
				mpz_mul_si(coeff,coeff,players[j]);
				mpz_mod(coeff,coeff,q);
			}
		}
		
		// Multiply it to the share and add it to element
		fmpz_poly_scalar_mul_mpz(aux_share,shares[players[i]-1],coeff);
		fmpz_poly_add(element,element,aux_share);
		
		// Reduce modulo q
		fmpz_p_mod(element,n,element,q);
	}
	
	mpz_clear(inverse);
	mpz_clear(coeff);
	mpz_clear(aux_inv);
	fmpz_poly_clear(aux_share);
}

/*DONE********************************************************************
*Name: compute_shamir_Zq
*
*Description: Given exactly t+1 shares of the secret it reconstructs it
*
*Arguments:
****************************************************************************/
void compute_shamir_Zq(fmpz_t element, fmpz_t shares[], int players[],
												 int t, int u, int n, mpz_t q) {
	mpz_t inverse;
	mpz_init(inverse);
	mpz_t coeff;
	mpz_init(coeff);
	mpz_t aux_inv;
	mpz_init(aux_inv);
	fmpz_t aux_share;
	fmpz_init(aux_share);
	fmpz_zero(element);
	fmpz_t qfmpz;
	fmpz_init(qfmpz);
	fmpz_set_mpz(qfmpz,q);
	fmpz_t aux_coeff;
	fmpz_init(aux_coeff);
	
	
	// Adding for all players
	for(int i=0; i<=t; ++i) {
		// Compute the Lagrange coefficient modulo q
		mpz_set_si(coeff,1);
		for(int j=0; j<=t; ++j) {
			if(j!=i) {
				mpz_set_si(aux_inv,players[j]-players[i]);
				if(mpz_cmp_ui(aux_inv,0)<0) mpz_add(aux_inv,aux_inv,q);
				mpz_set_ui(inverse,0);
				inverse_Zq(inverse,aux_inv,q);
				mpz_mul(coeff,coeff,inverse);
				mpz_mul_si(coeff,coeff,players[j]);
				mpz_mod(coeff,coeff,q);
			}
		}
		fmpz_set_mpz(aux_coeff,coeff);
		
		// Multiply it to the share and add it to element
		fmpz_mul(aux_share,shares[players[i]-1],aux_coeff);
		fmpz_add(element,element,aux_share);
		
		// Reduce modulo q
		fmpz_mod(element,element,qfmpz);
	}
	
	mpz_clear(inverse);
	mpz_clear(coeff);
	mpz_clear(aux_inv);
	fmpz_clear(aux_share);
	fmpz_clear(qfmpz);
	fmpz_clear(aux_coeff);
}

/*DONE********************************************************
 * Name: compute_polynomials
 * 
 * Definition: Computes evaluations for the polynomials in the PRSS
 * 
 * Arguments: mpz_t eval: where the result is stored
 * 			  int *players: vectors of players relating to the polynomial
 * 			  int z: place where the polynomial will be evaluated
 * 			  int t: threshold
 * 			  mpz_t q: modulo
 * ******************************************************/
void compute_polynomials(mpz_t eval, int *players, int z, int t, mpz_t q) {
	// We take the value x and compute its evaluation as product in Zq
	mpz_t inverse;
	mpz_t a;
	mpz_init(inverse);
	mpz_init(a);
	mpz_set_si(eval,1);
	for (int i=0; i<t; ++i) {
		int aux = players[i] - z;
		mpz_set_si(a,players[i]);
		inverse_Zq(inverse,a,q);
		mpz_mul_si(eval,eval,aux);
		mpz_mul(eval,eval,inverse);
		mpz_mod(eval,eval,q);
	}
	mpz_clear(inverse);
	mpz_clear(a);
}

/*DONE******************************************************************
*Name: mx_hmac_sha3_512
*
*Description: Computes HMAC-Sha3_512 for a given key and given data
*
*Arguments:   const void *key: key
*             int keylen: key length
*             const unsigned char *data: data
*             int datalen: data length
*             unsigned char *result: storing the HMAC
*             unsigned int *resultlen: result length
**************************************************************************/
unsigned char *mx_hmac_sha3_512(const void *key, int keylen,
                              const unsigned char *data, int datalen,
                              unsigned char *result, unsigned int *resultlen) {
    return HMAC(EVP_sha3_512(), key, keylen, data, datalen, result, resultlen);
}

/*DONE**********************************************************************
*Name: PRF
*
*Description: Computes pseudo-random polynomial (as in a pseudo-random funciton)
*             in the interval [-inter, inter] (inter<q) with key k
*
*Arguments:   fmpz_t key: element in Zq that acts as the seed
*             mpz_poly_t lambda: argument of the PRF
*             fmpz_poly_t result: evaluation of the PRF
*             mpz_t inter: limit of the interval of the range (inter<q/4)
*             int n: security parameter (number of coefficients in Rq)
**************************************************************************/
void PRF(fmpz_t key, fmpz_poly_t lambda, fmpz_poly_t result, mpz_t inter, int n) {
  char *keystring=fmpz_get_str(NULL,2,key);
  int keylen=strlen(keystring);
  /*printf("Keystring:\n");
  for (unsigned int i = 0; i < keylen; i++){
    printf("%c", keystring[i]);
  }
  printf("\n");
  printf("\n");*/

  // We compute for every coefficient
  fmpz_t coeff;
  fmpz_t aux;
  mpz_t mod;
  mpz_t inter2;
  mpz_t two;
  mpz_init(inter2);
  mpz_init(two);
  mpz_set_si(two,2);
  mpz_mul(inter2,two,inter);
  for (int i=0; i<n; ++i) {
    fmpz_init(coeff);
    fmpz_init(aux);
    fmpz_poly_get_coeff_fmpz(coeff,lambda,i);
    /*printf("Coeff ");
    fmpz_print(coeff);
    printf("\n");*/
    unsigned char *coeffstring=fmpz_get_str(NULL,2,coeff);
    int coefflen=strlen((char*)coeffstring);
    /*printf("Coeffstr:\n");
    for (unsigned int i = 0; i < coefflen; i++){
      printf("%c", coeffstring[i]);
    }
    printf("\n");*/

    //Compute the HMAC of the coefficient using the key
    unsigned char *hash = NULL;
    unsigned int hashlen = -1;
    hash = mx_hmac_sha3_512((const void *)keystring, keylen, coeffstring, coefflen, hash, &hashlen);
    
    flint_free(coeffstring);

    char buffer[2*hashlen+1];
    char *ptr= &buffer[0];
    for (unsigned int i = 0; i < hashlen; i++){
      ptr +=sprintf(ptr,"%02X",(unsigned char) hash[i]);
    }

    //Transform the word in hexadecimal into an arbitrarily big integer
    int verify=fmpz_set_str(aux,buffer,16);

    //Reduce modulo 2inter and then change it into [-inter,inter]
    mpz_init(mod);
    fmpz_get_mpz(mod,aux);
    mpz_mod(mod,mod,inter2);
    if(mpz_cmp(mod,inter)>=0) mpz_sub(mod,mod,inter2);

    //Set mod as the ith coefficient of result
    //fmpz_init(aux);
    fmpz_set_mpz(aux,mod);
    fmpz_poly_set_coeff_fmpz(result,i,aux);
    fmpz_clear(coeff);
	fmpz_clear(aux);
	mpz_clear(mod);
	//free(ptr);
  }
  mpz_clear(inter2);
  mpz_clear(two);
}

/*DONE***********************************************************
 * Name: PRSS_share
 * 
 * Definition: computes the PRSS share of a player
 * 
 * Arguments:	fmpz_poly_t result: where the PRSS share will be stored
 * 				int **key_map: map storing the keys
 * 				int **allowedt: matrix of allowed t-subsets per player
 * 				int **matrixt: matrix of allowed t-subsets in general
 * 				fmpz_t KH[]: PRSS keys
 * 				fmpz_poly_t lambda: argument for the PRF
 * 				mpz_t inter: image interval for the PRF
 * 				int n: dimension
 * 				int u: number of players
 * 				int t: threshold
 * 				int player: number of the player computing the share
 * 				mpz_t q: modulo
 * **********************************************************/
 void PRSS_share(fmpz_poly_t result, int **key_map, int **allowedt, 
			     int **matrixt, fmpz_t KH[], fmpz_poly_t lambda, 
			     mpz_t inter, int n, int u, int t, int player, mpz_t q) {
	 // Set result to zero
	 fmpz_poly_zero(result);
	 
	 // For each allowed subgroup of t players containing player add the
	 // contribution to the result
	 int binomut1 = binomial(u-1,u-t-1);
	 int *players = malloc(t*sizeof(int));
	 mpz_t eval;
	 mpz_init(eval);
	 fmpz_poly_t aux;
	 fmpz_poly_init(aux);
	 for(int i=0; i<binomut1; ++i) {
		 fmpz_poly_zero(aux);
		 
		// We get the vector of players
		int k=t-1;
		for(int j=0; j<u; ++j) {
			if(matrixt[j][allowedt[player-1][i]]==1){
				players[k] = j+1;
				k=k-1;
			}
		}
		
		// We compute the evaluation of the polynomial
		compute_polynomials(eval, players, player, t, q);
		
		// We compute the PRF value corresponding to the subset of players
		int binomt = binomial(u,t);
		int r = find_keys(key_map,players,t,0,binomt-1);
		PRF(KH[r],lambda,aux,inter,n);
		
		// We compute the scalar product of the PRF value with eval and 
		// add it to result
		fmpz_poly_scalar_mul_mpz(aux,aux,eval);
		fmpz_poly_add(result,result,aux);
	 }
	 
	 fmpz_p_mod(result, n, result, q);
	 
	 mpz_clear(eval);
	 fmpz_poly_clear(aux);
	 free(players);
 }

/*DONE*****************************************************************
*Name: commit_sha2
*
*Descrpition: computes commitment of a message as a sha512 hash
*
*Arguments: 	char *message: message to commit
*							int message_len: message length
*							unsigned char hash[]: where the raw hash will be stored
*																		(size SHA512_DIGEST_LENGTH)
*							char shash[]: where the hexadecimal hash string will be stored
*														(size 2*SHA512_DIGEST_LENGTH+1)
************************************************************************/
void commit_sha2(char *message, int message_len, unsigned char hash[], char shash[]) {
	SHA512_CTX ctx;

	SHA512_Init(&ctx);
	SHA512_Update(&ctx, message, message_len);
	SHA512_Final(hash, &ctx);

	char *ptr= &shash[0];
	for (unsigned int i = 0; i < SHA512_DIGEST_LENGTH; i++){
		ptr +=sprintf(ptr,"%02X",(unsigned char) hash[i]);
	}
}

/*DONE****************************************************************
*Name: verify_commit_sha2
*
*Description: returns 0 if the commitment is the hexadecimal hash sha2 string of
*							message and 1 otherwise
*
*Arguments:		char commitment[]: hexadecimal hash string (size2*SHA512_DIGEST_LENGTH+1)
*							char *message: message to verify
*							int message_len: message length
************************************************************************/
int verify_commit_sha2(char commitment[], char *message, int message_len){
	//Compute hash of message
	unsigned char hash[SHA512_DIGEST_LENGTH];
	char shash[2*SHA512_DIGEST_LENGTH+1];
	commit_sha2(message,message_len,hash,shash);

	if(strcmp(commitment,shash)==0) return 1;
	return 0;
}

/*DONE**************************************************************
*Name: encrypt
*
*Description: Encrypts a given message (element in R_q of all 0 and 1)
*
*Arguments:		fmpz_poly_t m: message as polynomial of 0 and 1 coeffs
*							fmpz_poly_t u: first part of the ciphertext (initialized at 0)
*							fmpz_poly_t v: second part of the ciphertext (initialized at 0)
*							fmpz_poly_t aE: first part of the public key
*							fmpz_poly_t bE: second part pf the public key
*							mpf_t sigma: standard deviation used to encrypt
*							int n: size of the elements in Rq
*							mpz_t q: modulo
****************************************************************/
void encrypt(fmpz_poly_t m, fmpz_poly_t u, fmpz_poly_t v, fmpz_poly_t aE,
							fmpz_poly_t bE, mpf_t sigma, int n, mpz_t q,fmpz_poly_t s) {
	//rE,eu,ev following the Gaussian distribution
	fmpz_poly_t rE, eu, ev;
	fmpz_poly_init(rE);
	fmpz_poly_init(eu);
	fmpz_poly_init(ev);
	disc_gauss_Rq(sigma,rE,n,q);
	disc_gauss_Rq(sigma,eu,n,q);
	disc_gauss_Rq(sigma,ev,n,q);

	//u=aE*rE+eu
	fmpz_poly_mul_Rq(u,aE,rE,n,q);
	fmpz_poly_add(u,u,eu);
	fmpz_p_mod(u,n,u,q);

	//v=bE*rE+ev+m*q/2
	mpz_t q2;
	mpz_init(q2);
	mpz_tdiv_q_ui(q2,q,2);
	fmpz_poly_scalar_mul_mpz(v,m,q2);
	fmpz_poly_mul_Rq(rE,bE,rE,n,q);
	fmpz_poly_add(v,v,rE);
	fmpz_poly_add(v,v,ev);
	fmpz_p_mod(v,n,v,q);

	//clear extra polys and mpz
	fmpz_poly_clear(rE);
	fmpz_poly_clear(eu);
	fmpz_poly_clear(ev);
	mpz_clear(q2);
}

/*DONE***************************************************************
*Name: round_message
*
*Description: Given an element in Rq returns an element in Rq rounding to 0 or
*							q/2 and then assigning q/2 to 1
*
*Arguments:	fmpz_poly_t e: element to round
* 			fmpz_poly_t m: result
* 			mpz_t q: modulo
* 			int n: dimension
*******************************************************************/
void round_message(fmpz_poly_t e, fmpz_poly_t m, mpz_t q, int n)  {
	//For every coefficient we make the rounding
	fmpz_t fq,q2,q4,nq4;
	fmpz_init(fq);
	fmpz_init(q2);
	fmpz_init(q4);
	fmpz_init(nq4);
	fmpz_set_mpz(fq,q);
	fmpz_tdiv_q_ui(q2,fq,2);
	fmpz_tdiv_q_ui(q4,fq,4);
	fmpz_mul_si(nq4,q4,-1);
	fmpz_t coeff;
	for(int i=0; i<n; ++i) {
		fmpz_init(coeff);
		fmpz_poly_get_coeff_fmpz(coeff,e,i);
		if(fmpz_cmp_ui(coeff,0)>=0) {
			if(fmpz_cmp(coeff,q4)<=0) fmpz_poly_set_coeff_ui(m,i,0);
			else fmpz_poly_set_coeff_ui(m,i,1);
		}
		else {
			if(fmpz_cmp(coeff,nq4)>=0) fmpz_poly_set_coeff_ui(m,i,0);
			else fmpz_poly_set_coeff_ui(m,i,1);
		}
		fmpz_clear(coeff);
	}
	fmpz_clear(fq);
	fmpz_clear(q2);
	fmpz_clear(q4);
	fmpz_clear(nq4);
}
