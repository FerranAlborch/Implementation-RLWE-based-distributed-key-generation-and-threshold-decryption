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
#include"functions.h"


/*NECESSARY LIBRARIES
*Openssl
*GMP
*MPFR
*FLINT
*/


/*TO COMPILE USE
*gcc encryption_sim.c functions.c -lm -lgmp -lssl -lcrypto -lflint -lmpfr -o encrypt.out -O2
*/

/*DONE********************************************
*Name: key_gen
*
*Description: generates private and public key of LPR cryptosystem and KH
*
*Arguments:   mpz_t aE[]: first part of public key
*             mpz_t bE[]: second part of public key
*             mpz_t s[]: secret key
*							fmpz_t KH[]: vector with all keys KH
*             mpf_t sigma: standard deviation of the discrete gaussian
*             int n: security parameter
*							int keys: number of keys
*             mpz_t q: modulo
***************************************************************/
void key_gen(fmpz_poly_t aE, fmpz_poly_t bE, fmpz_poly_t s, fmpz_t KH[],
	 						mpf_t sigma,int n,int keys, mpz_t q) {
	//Compute a random element in Zq for every key KH
	double rando;
	mpf_t r;
	mpf_t q2;
	mpf_t aux;
	mpz_t aux2;
	mpf_init(q2);
	mpf_set_z(q2,q);
	for(int i=0; i<keys; ++i) {
		mpf_init(r);
		mpf_init(aux);
		mpz_init(aux2);
		rando = rand_gen();
		mpf_set_d(r,rando);
		mpf_mul(aux,r,q2);
		round_mpf(aux2,aux);
		fmpz_set_mpz(KH[i],aux2);
		mpf_clear(r);
		mpf_clear(aux);
		mpz_clear(aux2);
	}
	mpf_clear(q2);

	// Choose s following the gaussian distribution
  disc_gauss_Rq(sigma,s,n,q);

  // Choose aE uniformly at random
  fmpz_poly_init(aE);
	rand_Rq(aE,n,q);

  //Compute bE
  fmpz_poly_t e;
  fmpz_poly_init(e);
  disc_gauss_Rq(sigma,e,n,q);
	fmpz_poly_t aux4;
	fmpz_poly_init(aux4);
	fmpz_poly_mul_Rq(aux4,aE,s,n,q);
  fmpz_poly_add(bE,aux4,e);
	fmpz_p_mod(bE,n,bE,q);
	fmpz_poly_clear(aux4);
	fmpz_poly_clear(e);
}

int main(int argc, char *argv[]) {
	//We set default float precission to 2048
	mpf_set_default_prec(2048);
	
	//Changeable parameters: n,q,u,t,repetitions
	int n = atoi(argv[1]);
	mpz_t q;
	mpz_init_set_str(q, "713623846352979940529142984724747568191373381", 10);
	int u= atoi(argv[2]);
	int t= atoi(argv[3]);
	int repetitions=10;

	//We set the general parameters sigmaenc,interdec and seed
	int keys = binomial(u,t);
	
	mpf_t auxf;
	mpf_init(auxf);
	mpf_t aux2f;
	mpf_init(aux2f);
	mpf_t qf;
	mpf_init(qf);
	mpf_set_z(qf,q);
	mpz_t auxz;
	mpz_init(auxz);
	// Compute kappa
	mpz_t kappa;
	mpz_init(kappa);
	mpz_ui_pow_ui(auxz,2,99);
	mpf_set_z(auxf,auxz);
	mpf_ui_div(auxf,n*u,auxf);
	mpf_mul(auxf,qf,auxf);
	mpf_add_ui(auxf,auxf,1);
	mpf_sqrt(auxf,auxf);
	mpf_sub_ui(auxf,auxf,1);
	mpf_div_ui(auxf,auxf,4*n*u);
	mpf_trunc(auxf,auxf);
	round_mpf(kappa,auxf);
	//gmp_printf("kappa = %Zd\n",kappa);
	
	// Compute xi
	mpf_t sigma;
	mpf_init(sigma);
	mpf_set_z(sigma,kappa);
	mpf_set_d(auxf,0.5);
	mpf_add(sigma,sigma,auxf);
	mpz_ui_pow_ui(auxz,2,120);
	mpf_set_z(auxf,auxz);
	mpf_div(sigma,sigma,auxf);
	mpf_set_d(auxf,sqrt(M_PI/2));
	mpf_mul(sigma,sigma,auxf);
	mpf_log(sigma,sigma);
	mpf_set_d(auxf,-2.0);
	mpf_mul(sigma,sigma,auxf);
	mpf_set_z(auxf,kappa);
	mpf_set_d(aux2f,0.5);
	mpf_add(auxf,auxf,aux2f);
	mpf_mul(auxf,auxf,auxf);
	mpf_div(sigma,auxf,sigma);
	mpf_sqrt(sigma,sigma);
	//gmp_printf("xi = %.Ff\n",sigma);	
	
	mpf_clear(auxf);
	mpf_clear(aux2f);
	mpf_clear(qf);
	mpz_clear(auxz);
	mpz_clear(kappa);
	
	
	// Generate the keys
	fmpz_poly_t s;
	fmpz_poly_t aE;
	fmpz_poly_t bE;
	fmpz_t KH[keys];

	fmpz_poly_init(s);
	fmpz_poly_init(aE);
	fmpz_poly_init(bE);
	for(int i=0; i<keys; ++i) fmpz_init(KH[i]);

	key_gen(aE,bE,s,KH,sigma,n,keys,q);
	
	
	// Simulate encryption
	fmpz_poly_t m1, uenc, v;
	fmpz_poly_init(m1);
	fmpz_poly_init(uenc);
	fmpz_poly_init(v);
	double avg = 0.0;
	for (int l=0; l<repetitions; ++l) {
		//Set all variables to zero again: mshare, times, m1, m2, uenc, v, e
		fmpz_poly_zero(m1);
		fmpz_poly_zero(uenc);
		fmpz_poly_zero(v);
		
		// Generate a random message m1
		for(int i=0; i<n; ++i) {
			double random = rand_gen();
			if(random<0.5) fmpz_poly_set_coeff_ui(m1,i,0);
			else fmpz_poly_set_coeff_ui(m1,i,1);
		}
		
		// Encrypt random message m1
		clock_t begin = clock();
		encrypt(m1,uenc,v,aE,bE,sigma,n,q,s);
		clock_t end = clock();
		
		// Compute time spent in player i
		double times = (double)(end-begin)/CLOCKS_PER_SEC;
		
		avg = avg + times;
	}
	
	avg = avg/repetitions*1000;
	//printf("Average encryption time: %fms\n", avg);
	FILE *f = fopen("encrypt.csv","a");
	fprintf(f,"%f\n",avg);
	fclose(f);
	
	mpz_clear(q);
	mpf_clear(sigma);
	fmpz_poly_clear(s);
	fmpz_poly_clear(aE);
	fmpz_poly_clear(bE);
	for(int i=0; i<keys; ++i) fmpz_clear(KH[i]);
	fmpz_poly_clear(m1);
	fmpz_poly_clear(uenc);
	fmpz_poly_clear(v);
}
