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
*gcc decryption_sim.c functions.c -lm -lgmp -lssl -lcrypto -lflint -lmpfr -o decrypt.out -O2
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

/*DONE******************************************************************
*Name: decrypt_sim
*
*Decription: Simulates the threshold decryption
*
*Arguments:
*************************************************************************/
void decrypt_sim(fmpz_poly_t mshare[], double times[], fmpz_poly_t sshare[],
					fmpz_poly_t uenc, fmpz_poly_t v, int **key_map,
					int **allowedt, int **matrixt, fmpz_t KH[], mpz_t inter,
					int n, int u, int t, mpz_t q) {
	// We set lambda as v+uenc
	fmpz_poly_t lambda;
	fmpz_poly_init(lambda);
	fmpz_poly_add(lambda,uenc,v);
	fmpz_p_mod(lambda,n,lambda,q);

	fmpz_poly_t PRSS;
	fmpz_poly_init(PRSS);
	for(int i =0; i<u; ++i) {
		// Begin timer
		clock_t begin = clock();
		
		// We compute the PRSS share of the player
		fmpz_poly_zero(PRSS);
		PRSS_share(PRSS,key_map,allowedt,matrixt,KH,lambda,inter,n,u,t,i+1,q);
		
		// We compute v-sshare[i]*u
		fmpz_poly_mul_Rq(mshare[i],uenc,sshare[i],n,q);
		fmpz_poly_sub(mshare[i],v,mshare[i]);
		fmpz_p_mod(mshare[i],n,mshare[i],q);
		
		// We add it to the PRSS share and finished
		fmpz_poly_add(mshare[i],mshare[i],PRSS);

		// End timer
		clock_t end = clock();

		// Compute time spent in player i
		times[i] = (double)(end-begin)/CLOCKS_PER_SEC;
	}
	
	fmpz_poly_clear(lambda);
	fmpz_poly_clear(PRSS);
}


int main(int argc, char *argv[]) {
	//We set default float precission to 2048
	mpf_set_default_prec(2048);
	
	//Changeable parameters: n,q,u,t,repetitions
	int n = atoi(argv[1]);
	int np=4;
	mpz_t q;
	mpz_t qp;
	mpz_init_set_str(qp, "17", 10);
	mpz_init_set_str(q, "713623846352979940529142984724747568191373381", 10);
	int u= atoi(argv[2]);
	int t= atoi(argv[3]);
	int repetitions=10;

	//We set the general parameters sigmaenc,interdec and seed
	int binomt = binomial(u,t);
	
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
	
	// Compute interdec
	mpz_t interdec;
	mpz_init(interdec);
	mpz_tdiv_q_ui(interdec,q,4);
	mpz_sub(interdec,interdec,kappa);
	mpz_mul(auxz,kappa,kappa);
	mpz_mul_si(auxz,auxz,2*n*u);
	mpz_sub(interdec,interdec,auxz);
	mpz_tdiv_q_ui(interdec,interdec,binomt);
	//gmp_printf("interdec = %Zd\n",interdec);
	
	mpf_clear(auxf);
	mpf_clear(aux2f);
	mpf_clear(qf);
	mpz_clear(auxz);
	
	//Set seed
	srand(time(NULL));

	
	//Compute the matrix of t+1 subsets of players
	int binomt1 = binomial(u,t+1);
	int** matrixt1 = (int**)malloc(u * sizeof(int*));
	for (int index=0;index<u;++index){
	matrixt1[index] = (int*)malloc(binomt1 * sizeof(int));
	for(int j=0; j<binomt1;++j) matrixt1[index][j]=0;
	}
	matrix(u,t+1,binomt1,matrixt1);
	
	//Compute the matrix of t subsets of players
	int** matrixt = (int**)malloc(u * sizeof(int*));
	for (int index=0;index<u;++index){
	matrixt[index] = (int*)malloc(binomt * sizeof(int));
	for(int j=0; j<binomt;++j) matrixt[index][j]=0;
	}
	matrix(u,t,binomt,matrixt);
	
	// Compute allowedt
	int binomut1 = binomial(u-1,u-t-1);
	int **allowedt = malloc(u*sizeof(int*));
	for(int i=0; i<u; ++i) {
		allowedt[i] = malloc(binomut1*sizeof(int));
		for(int j=0; j<binomut1; ++j) allowedt[i][j]=0;
	}
	for(int i=0; i<u; ++i) {
		int k=0;
		for(int j=0; j<binomt; ++j) {
			if(matrixt[i][j]==0) {
				allowedt[i][k]=j;
				k=k+1;
			}
		}
	}

	//Create the keys "map"
	int keys=binomt;
	int** key_map1 = (int**)malloc((u+keys) * sizeof(int*));
	for (int index=0;index<u+keys;++index){
		key_map1[index] = (int*)malloc(t * sizeof(int));
		for(int j=0; j<t;++j) key_map1[index][j]=0;
	}
	for(int j=0; j<t;++j) key_map1[0][j]=0;
	for(int j=0; j<t;++j) key_map1[1][j]=0;
	int* state = (int*)malloc(0 * sizeof(int));
	int* indexes = (int*)malloc(u * sizeof(int));
	for(int i=0; i<u;++i) indexes[i]=i+1;
	create_keys(state,indexes,t,key_map1,0,u,keys,t);
	free(state);
	free(indexes);

	//Compute the keys
	fmpz_poly_t s;
	fmpz_poly_t aE;
	fmpz_poly_t bE;
	fmpz_poly_t ekg;
	fmpz_t KH[keys];

	fmpz_poly_init(s);
	fmpz_poly_init(aE);
	fmpz_poly_init(bE);
	for(int i=0; i<keys; ++i) fmpz_init(KH[i]);

	key_gen(aE,bE,s,KH,sigma,n,keys,q);

	// Create Shamir shares of s
	fmpz_poly_t sshare[u];
	for(int i=0; i<u; ++i) fmpz_poly_init(sshare[i]);
	gen_shamir_Rq(sshare,s,n,u,t,q);



	// Simulation of decryption
	// Create storage for the shares of the message
	fmpz_poly_t mshare[u];
	for(int i=0; i<u; ++i) fmpz_poly_init(mshare[i]);
	
	//Create time vector
	double timesdec[u];
	for(int i=0; i<u; ++i) timesdec[i] = 0.0;

	// Encrypt and decrypt 200 random messages, verify correctness and print average
	// maximum time and average minimum time
	fmpz_poly_t m1, m2, uenc, v, e[binomt1];
	fmpz_poly_init(m1);
	fmpz_poly_init(m2);
	fmpz_poly_init(uenc);
	fmpz_poly_init(v);
	for(int i=0; i<binomt1; ++i) fmpz_poly_init(e[i]);
	double avgmin = 0.0;
	double avgmax = 0.0;
	int correctsim;
	for(int l=0; l<repetitions; ++l) {
		//Set all variables to zero again: mshare, times, m1, m2, uenc, v, e
		for(int i=0; i<u; ++i) fmpz_poly_zero(mshare[i]);
		fmpz_poly_zero(m1);
		fmpz_poly_zero(m2);
		fmpz_poly_zero(uenc);
		fmpz_poly_zero(v);
		for(int i=0; i<u; ++i) timesdec[i] = 0.0;
		for(int i=0; i<binomt1; ++i) fmpz_poly_zero(e[i]);

		// Generate a random message m1
		for(int i=0; i<n; ++i) {
			double random = rand_gen();
			if(random<0.5) fmpz_poly_set_coeff_ui(m1,i,0);
			else fmpz_poly_set_coeff_ui(m1,i,1);
		}

		// Encrypt random message m1
		encrypt(m1,uenc,v,aE,bE,sigma,n,q,s);

		// Decrypt random message
		decrypt_sim(mshare,timesdec,sshare,uenc,v,key_map1,allowedt,matrixt,KH,interdec,n,u,t,q);

		// Compute the minimum and maximum time and add it to the average
		double min = 10000000000.0;
		double max = 0.0;
		for(int i=0; i<u; ++i) {
			if(timesdec[i]<min) min=timesdec[i];
			if(timesdec[i]>max) max=timesdec[i];
		}
		avgmin=avgmin+min;
		avgmax=avgmax+max;

		// Verify that the threshold decryption is correct
		for(int i=0; i<binomt1; ++i) {
			// Find the vector of players
			int players[t+1];
			int count = 0;
			for(int k=0; k<u; ++k) {
				if(matrixt1[k][i]==1) {
					players[count]=k+1;
					count=count+1;
				}
			}
			compute_shamir_Rq(e[i],mshare,players,t,u,n,q);
			fmpz_p_mod(e[i],n,e[i],q);
			round_message(e[i],m2,q,n);
		 	correctsim = fmpz_poly_equal(m1,m2);
		}
		int compare=1;
		for(int i=0; i<binomt1; ++i) {
			int compaux = fmpz_poly_equal(e[0],e[i]);
			if (compaux==0) compare =0;
		}

		// Print whether the the threshold decryption is correct or not
		//printf("Returns 1 if all the threshold decryptions are equal, 0 otherwise: %d\n", compare);
		//printf("Returns 1 if the threshold decryption is correct, 0 otherwise: %d\n", correctsim);
	}

	// Compute the average times in miliseconds and print them
	avgmin = avgmin/repetitions*1000;
	avgmax = avgmax/repetitions*1000;
	FILE *f = fopen("decrypt.csv","a");
	fprintf(f,"%f\n",avgmax);
	fclose(f);
	//printf("Average minimum decryption time: %fms\n", avgmin);
	//printf("Average maximum decryption time: %fms\n", avgmax);


	// Free all storage
	mpz_clear(q);
	mpz_clear(qp);
	mpz_clear(interdec);
	mpf_clear(sigma);
	for(int index=0; index<u; ++index) free(matrixt[index]);
	for(int index=0; index<u; ++index) free(allowedt[index]);
	for(int index=0; index<u+keys; ++index) free(key_map1[index]);
	free(key_map1);
	fmpz_poly_clear(s);
	fmpz_poly_clear(aE);
	fmpz_poly_clear(bE);
	for(int i=0; i<keys; ++i) fmpz_clear(KH[i]);
	for(int i=0; i<u; ++i) fmpz_poly_clear(sshare[i]);
	for(int i=0; i<u; ++i) fmpz_poly_init(mshare[i]);
	fmpz_poly_clear(m1);
	fmpz_poly_clear(m2);
	fmpz_poly_clear(uenc);
	fmpz_poly_clear(v);
	for(int i=0; i<binomt1; ++i) fmpz_poly_clear(e[i]);
}
