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
*gcc keygen_sim.c functions.c -lm -lgmp -lssl -lcrypto -lflint -lmpfr -o keygen.out -O2
*/



/*DONE***************************************************************
*Name: decrypt
*
*Description: Returns m given u,v,s
*
*Arguments:	fmpz_poly_t m: result message
* 			fmpz_poly_t s: secret key
* 			fmpz_poly_t u: part of the encryption
* 			fmpz_poly_t v: part of the encryption
* 			mpz_t q: modulo
* 			int n: dimension
********************************************************************/
void decrypt(fmpz_poly_t m, fmpz_poly_t s, fmpz_poly_t u, fmpz_poly_t v,
				mpz_t q, int n) {
	fmpz_poly_t e;
	fmpz_poly_init(e);
	fmpz_poly_mul_Rq(e,s,u,n,q);
	fmpz_poly_sub(e,v,e);
	fmpz_p_mod(e,n,e,q);
	round_message(e,m,q,n);
	fmpz_poly_clear(e);
}

/*DONE******************************************************************
*Name: keygen_sim
*
*Decription: Simulates the distributed key generation
*
*Arguments:
*************************************************************************/
void keygen_sim(double **times, fmpz_poly_t sshare[],
				fmpz_poly_t aE, fmpz_poly_t bE, fmpz_t KH[], int** allowedt, fmpz_poly_t lambda,
				int** matrixt1, int **matrixt, int **key_map1, int binomt, mpf_t sigmakg,
				int binomt1, int u, int n, int t, mpz_t q, mpz_t interkg) {
	// Generate storage: sNIVSS_keys, eNIVSS_keys, si, commit_si, ei, 
	// commit_ei, aux_KH, KHj_shares, commit_KHj_shares, aux_aE,
	// aEj_shares, commit_aEj_shares, auxNIVSS
	fmpz_t **sNIVSS_keys;
	sNIVSS_keys = flint_malloc(u*sizeof(fmpz_t*));
	for(int i=0; i<u; ++i) {
		sNIVSS_keys[i] = flint_malloc(binomt*sizeof(fmpz_t));
		for(int j=0; j<binomt; ++j) fmpz_init(sNIVSS_keys[i][j]);
	}
	fmpz_t **eNIVSS_keys;
	eNIVSS_keys = flint_malloc(u*sizeof(fmpz_t*));
	for(int i=0; i<u; ++i) {
		eNIVSS_keys[i] = flint_malloc(binomt*sizeof(fmpz_t));
		for(int j=0; j<binomt; ++j) fmpz_init(eNIVSS_keys[i][j]);
	}
	char **commit_si;
	commit_si = malloc(u*sizeof(char*));
	for(int i=0; i<u; ++i) {
		commit_si[i] = malloc((2*SHA512_DIGEST_LENGTH+1)*sizeof(char));
	}
	char **commit_ei;
	commit_ei = malloc(u*sizeof(char*));
	for(int i=0; i<u; ++i) {
		commit_ei[i] = malloc((2*SHA512_DIGEST_LENGTH+1)*sizeof(char));
	}
	fmpz_t ***KHj_shares;
	KHj_shares = flint_malloc(u*sizeof(fmpz_t**));
	for(int i=0; i<u; ++i) {
		KHj_shares[i] = flint_malloc(binomt*sizeof(fmpz_t*));
		for(int j=0; j<binomt; ++j) {
			KHj_shares[i][j] = flint_malloc(u*sizeof(fmpz_t));
			for(int k=0; k<u; ++k) fmpz_init(KHj_shares[i][j][k]);
		}
	}
	char ****commit_KHj_shares;
	commit_KHj_shares = malloc(u*sizeof(char***));
	for(int i=0; i<u; ++i) {
		commit_KHj_shares[i] = malloc(binomt*sizeof(char**));
		for(int j=0; j<binomt; ++j) {
			commit_KHj_shares[i][j] = malloc(u*sizeof(char*));
			for(int k=0; k<u; ++k) {
				commit_KHj_shares[i][j][k] = malloc((2*SHA512_DIGEST_LENGTH+1)*sizeof(char));
			}
		}
	}
	fmpz_poly_t **aEj_shares;
	aEj_shares = flint_malloc(u*sizeof(fmpz_poly_t*));
	for(int i=0; i<u; ++i) {
		aEj_shares[i] = flint_malloc(u*sizeof(fmpz_poly_t));
		for(int j=0; j<u; ++j) fmpz_poly_init(aEj_shares[i][j]);
	}
	char ***commit_aEj_shares;
	commit_aEj_shares = malloc(u*sizeof(char**));
	for(int i=0; i<u; ++i) {
		commit_aEj_shares[i] = malloc(u*sizeof(char*));
		for(int j=0; j<u; ++j) {
			commit_aEj_shares[i][j] = malloc((2*SHA512_DIGEST_LENGTH+1)*sizeof(char));
		}
	}
	fmpz_poly_t si[u],ei[u],aux_aE;
	fmpz_poly_init(aux_aE);
	fmpz_t aux_KH;
	fmpz_init(aux_KH);
	fmpz_poly_t auxNIVSS;
	fmpz_poly_init(auxNIVSS);


	//printf("Step 1\n");
	// First step. For every player we sample all values from their respective distributions
	// and then commit them to send them.
	for(int i=0; i<u; ++i) {
		// Begin timer
		clock_t begin = clock();

		//Sample contribution to s and to e from the discrete Gaussian
		fmpz_poly_init(si[i]);
		fmpz_poly_init(ei[i]);
		disc_gauss_Rq(sigmakg,si[i],n,q);
		disc_gauss_Rq(sigmakg,ei[i],n,q);

		// Produce the keys for the NIVSS for both of them
		for(int j=0; j<binomt; ++j) {
			rand_Zq(sNIVSS_keys[i][j],n,q);
			rand_Zq(eNIVSS_keys[i][j],n,q);
		}

		// Compute and commit sj-r and ej-r
		for(int j=0; j<binomt; ++j) {
			fmpz_poly_zero(auxNIVSS);
			PRF(sNIVSS_keys[i][j],lambda,auxNIVSS,interkg,n);
			fmpz_poly_sub(si[i],si[i],auxNIVSS);
			
			fmpz_poly_zero(auxNIVSS);
			PRF(eNIVSS_keys[i][j],lambda,auxNIVSS,interkg,n);
			fmpz_poly_sub(ei[i],ei[i],auxNIVSS);
		}
		
		char *stringsi=fmpz_poly_get_str(si[i]);
		int silen=strlen(stringsi);
		unsigned char hashsi[SHA512_DIGEST_LENGTH];
		commit_sha2(stringsi,silen,hashsi,commit_si[i]);
		free(stringsi);
		char *stringei=fmpz_poly_get_str(ei[i]);
		int eilen=strlen(stringei);
		unsigned char hashei[SHA512_DIGEST_LENGTH];
		commit_sha2(stringei,eilen,hashei,commit_ei[i]);
		free(stringei);
		

		// Sample the KHj, compute Shamir shares and commit them
		for(int j=0; j<binomt; ++j) {
			fmpz_zero(aux_KH);
			rand_Zq(aux_KH,n,q);
			gen_shamir_Zq(KHj_shares[i][j],aux_KH,n,u,t,q);
			for(int k=0; k<u; ++k) {
				char *KHjsstring = fmpz_get_str(NULL,2,KHj_shares[i][j][k]);
				int KHjlen = strlen(KHjsstring);
				unsigned char hashKH[SHA512_DIGEST_LENGTH];
				commit_sha2(KHjsstring,KHjlen,hashKH,commit_KHj_shares[i][j][k]);
				free(KHjsstring);
			}
		}


		// Sample aEj, compute Shamir shares and commit them
		fmpz_poly_zero(aux_aE);
		rand_Rq(aux_aE,n,q);
		gen_shamir_Rq(aEj_shares[i],aux_aE,n,u,t,q);
		for(int j=0; j<u; ++j) {
			char *aEjsstring = fmpz_poly_get_str(aEj_shares[i][j]);
			int aEjlen = strlen(aEjsstring);
			unsigned char hashaEj[SHA512_DIGEST_LENGTH];
			commit_sha2(aEjsstring,aEjlen,hashaEj,commit_aEj_shares[i][j]);
			free(aEjsstring);
		}

		// End timer
		clock_t end = clock();

		// Compute time spent in player i
		times[0][i] = (double)(end-begin)/CLOCKS_PER_SEC;
	}
	

	// Free storage for: aux_KH, aux_aE, auxNIVSS
	fmpz_clear(aux_KH);
	fmpz_poly_clear(aux_aE);
	fmpz_poly_clear(auxNIVSS);



	//printf("Step 2\n");
	// Second step. Every player verifies all commitments sent to him
	for(int i=0; i<u; ++i) {
		// Begin timer
		clock_t begin = clock();

		for(int l=0; l<u; ++l) {
			// Verify commitments of si and ei
			char *stringsi=fmpz_poly_get_str(si[l]);
			int silen=strlen(stringsi);
			int verify_si=verify_commit_sha2(commit_si[l],stringsi,silen);
			if (verify_si==0) printf("ERROR: Commitment of sj done by player %d for player %d does not match\n",l+1,i+1);
			free(stringsi);
			char *stringei=fmpz_poly_get_str(ei[l]);
			int eilen=strlen(stringei);
			int verify_ei=verify_commit_sha2(commit_ei[l],stringei,eilen);
			if (verify_si==0) printf("ERROR: Commitment of ej done by player %d for player %d does not match\n",l+1,i+1);
			free(stringei);


			// Verify commitments of KHj shares
			for(int j=0; j<binomt; ++j) {
				char *KHjsstring = fmpz_get_str(NULL,2,KHj_shares[l][j][i]);
				int KHjlen = strlen(KHjsstring);
				int verify_KHj=verify_commit_sha2(commit_KHj_shares[l][j][i],KHjsstring,KHjlen);
				free(KHjsstring);
				if(verify_KHj==0) printf("ERROR: Commitment of KH number %d done by player %d for player %d does not match\n",j+1,l+1,i+1);
			}


			// Verify commitments of aEj shares
			char *aEjsstring = fmpz_poly_get_str(aEj_shares[l][i]);
			int aEjlen = strlen(aEjsstring);
			int verify_aEj=verify_commit_sha2(commit_aEj_shares[l][i],aEjsstring,aEjlen);
			free(aEjsstring);
			if(verify_aEj==0) printf("ERROR: Commitment of aEj done by player %d for player %d does not match\n",l+1,i+1);
		}

		// End timer
		clock_t end = clock();

		// Compute time spent in player i
		times[1][i] = (double)(end-begin)/CLOCKS_PER_SEC;
	}


	// Free storage for: commit_si, commit_ei, commit_KHj_shares, commit_aEj_shares
	for(int i=0; i<u; ++i) {
		free(commit_si[i]);
		free(commit_ei[i]);
		for(int j=0; j<binomt; ++j) {
			for(int k=0; k<u; ++k) {
				free(commit_KHj_shares[i][j][k]);
			}
			free(commit_KHj_shares[i][j]);
		}
		free(commit_KHj_shares[i]);
		for(int j=0; j<u; ++j) {
			free(commit_aEj_shares[i][j]);
		}
		free(commit_aEj_shares[i]);
	}
	free(commit_si);
	free(commit_ei);
	free(commit_KHj_shares);
	free(commit_aEj_shares);



	// Generate storage for: KH_shares, aE_shares, max, coeff, qfmpz, sNIVSS
	// eNIVSS, eshare
	fmpz_t **KH_shares;
	KH_shares = flint_malloc(binomt*sizeof(fmpz_t*));
	for(int i=0; i<binomt; ++i) {
		KH_shares[i] = flint_malloc(u*sizeof(fmpz_t));
		for(int j=0; j<u; ++j) fmpz_init(KH_shares[i][j]);
	}
	fmpz_poly_t aE_shares[u];
	fmpz_t max,coeff,qfmpz;
	fmpz_init(max);
	fmpz_set_mpz(max,interkg);
	fmpz_mul_ui(max,max,binomt);
	fmpz_init(coeff);
	fmpz_init(qfmpz);
	fmpz_set_mpz(qfmpz,q);
	fmpz_poly_t sNIVSS,eNIVSS;
	fmpz_poly_init(sNIVSS);
	fmpz_poly_init(eNIVSS);
	fmpz_poly_t eshare[u];
	for(int i=0; i<u; ++i) fmpz_poly_init(eshare[i]);

	//printf("Step 3\n");
	// Third step. For s and e compute the shamir shares of each contribution 
	// following the NIVSS protocol and add them all to get their share 
	// on s and e. For KH and aE every player adds all the shares he received 
	// and sends the share to the appropriate players
	for(int i=0; i<u; ++i) {
		// Begin timer
		clock_t begin = clock();

		// Verify that every sj-r and ej-r are in the interval needed
		for(int k=0; k<u; ++k) {
			for(int l=0; l<n; ++l) {
				fmpz_zero(coeff);
				fmpz_poly_get_coeff_fmpz(coeff,si[k],l);
				if(fmpz_cmp(coeff,max)>0) printf("ERROR: coefficient %d of masked contribution for s of player %d is too big.\n",l+1,k+1);

				fmpz_zero(coeff);
				fmpz_poly_get_coeff_fmpz(coeff,ei[k],l);
				if(fmpz_cmp(coeff,max)>0) printf("ERROR: coefficient %d of masked contribution for e of player %d is too big.\n",l+1,k+1);
			}
		}

		// For e and s we compute their Shamir share
		for(int k=0; k<u; ++k) {
			// For every player compute the NIVSS contribution for s and e
			fmpz_poly_zero(sNIVSS);
			PRSS_share(sNIVSS,key_map1,allowedt,matrixt,sNIVSS_keys[k],lambda,interkg,n,u,t,i+1,q);
			
			fmpz_poly_zero(eNIVSS);
			PRSS_share(eNIVSS,key_map1,allowedt,matrixt,eNIVSS_keys[k],lambda,interkg,n,u,t,i+1,q);
			
			// We add that to si[k] and ei[k] to get our share on the contribution 
			// and add it to the total to get the general share
			fmpz_poly_add(sshare[i],sshare[i],sNIVSS);
			fmpz_poly_add(sshare[i],sshare[i],si[k]);
			fmpz_p_mod(sshare[i],n,sshare[i],q);
			
			fmpz_poly_add(eshare[i],eshare[i],eNIVSS);
			fmpz_poly_add(eshare[i],eshare[i],ei[k]);
			fmpz_p_mod(eshare[i],n,eshare[i],q);
		}


		// We add all the KHj shares received
		for(int j=0; j<binomt; ++j) {
			for(int k=0; k<u; ++k) {
				fmpz_add(KH_shares[j][i],KH_shares[j][i],KHj_shares[k][j][i]);
				fmpz_mod(KH_shares[j][i],KH_shares[j][i],qfmpz);
			}
		}


		// We add all the aEj shares received
		fmpz_poly_init(aE_shares[i]);
		for(int j=0; j<u; ++j) {
			fmpz_poly_add(aE_shares[i],aE_shares[i],aEj_shares[j][i]);
			fmpz_p_mod(aE_shares[i],n,aE_shares[i],q);
		}

		// End timer
		clock_t end = clock();

		// Compute time spent in player i
		times[2][i] = (double)(end-begin)/CLOCKS_PER_SEC;
	}

	
	// Free storage for: max, coeff, qfmpz, si, ei, sNIVSS, eNIVSS, sNIVSS_keys, 
	// eNIVSS_keys, KHj_shares, aEj_shares
	fmpz_clear(max);
	fmpz_clear(coeff);
	fmpz_clear(qfmpz);
	fmpz_poly_clear(sNIVSS);
	fmpz_poly_clear(eNIVSS);
	for(int i=0; i<u; ++i) {
		fmpz_poly_clear(si[i]);
		fmpz_poly_clear(ei[i]);
		for(int j=0; j<binomt; ++j) {
			fmpz_clear(sNIVSS_keys[i][j]);
			fmpz_clear(eNIVSS_keys[i][j]);
			for(int k=0; k<u; ++k) {
				fmpz_clear(KHj_shares[i][j][k]);
			}

			flint_free(KHj_shares[i][j]);
		}
		flint_free(sNIVSS_keys[i]);
		flint_free(eNIVSS_keys[i]);
		flint_free(KHj_shares[i]);
		for(int j=0; j<u; ++j) {
			fmpz_poly_clear(aEj_shares[i][j]);
		}
		flint_free(aEj_shares[i]);
			
	}

	flint_free(sNIVSS_keys);
	flint_free(eNIVSS_keys);
	flint_free(KHj_shares);
	flint_free(aEj_shares);


	// Generate storage for: bE_share, auxKH, auxaE
	fmpz_poly_t bE_share[u];
	for(int i=0; i<u; ++i) fmpz_poly_init(bE_share[i]);
	fmpz_t **auxKH;
	auxKH = flint_malloc(binomt*sizeof(fmpz_t*));
	for(int j=0; j<binomt; ++j) {
		auxKH[j] = flint_malloc(binomt1*sizeof(fmpz_t));
		for(int l=0; l<binomt1; ++l) fmpz_init(auxKH[j][l]);
	}
	fmpz_poly_t auxaE[binomt1];
	for(int j=0; j<binomt1; ++j) fmpz_poly_init(auxaE[j]);


	//printf("Step 4\n");
	// Fourth step. We reconstruct KH and aE and we compute shares for bE
	for(int i=0; i<u; ++i) {
		// Begin timer
		clock_t begin = clock();

		// Retrieve KH
		int binomut1 = binomial(u-1,u-t-1);
		for(int j=0; j<binomut1; ++j) {
			int column = allowedt[i][j];
			int compare_KH = 1;
			for(int l=0; l<binomt1; ++l) {
				int players[t+1];
				int count = 0;
				for(int k=0; k<u; ++k) {
					if(matrixt1[k][l]==1) {
						players[count]=k+1;
						count=count+1;
					}
				}
				compute_shamir_Zq(auxKH[column][l],KH_shares[column],players,t,u,n,q);
				if(fmpz_equal(auxKH[column][l],auxKH[column][0])==0) compare_KH = 0;
			}
			//printf("Returns 1 if all KH are equal, 0 otherwise: %d\n",compare_KH);
			fmpz_zero(KH[column]);
			if(compare_KH==1) fmpz_set(KH[column],auxKH[column][0]);
		}


		// Retrieve aE
		int compare_aE = 1;
		for(int j=0; j<binomt1; ++j) {
			int players[t+1];
			int count = 0;
			for(int k=0; k<u; ++k) {
				if(matrixt1[k][j]==1) {
					players[count]=k+1;
					count=count+1;
				}
			}
			compute_shamir_Rq(auxaE[j],aE_shares,players,t,u,n,q);
			if(fmpz_poly_equal(auxaE[j],auxaE[0])==0) compare_aE = 0;
		}
		//printf("Returns 1 if all aE are equal, 0 otherwise: %d\n",compare_aE);
		if(compare_aE==1) fmpz_poly_set(aE,auxaE[0]);
		fmpz_p_mod(aE,n,aE,q);


		// Compute the shares of bE
		fmpz_poly_mul_Rq(bE_share[i],aE,sshare[i],n,q);
		fmpz_poly_add(bE_share[i],bE_share[i],eshare[i]);
		fmpz_p_mod(bE_share[i],n,bE_share[i],q);

		// End timer
		clock_t end = clock();

		// Compute time spent in player i
		times[3][i] = (double)(end-begin)/CLOCKS_PER_SEC;
	}


	// Free storage for: KH_shares, auxKH, auxaE, eshare
	for(int i=0; i<binomt; ++i) {
		for(int j=0; j<u; ++j) {
			fmpz_clear(KH_shares[i][j]);
		}
		flint_free(KH_shares[i]);
	}
	flint_free(KH_shares);
	for(int j=0; j<binomt; ++j) {
		for(int k=0; k<binomt1; ++k) fmpz_clear(auxKH[j][k]);
		flint_free(auxKH[j]);
	}
	flint_free(auxKH);
	for(int i=0; i<u; ++i) {
		fmpz_poly_clear(aE_shares[i]);
	}
	for(int j=0; j<binomt1; ++j) fmpz_poly_clear(auxaE[j]);
	for(int i=0; i<u; ++i) fmpz_poly_clear(eshare[i]);

	
	// Create storage: auxbE
	fmpz_poly_t auxbE[binomt1];
	for(int i=0; i<binomt1; ++i) fmpz_poly_init(auxbE[i]);

	//printf("Step 5\n");
	// Fifth step. Recover bE
	for(int i=0; i<u; ++i) {
		// Begin timer
		clock_t begin = clock();

		int compare_bE = 1;
		for(int j=0; j<binomt1; ++j) {
			int players[t+1];
			int count = 0;
			for(int k=0; k<u; ++k) {
				if(matrixt1[k][j]==1) {
					players[count]=k+1;
					count=count+1;
				}
			}
			compute_shamir_Rq(auxbE[j],bE_share,players,t,u,n,q);
			if(fmpz_poly_equal(auxbE[j],auxbE[0])==0) compare_bE = 0;
		}
		//printf("Returns 1 if all bE are equal, 0 otherwise: %d\n",compare_bE);
		if(compare_bE==1) fmpz_poly_set(bE,auxbE[0]);
		fmpz_p_mod(bE,n,bE,q);

		// End timer
		clock_t end = clock();

		// Compute time spent in player i
		times[4][i] = (double)(end-begin)/CLOCKS_PER_SEC;
	}


	// We free storage for:  bE_share, auxbE
	for(int i=0; i<u; ++i) fmpz_poly_clear(bE_share[i]);
	for(int i=0; i<binomt1; ++i) fmpz_poly_clear(auxbE[i]);
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

	//We set the general parameters sigma, interkg and seed
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
	
	// Compute interkg
	mpz_t interkg;
	mpz_init(interkg);
	mpz_ui_pow_ui(auxz,2,100);
	mpz_mul(interkg,kappa,auxz);
	mpz_tdiv_q_ui(interkg,interkg,binomt);
	//gmp_printf("interkg = %Zd\n",interkg);
	
	//printf("\n");
	
	mpf_clear(auxf);
	mpf_clear(aux2f);
	mpf_clear(qf);
	mpz_clear(auxz);
	
	//Set seed
	srand(time(NULL));

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

	//Compute the matrix of t+1 subsets of players
	int binomt1 = binomial(u,t+1);
	int** matrixt1 = (int**)malloc(u * sizeof(int*));
	for (int index=0;index<u;++index){
	matrixt1[index] = (int*)malloc(binomt1 * sizeof(int));
	for(int j=0; j<binomt1;++j) matrixt1[index][j]=0;
	}
	matrix(u,t+1,binomt1,matrixt1);


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


	//printf("Beginning of the for\n");
	// Simulation of Key Generation
	fmpz_poly_t lambda;
	fmpz_poly_init(lambda);
	//Set the argument of the PRF as 0,1,2,3,....
	for(int i=0; i<n;++i) fmpz_poly_set_coeff_si(lambda,i,i);
	fmpz_poly_t s;
	fmpz_poly_t aE;
	fmpz_poly_t bE;
	fmpz_poly_init(s);
	fmpz_poly_init(aE);
	fmpz_poly_init(bE);
	fmpz_t KH[keys];
	for(int i=0; i<keys; ++i) fmpz_init(KH[i]);
	fmpz_poly_t sshare[u];
	for(int i=0; i<u; ++i) fmpz_poly_init(sshare[i]);
	double **timeskg = (double**)malloc(5 * sizeof(double*));
	for (int index=0;index<5;++index){
		timeskg[index] = (double*)malloc(u * sizeof(double));
		for(int j=0; j<u;++j) timeskg[index][j]=0.0;
	}


	// Simulate the key generation repetitions times
	fmpz_poly_t m1, m2, uenc, v;
	fmpz_poly_init(m1);
	fmpz_poly_init(m2);
	fmpz_poly_init(uenc);
	fmpz_poly_init(v);
	double avgmin=0.0;
	double avgmax=0.0;
	for(int rep=0; rep<repetitions; ++rep) {
		// Initialize again all needed variables: s,aE,bE,KH,sshares,timeskg
		fmpz_poly_zero(s);
		fmpz_poly_zero(aE);
		fmpz_poly_zero(bE);
		for(int i=0; i<keys; ++i) fmpz_zero(KH[i]);
		for (int index=0;index<5;++index){
			for(int j=0; j<u;++j) timeskg[index][j]=0.0;
		}

		// Simulate Key Generation
		keygen_sim(timeskg,sshare,aE,bE,KH,allowedt,lambda,matrixt1,matrixt,key_map1,keys,sigma,binomt1,u,n,t,q,interkg);

		// Recover s
		int players[t+1];
		int count = 0;
		for(int k=0; k<u; ++k) {
			if(matrixt1[k][0]==1) {
				players[count]=k+1;
				count=count+1;
			}
		}
		compute_shamir_Rq(s,sshare,players,t,u,n,q);

		// Verify several messages are correctly encrypted and decrypted
		fmpz_poly_zero(m1);
		fmpz_poly_zero(m2);
		fmpz_poly_zero(uenc);
		fmpz_poly_zero(v);
		for(int rep=0; rep<20; ++rep) {
			// Generate random message m1
			for(int i=0; i<n; ++i) {
				double random = rand_gen();
				if(random<0.5) fmpz_poly_set_coeff_ui(m1,i,0);
				else fmpz_poly_set_coeff_ui(m1,i,1);
			}
			//Encrypt message
			encrypt(m1,uenc,v,aE,bE,sigma,n,q,s);

			// Decrypt message
			decrypt(m2,s,uenc,v,q,n);
			int compare_messages = fmpz_poly_equal(m1,m2);
			if(compare_messages==0) printf("ERROR: message incorrectly decrypted\n\n");
		}

		// Compute minimum and maximum time and add it to the Average
		for(int i = 0; i<5; ++i) {
			double min = 1000000000000000000.0;
			double max = 0.0;
			for(int j=0; j<u; ++j) {
				if(timeskg[i][j]<min) min=timeskg[i][j];
				if(timeskg[i][j]>max) max=timeskg[i][j];
			}
			avgmin=avgmin+min;
			avgmax=avgmax+max;
		}
		//printf("The maximum average time in Key Generation was %fms in iteration %d\n", avgmax*1000/(rep+1), rep+1);
	}

	// Compute the average and print it in miliseconds
	avgmin = avgmin/repetitions*1000;
	avgmax = avgmax/repetitions*1000;
	FILE *f = fopen("keygen.csv","a");
	fprintf(f,"%f\n",avgmax);
	fclose(f);
	//printf("The minimum average time in Key Generation was %fms\n", avgmin);
	//printf("The maximum average time in Key Generation was %fms\n", avgmax);


	// Free all storage
	mpz_clear(q);
	mpz_clear(qp);
	mpz_clear(interkg);
	mpf_clear(sigma);
	for(int index=0; index<u; ++index) free(matrixt1[index]);
	free(matrixt1);
	for(int i=0; i<u; ++i) free(matrixt[i]);
	free(matrixt);
	for(int index=0; index<u; ++index) free(allowedt[index]);
	free(allowedt);
	for(int index=0; index<u+keys; ++index) free(key_map1[index]);
	free(key_map1);
	for(int i=0; i<keys; ++i) fmpz_clear(KH[i]);
	for(int i=0; i<u; ++i) fmpz_poly_clear(sshare[i]);
	fmpz_poly_clear(lambda);
	for (int index=0;index<5;++index) free(timeskg[index]);
	free(timeskg);
}
