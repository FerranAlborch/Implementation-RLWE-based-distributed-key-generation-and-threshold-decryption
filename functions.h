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


void mpf_log(mpf_t result, mpf_t element);

void fmpz_p_mod(fmpz_poly_t a, int l, fmpz_poly_t b, mpz_t q);

void fmpz_poly_mul_Rq(fmpz_poly_t c, fmpz_poly_t a, fmpz_poly_t b, int n,
												mpz_t q);

void pop_back(int vector[], int new_vector[], int length);

void push_back(int vector[], int new_vector[], int value, int length);

int cmpfunc (const void * a, const void * b);

int equal_vec(int* vec1, int* vec2, int size);

int smaller_vec(int* vec1, int* vec2, int size);

int find_keys(int** key_map1, int* key, int size, int start, int end);

void create_keys(int state[], int indexes[], int size, int** M,
              int length_state, int length_indexes, int keys, int t);

int binomial(int n, int k);

void subsets(int state[], int indexes[], int size, int** M,
              int length_state, int length_indexes, int u, int binom);

void matrix(int u, int t, int binom, int** M);

double rand_gen();

void normal_sampler(mpf_t sigma, mpf_t sample);

void round_mpf(mpz_t near, mpf_t x);

void disc_gauss_Zq(mpf_t sigma, mpz_t q, fmpz_t sample);

void disc_gauss_Rq(mpf_t sigma, fmpz_poly_t sample, int n, mpz_t q);

void rand_Rq(fmpz_poly_t c, int n, mpz_t q);

void rand_Zq(fmpz_t c, int n, mpz_t q);

void inverse_Zq(mpz_t inverse, mpz_t a, mpz_t q);

void gen_shamir_Rq(fmpz_poly_t shares[], fmpz_poly_t element, int n, int u, int t,
							mpz_t q);
							
void gen_shamir_Zq(fmpz_t shares[], fmpz_t element, int n, int u, int t,
							mpz_t q);

void compute_shamir_Rq(fmpz_poly_t element, fmpz_poly_t shares[], int players[],
												 int t, int u, int n, mpz_t q);
              				  
void compute_shamir_Zq(fmpz_t element, fmpz_t shares[], int players[],
												 int t, int u, int n, mpz_t q);
												 
void compute_polynomials(mpz_t eval, int *players, int z, int t, mpz_t q);

unsigned char *mx_hmac_sha3_512(const void *key, int keylen,
                                const unsigned char *data, int datalen,
                                unsigned char *result, unsigned int *resultlen);

void PRF(fmpz_t key, fmpz_poly_t lambda, fmpz_poly_t result, mpz_t inter, int n);

void PRSS_share(fmpz_poly_t result, int **key_map, int **allowedt, 
			     int **matrixt, fmpz_t KH[], fmpz_poly_t lambda, 
			     mpz_t inter, int n, int u, int t, int player, mpz_t q);

void commit_sha2(char *message, int message_len, unsigned char hash[],
                  char shash[]);

int verify_commit_sha2(char commitment[], char *message, int message_len);

void encrypt(fmpz_poly_t m, fmpz_poly_t u, fmpz_poly_t v, fmpz_poly_t aE,
							fmpz_poly_t bE, mpf_t sigma, int n, mpz_t q,fmpz_poly_t s);
							
void round_message(fmpz_poly_t e, fmpz_poly_t m, mpz_t q, int n);
