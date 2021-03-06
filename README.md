# Implementation: RLWE-based distributed key generation and threshold decryption

This repository contains source code of the protocols from the paper *RLWE-based distributed key generation and threshold decryption* by Ferran Alborch, Ramiro Martínez and Paz Morillo, [https://doi.org/10.3390/math10050728](https://doi.org/10.3390/math10050728).

## Run

There are three C files, one for each simulation of encryption, decryption and key generation, together with auxiliary functions. To compile and execute these files on Linux, run
```
gcc encryption_sim.c functions.c -lm -lgmp -lssl -lcrypto -lflint -lmpfr -o encrypt.out -O2
./encrypt.out n u t

gcc decryption_sim.c functions.c -lm -lgmp -lssl -lcrypto -lflint -lmpfr -o decrypt.out -O2
./decrypt.out n u t

gcc keygen_sim.c functions.c -lm -lgmp -lssl -lcrypto -lflint -lmpfr -o keygen.out -O2
./keygen.out n u t
```
respectively, where n, u, t are the non-negative integer values you want to give to these parameters. To run the codes a large number of times the use of a python script is advisable.

To be able to run these programs the following libraries are needed:
- FLINT (Fast Library for Number Theory), found [here](https://www.flintlib.org/downloads.html). The version used is 2.7.1. FLINT requires two other libraries:
  - GMP, found [here](https://gmplib.org/). The version used is 6.2.1.
  - MPFR, found [here](https://www.mpfr.org/). The version used is 4.1.0.
- OpenSSL, found [here](https://www.openssl.org/). The version used is 1.1.1.

<br/>
There are also two files that allow for the parameter analysis in python. One is an ipynb file and the other is needed to be able to execute the code in the notebook.

<br/><br/>
WARNING: This is an academic proof of concept, and in particular has not received code review. This implementation is NOT ready for any type of production use.
