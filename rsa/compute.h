#ifndef __RSA_COMPUTE_H__
#define __RSA_COMPUTE_H__

#include "../rsa.h"
#include "../bignum.h"

void rsa_bn_mod_exp(bn_t *a, bn_t *b, bn_t *c, uint32_t cdigits, bn_t *d, uint32_t ddigits);
int rsa_private_compute(bn_t *out, bn_t *in, rsa_sk_t *sk, bn_t *n, uint32_t ndigits);
int rsa_public_compute(bn_t *out, bn_t *in, rsa_pk_t *pk, bn_t *n, uint32_t ndigits);

#endif  // __RSA_COMPUTE_H__
