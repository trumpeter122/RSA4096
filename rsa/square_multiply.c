#include <string.h>
#include "compute.h"

static uint32_t highest_bit_index(bn_t *a, uint32_t digits);

void rsa_bn_mod_exp(bn_t *a, bn_t *b, bn_t *c, uint32_t cdigits, bn_t *d, uint32_t ddigits)
{
    bn_t t[BN_MAX_DIGITS], base[BN_MAX_DIGITS];
    int bit;

    BN_ASSIGN_DIGIT(t, 1, ddigits);
    bn_assign(base, b, ddigits);

    cdigits = bn_digits(c, cdigits);
    if(cdigits != 0) {
        bit = (int)highest_bit_index(c, cdigits);
        for(; bit>=0; bit--) {
            uint32_t word = (uint32_t)bit / BN_DIGIT_BITS;
            uint32_t shift = (uint32_t)bit % BN_DIGIT_BITS;

            bn_mod_mul(t, t, t, d, ddigits);
            if((c[word] >> shift) & 1) {
                bn_mod_mul(t, t, base, d, ddigits);
            }
        }
    }

    bn_assign(a, t, ddigits);

    // Clear potentially sensitive information
    memset((uint8_t *)t, 0, sizeof(t));
    memset((uint8_t *)base, 0, sizeof(base));
}

int rsa_private_compute(bn_t *out, bn_t *in, rsa_sk_t *sk, bn_t *n, uint32_t ndigits)
{
    bn_t d[BN_MAX_DIGITS];
    uint32_t ddigits;

    bn_decode(d, BN_MAX_DIGITS, sk->exponent, RSA_MAX_MODULUS_LEN);
    ddigits = bn_digits(d, BN_MAX_DIGITS);
    rsa_bn_mod_exp(out, in, d, ddigits, n, ndigits);

    // Clear potentially sensitive information
    memset((uint8_t *)d, 0, sizeof(d));

    return 0;
}

int rsa_public_compute(bn_t *out, bn_t *in, rsa_pk_t *pk, bn_t *n, uint32_t ndigits)
{
    bn_t e[BN_MAX_DIGITS];
    uint32_t edigits;

    bn_decode(e, BN_MAX_DIGITS, pk->exponent, RSA_MAX_MODULUS_LEN);
    edigits = bn_digits(e, BN_MAX_DIGITS);
    rsa_bn_mod_exp(out, in, e, edigits, n, ndigits);

    // Clear potentially sensitive information
    memset((uint8_t *)e, 0, sizeof(e));

    return 0;
}

static uint32_t highest_bit_index(bn_t *a, uint32_t digits)
{
    bn_t top;
    uint32_t bit;

    top = a[digits-1];
    bit = (digits - 1) * BN_DIGIT_BITS;
    while(top >>= 1) {
        bit++;
    }

    return bit;
}
