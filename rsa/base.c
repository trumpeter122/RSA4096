#include <string.h>
#include "compute.h"

void rsa_bn_mod_exp(bn_t *a, bn_t *b, bn_t *c, uint32_t cdigits, bn_t *d, uint32_t ddigits)
{
    bn_t bpower[3][BN_MAX_DIGITS], ci, t[BN_MAX_DIGITS];
    int i;
    uint32_t ci_bits, j, s;

    bn_assign(bpower[0], b, ddigits);
    bn_mod_mul(bpower[1], bpower[0], b, d, ddigits);
    bn_mod_mul(bpower[2], bpower[1], b, d, ddigits);

    BN_ASSIGN_DIGIT(t, 1, ddigits);

    cdigits = bn_digits(c, cdigits);
    i = cdigits - 1;
    for(; i>=0; i--) {
        ci = c[i];
        ci_bits = BN_DIGIT_BITS;

        if(i == (int)(cdigits - 1)) {
            while(!DIGIT_2MSB(ci)) {
                ci <<= 2;
                ci_bits -= 2;
            }
        }

        for(j=0; j<ci_bits; j+=2) {
            bn_mod_mul(t, t, t, d, ddigits);
            bn_mod_mul(t, t, t, d, ddigits);
            if((s = DIGIT_2MSB(ci)) != 0) {
                bn_mod_mul(t, t, bpower[s-1], d, ddigits);
            }
            ci <<= 2;
        }
    }

    bn_assign(a, t, ddigits);

    // Clear potentially sensitive information
    memset((uint8_t *)bpower, 0, sizeof(bpower));
    memset((uint8_t *)t, 0, sizeof(t));
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
