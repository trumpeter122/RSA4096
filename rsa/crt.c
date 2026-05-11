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
    bn_t p[BN_MAX_DIGITS], q[BN_MAX_DIGITS], dp[BN_MAX_DIGITS], dq[BN_MAX_DIGITS];
    bn_t qinv[BN_MAX_DIGITS], cp[BN_MAX_DIGITS], cq[BN_MAX_DIGITS];
    bn_t m1[BN_MAX_DIGITS], m2[BN_MAX_DIGITS], diff[BN_MAX_DIGITS], h[BN_MAX_DIGITS];
    bn_t qh[2*BN_MAX_DIGITS], m2n[BN_MAX_DIGITS];
    uint32_t pdigits, qdigits;

    bn_decode(p, BN_MAX_DIGITS, sk->prime1, RSA_MAX_PRIME_LEN);
    bn_decode(q, BN_MAX_DIGITS, sk->prime2, RSA_MAX_PRIME_LEN);
    bn_decode(dp, BN_MAX_DIGITS, sk->prime_exponent1, RSA_MAX_PRIME_LEN);
    bn_decode(dq, BN_MAX_DIGITS, sk->prime_exponent2, RSA_MAX_PRIME_LEN);
    bn_decode(qinv, BN_MAX_DIGITS, sk->coefficient, RSA_MAX_PRIME_LEN);

    pdigits = bn_digits(p, BN_MAX_DIGITS);
    qdigits = bn_digits(q, BN_MAX_DIGITS);

    bn_mod(cp, in, ndigits, p, pdigits);
    bn_mod(cq, in, ndigits, q, qdigits);
    rsa_bn_mod_exp(m1, cp, dp, bn_digits(dp, pdigits), p, pdigits);
    rsa_bn_mod_exp(m2, cq, dq, bn_digits(dq, qdigits), q, qdigits);

    if(bn_cmp(m1, m2, pdigits) >= 0) {
        bn_sub(diff, m1, m2, pdigits);
    } else {
        bn_sub(diff, m2, m1, pdigits);
        bn_sub(diff, p, diff, pdigits);
    }

    bn_mod_mul(h, diff, qinv, p, pdigits);
    bn_assign_zero(qh, 2*BN_MAX_DIGITS);
    bn_mul(qh, q, h, pdigits);

    bn_assign_zero(m2n, ndigits);
    bn_assign(m2n, m2, qdigits);
    bn_add(out, qh, m2n, ndigits);
    if(bn_cmp(out, n, ndigits) >= 0) {
        bn_sub(out, out, n, ndigits);
    }

    // Clear potentially sensitive information
    memset((uint8_t *)p, 0, sizeof(p));
    memset((uint8_t *)q, 0, sizeof(q));
    memset((uint8_t *)dp, 0, sizeof(dp));
    memset((uint8_t *)dq, 0, sizeof(dq));
    memset((uint8_t *)qinv, 0, sizeof(qinv));
    memset((uint8_t *)cp, 0, sizeof(cp));
    memset((uint8_t *)cq, 0, sizeof(cq));
    memset((uint8_t *)m1, 0, sizeof(m1));
    memset((uint8_t *)m2, 0, sizeof(m2));
    memset((uint8_t *)diff, 0, sizeof(diff));
    memset((uint8_t *)h, 0, sizeof(h));
    memset((uint8_t *)qh, 0, sizeof(qh));
    memset((uint8_t *)m2n, 0, sizeof(m2n));

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
