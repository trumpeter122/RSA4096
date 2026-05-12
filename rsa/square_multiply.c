#include <string.h>
#include "compute.h"

#define WINDOW_BITS 4
#define WINDOW_SIZE (1U << WINDOW_BITS)

static uint32_t highest_bit_index(bn_t *a, uint32_t digits);
static uint32_t exponent_bit(bn_t *a, uint32_t bit);
static void binary_mod_exp(bn_t *a, bn_t *base, bn_t e, bn_t *d, uint32_t ddigits);

void rsa_bn_mod_exp(bn_t *a, bn_t *b, bn_t *c, uint32_t cdigits, bn_t *d, uint32_t ddigits)
{
    bn_t powers[WINDOW_SIZE-1][BN_MAX_DIGITS], t[BN_MAX_DIGITS], base[BN_MAX_DIGITS];
    int bit, top_bit, window_bits;
    uint32_t i, k, value;

    BN_ASSIGN_DIGIT(t, 1, ddigits);
    if(bn_cmp(b, d, ddigits) >= 0) {
        bn_mod(base, b, ddigits, d, ddigits);
    } else {
        bn_assign(base, b, ddigits);
    }

    cdigits = bn_digits(c, cdigits);
    if(cdigits == 0) {
        bn_assign(a, t, ddigits);
        memset((uint8_t *)t, 0, sizeof(t));
        memset((uint8_t *)base, 0, sizeof(base));
        return;
    }

    if(cdigits == 1) {
        binary_mod_exp(a, base, c[0], d, ddigits);
        memset((uint8_t *)t, 0, sizeof(t));
        memset((uint8_t *)base, 0, sizeof(base));
        return;
    }

    bn_assign(powers[0], base, ddigits);
    for(i=1; i<WINDOW_SIZE-1; i++) {
        bn_mod_mul(powers[i], powers[i-1], base, d, ddigits);
    }

    top_bit = (int)highest_bit_index(c, cdigits);
    bit = top_bit;
    window_bits = bit % WINDOW_BITS + 1;
    while(bit >= 0) {
        value = 0;
        for(k=0; k<(uint32_t)window_bits; k++) {
            value = (value << 1) | exponent_bit(c, (uint32_t)(bit - (int)k));
        }

        if(bit == top_bit) {
            bn_assign(t, powers[value-1], ddigits);
        } else {
            for(k=0; k<(uint32_t)window_bits; k++) {
                bn_mod_mul(t, t, t, d, ddigits);
            }
            if(value != 0) {
                bn_mod_mul(t, t, powers[value-1], d, ddigits);
            }
        }

        bit -= window_bits;
        window_bits = WINDOW_BITS;
    }

    bn_assign(a, t, ddigits);

    // Clear potentially sensitive information
    memset((uint8_t *)powers, 0, sizeof(powers));
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

static uint32_t exponent_bit(bn_t *a, uint32_t bit)
{
    return (a[bit / BN_DIGIT_BITS] >> (bit % BN_DIGIT_BITS)) & 1U;
}

static void binary_mod_exp(bn_t *a, bn_t *base, bn_t e, bn_t *d, uint32_t ddigits)
{
    bn_t t[BN_MAX_DIGITS];
    int bit;

    bit = BN_DIGIT_BITS - 1;
    while(bit > 0 && ((e >> bit) & 1U) == 0) {
        bit--;
    }

    bn_assign(t, base, ddigits);
    for(bit--; bit>=0; bit--) {
        bn_mod_mul(t, t, t, d, ddigits);
        if((e >> bit) & 1U) {
            bn_mod_mul(t, t, base, d, ddigits);
        }
    }

    bn_assign(a, t, ddigits);

    // Clear potentially sensitive information
    memset((uint8_t *)t, 0, sizeof(t));
}
