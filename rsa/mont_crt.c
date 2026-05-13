#include <string.h>
#include "compute.h"
#include "../mul/mul.h"

#define MONT_WINDOW_BITS 6
#define MONT_WINDOW_TABLE (1U << (MONT_WINDOW_BITS - 1))
#define MONT_WINDOW_MIN_BITS 256

typedef struct {
    rsa_sk_t *sk;
    bn_t p[BN_MAX_DIGITS], q[BN_MAX_DIGITS], dp[BN_MAX_DIGITS], dq[BN_MAX_DIGITS];
    bn_t qinv[BN_MAX_DIGITS], rp[BN_MAX_DIGITS], rrp[BN_MAX_DIGITS], rq[BN_MAX_DIGITS], rrq[BN_MAX_DIGITS];
    bn_t p_n0inv, q_n0inv;
    uint32_t pdigits, qdigits, dpdigits, dqdigits;
} mont_private_cache_t;

typedef struct {
    rsa_pk_t *pk;
    bn_t n[BN_MAX_DIGITS], e[BN_MAX_DIGITS], r[BN_MAX_DIGITS], rr[BN_MAX_DIGITS], one[BN_MAX_DIGITS];
    bn_t n0inv;
    uint32_t ndigits, edigits;
} mont_public_cache_t;

static bn_t montgomery_n0inv(bn_t n0);
static uint32_t exponent_bit(bn_t *a, uint32_t bit);
static uint32_t highest_bit_index(bn_t *a, uint32_t digits);
static void montgomery_mul(bn_t *a, bn_t *b, bn_t *c, bn_t *n, uint32_t digits, bn_t n0inv);
static void montgomery_precompute_constants(bn_t *r, bn_t *rr, bn_t *n, uint32_t digits);
static void montgomery_mod_mul_precomp(bn_t *a, bn_t *b, bn_t *c, bn_t *n,
                                       uint32_t digits, bn_t n0inv, bn_t *rr);
static void montgomery_mod_exp_precomp(bn_t *a, bn_t *b, bn_t *c, uint32_t cdigits,
                                       bn_t *n, uint32_t ndigits, bn_t n0inv, bn_t *r, bn_t *rr);
static void montgomery_mod_exp_binary(bn_t *a, bn_t *base_m, bn_t *c, uint32_t cdigits,
                                      bn_t *n, uint32_t ndigits, bn_t n0inv);
static void montgomery_mod_exp_window(bn_t *a, bn_t *base_m, bn_t *c, uint32_t cdigits,
                                      bn_t *n, uint32_t ndigits, bn_t n0inv, bn_t *r);
static void montgomery_mod_exp_65537(bn_t *a, bn_t *b, bn_t *n, uint32_t ndigits,
                                     bn_t n0inv, bn_t *rr, bn_t *one);
static void classic_mod_exp(bn_t *a, bn_t *b, bn_t *c, uint32_t cdigits, bn_t *n, uint32_t ndigits);
static mont_private_cache_t *get_private_cache(rsa_sk_t *sk);
static mont_public_cache_t *get_public_cache(rsa_pk_t *pk);

void rsa_bn_mod_exp(bn_t *a, bn_t *b, bn_t *c, uint32_t cdigits, bn_t *d, uint32_t ddigits)
{
    bn_t r[BN_MAX_DIGITS], rr[BN_MAX_DIGITS];
    bn_t n0inv;

    if(ddigits == 0 || d[0] == 0 || ((d[0] & 1) == 0)) {
        classic_mod_exp(a, b, c, cdigits, d, ddigits);
        return;
    }

    n0inv = montgomery_n0inv(d[0]);
    montgomery_precompute_constants(r, rr, d, ddigits);
    montgomery_mod_exp_precomp(a, b, c, cdigits, d, ddigits, n0inv, r, rr);

    // Clear potentially sensitive information
    memset((uint8_t *)r, 0, sizeof(r));
    memset((uint8_t *)rr, 0, sizeof(rr));
}

int rsa_private_compute(bn_t *out, bn_t *in, rsa_sk_t *sk, bn_t *n, uint32_t ndigits)
{
    mont_private_cache_t *cache;
    bn_t cp[BN_MAX_DIGITS], cq[BN_MAX_DIGITS];
    bn_t m1[BN_MAX_DIGITS], m2[BN_MAX_DIGITS], diff[BN_MAX_DIGITS], h[BN_MAX_DIGITS];
    bn_t qh[2*BN_MAX_DIGITS], m2n[BN_MAX_DIGITS];

    cache = get_private_cache(sk);

    bn_mod(cp, in, ndigits, cache->p, cache->pdigits);
    bn_mod(cq, in, ndigits, cache->q, cache->qdigits);
    montgomery_mod_exp_precomp(m1, cp, cache->dp, cache->dpdigits, cache->p, cache->pdigits,
                               cache->p_n0inv, cache->rp, cache->rrp);
    montgomery_mod_exp_precomp(m2, cq, cache->dq, cache->dqdigits, cache->q, cache->qdigits,
                               cache->q_n0inv, cache->rq, cache->rrq);

    if(bn_cmp(m1, m2, cache->pdigits) >= 0) {
        bn_sub(diff, m1, m2, cache->pdigits);
    } else {
        bn_sub(diff, m2, m1, cache->pdigits);
        bn_sub(diff, cache->p, diff, cache->pdigits);
    }

    montgomery_mod_mul_precomp(h, diff, cache->qinv, cache->p, cache->pdigits,
                               cache->p_n0inv, cache->rrp);
    bn_assign_zero(qh, 2*BN_MAX_DIGITS);
    bn_mul(qh, cache->q, h, cache->pdigits);

    bn_assign_zero(m2n, ndigits);
    bn_assign(m2n, m2, cache->qdigits);
    bn_add(out, qh, m2n, ndigits);
    if(bn_cmp(out, n, ndigits) >= 0) {
        bn_sub(out, out, n, ndigits);
    }

    // Clear potentially sensitive information
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
    mont_public_cache_t *cache;

    (void)n;
    (void)ndigits;

    cache = get_public_cache(pk);
    if(cache->edigits == 1 && cache->e[0] == 65537) {
        montgomery_mod_exp_65537(out, in, cache->n, cache->ndigits,
                                 cache->n0inv, cache->rr, cache->one);
    } else {
        montgomery_mod_exp_precomp(out, in, cache->e, cache->edigits, cache->n, cache->ndigits,
                                   cache->n0inv, cache->r, cache->rr);
    }

    return 0;
}

static void montgomery_mod_exp_precomp(bn_t *a, bn_t *b, bn_t *c, uint32_t cdigits,
                                       bn_t *n, uint32_t ndigits, bn_t n0inv, bn_t *r, bn_t *rr)
{
    bn_t base[BN_MAX_DIGITS], base_m[BN_MAX_DIGITS], one[BN_MAX_DIGITS];
    uint32_t bits;

    cdigits = bn_digits(c, cdigits);
    if(cdigits == 0) {
        BN_ASSIGN_DIGIT(a, 1, ndigits);
        return;
    }

    if(bn_cmp(b, n, ndigits) >= 0) {
        bn_mod(base, b, ndigits, n, ndigits);
    } else {
        bn_assign(base, b, ndigits);
    }

    montgomery_mul(base_m, base, rr, n, ndigits, n0inv);
    bits = highest_bit_index(c, cdigits) + 1;
    if(bits >= MONT_WINDOW_MIN_BITS) {
        montgomery_mod_exp_window(a, base_m, c, cdigits, n, ndigits, n0inv, r);
    } else {
        montgomery_mod_exp_binary(a, base_m, c, cdigits, n, ndigits, n0inv);
    }

    BN_ASSIGN_DIGIT(one, 1, ndigits);
    montgomery_mul(a, a, one, n, ndigits, n0inv);

    // Clear potentially sensitive information
    memset((uint8_t *)base, 0, sizeof(base));
    memset((uint8_t *)base_m, 0, sizeof(base_m));
    memset((uint8_t *)one, 0, sizeof(one));
}

static void montgomery_mod_exp_binary(bn_t *a, bn_t *base_m, bn_t *c, uint32_t cdigits,
                                      bn_t *n, uint32_t ndigits, bn_t n0inv)
{
    bn_t result[BN_MAX_DIGITS];
    int bit;

    bit = (int)highest_bit_index(c, cdigits);
    bn_assign(result, base_m, ndigits);
    for(bit--; bit>=0; bit--) {
        montgomery_mul(result, result, result, n, ndigits, n0inv);
        if(exponent_bit(c, (uint32_t)bit)) {
            montgomery_mul(result, result, base_m, n, ndigits, n0inv);
        }
    }

    bn_assign(a, result, ndigits);

    // Clear potentially sensitive information
    memset((uint8_t *)result, 0, sizeof(result));
}

static void montgomery_mod_exp_window(bn_t *a, bn_t *base_m, bn_t *c, uint32_t cdigits,
                                      bn_t *n, uint32_t ndigits, bn_t n0inv, bn_t *r)
{
    bn_t table[MONT_WINDOW_TABLE][BN_MAX_DIGITS], base2[BN_MAX_DIGITS], result[BN_MAX_DIGITS];
    int bit;

    bn_assign(table[0], base_m, ndigits);
    montgomery_mul(base2, base_m, base_m, n, ndigits, n0inv);
    for(uint32_t i=1; i<MONT_WINDOW_TABLE; i++) {
        montgomery_mul(table[i], table[i-1], base2, n, ndigits, n0inv);
    }

    bn_assign(result, r, ndigits);
    bit = (int)highest_bit_index(c, cdigits);
    while(bit >= 0) {
        if(!exponent_bit(c, (uint32_t)bit)) {
            montgomery_mul(result, result, result, n, ndigits, n0inv);
            bit--;
        } else {
            int low = bit - MONT_WINDOW_BITS + 1;
            uint32_t value = 0;

            if(low < 0) {
                low = 0;
            }
            while(!exponent_bit(c, (uint32_t)low)) {
                low++;
            }

            for(int j=bit; j>=low; j--) {
                value = (value << 1) | exponent_bit(c, (uint32_t)j);
                montgomery_mul(result, result, result, n, ndigits, n0inv);
            }
            montgomery_mul(result, result, table[value >> 1], n, ndigits, n0inv);
            bit = low - 1;
        }
    }

    bn_assign(a, result, ndigits);

    // Clear potentially sensitive information
    memset((uint8_t *)table, 0, sizeof(table));
    memset((uint8_t *)base2, 0, sizeof(base2));
    memset((uint8_t *)result, 0, sizeof(result));
}

static void montgomery_mod_exp_65537(bn_t *a, bn_t *b, bn_t *n, uint32_t ndigits,
                                     bn_t n0inv, bn_t *rr, bn_t *one)
{
    bn_t base[BN_MAX_DIGITS], base_m[BN_MAX_DIGITS], result[BN_MAX_DIGITS];

    if(bn_cmp(b, n, ndigits) >= 0) {
        bn_mod(base, b, ndigits, n, ndigits);
    } else {
        bn_assign(base, b, ndigits);
    }

    montgomery_mul(base_m, base, rr, n, ndigits, n0inv);
    bn_assign(result, base_m, ndigits);
    for(uint32_t i=0; i<16; i++) {
        montgomery_mul(result, result, result, n, ndigits, n0inv);
    }
    montgomery_mul(result, result, base_m, n, ndigits, n0inv);
    montgomery_mul(a, result, one, n, ndigits, n0inv);

    // Clear potentially sensitive information
    memset((uint8_t *)base, 0, sizeof(base));
    memset((uint8_t *)base_m, 0, sizeof(base_m));
    memset((uint8_t *)result, 0, sizeof(result));
}

static void montgomery_mod_mul_precomp(bn_t *a, bn_t *b, bn_t *c, bn_t *n,
                                       uint32_t digits, bn_t n0inv, bn_t *rr)
{
    bn_t bm[BN_MAX_DIGITS], cm[BN_MAX_DIGITS], tm[BN_MAX_DIGITS], one[BN_MAX_DIGITS];

    montgomery_mul(bm, b, rr, n, digits, n0inv);
    montgomery_mul(cm, c, rr, n, digits, n0inv);
    montgomery_mul(tm, bm, cm, n, digits, n0inv);
    BN_ASSIGN_DIGIT(one, 1, digits);
    montgomery_mul(a, tm, one, n, digits, n0inv);

    // Clear potentially sensitive information
    memset((uint8_t *)bm, 0, sizeof(bm));
    memset((uint8_t *)cm, 0, sizeof(cm));
    memset((uint8_t *)tm, 0, sizeof(tm));
    memset((uint8_t *)one, 0, sizeof(one));
}

static void montgomery_precompute_constants(bn_t *r, bn_t *rr, bn_t *n, uint32_t digits)
{
    bn_t tmp[2*BN_MAX_DIGITS];

    bn_assign_zero(tmp, 2*digits);
    tmp[digits] = 1;
    bn_mod(r, tmp, digits + 1, n, digits);

    mul_bn_mul(tmp, r, r, digits);
    bn_mod(rr, tmp, 2*digits, n, digits);

    // Clear potentially sensitive information
    memset((uint8_t *)tmp, 0, sizeof(tmp));
}

static mont_private_cache_t *get_private_cache(rsa_sk_t *sk)
{
    static mont_private_cache_t cache;

    if(cache.sk == sk) {
        return &cache;
    }

    memset((uint8_t *)&cache, 0, sizeof(cache));
    cache.sk = sk;

    bn_decode(cache.p, BN_MAX_DIGITS, sk->prime1, RSA_MAX_PRIME_LEN);
    bn_decode(cache.q, BN_MAX_DIGITS, sk->prime2, RSA_MAX_PRIME_LEN);
    bn_decode(cache.dp, BN_MAX_DIGITS, sk->prime_exponent1, RSA_MAX_PRIME_LEN);
    bn_decode(cache.dq, BN_MAX_DIGITS, sk->prime_exponent2, RSA_MAX_PRIME_LEN);
    bn_decode(cache.qinv, BN_MAX_DIGITS, sk->coefficient, RSA_MAX_PRIME_LEN);

    cache.pdigits = bn_digits(cache.p, BN_MAX_DIGITS);
    cache.qdigits = bn_digits(cache.q, BN_MAX_DIGITS);
    cache.dpdigits = bn_digits(cache.dp, cache.pdigits);
    cache.dqdigits = bn_digits(cache.dq, cache.qdigits);
    cache.p_n0inv = montgomery_n0inv(cache.p[0]);
    cache.q_n0inv = montgomery_n0inv(cache.q[0]);
    montgomery_precompute_constants(cache.rp, cache.rrp, cache.p, cache.pdigits);
    montgomery_precompute_constants(cache.rq, cache.rrq, cache.q, cache.qdigits);

    return &cache;
}

static mont_public_cache_t *get_public_cache(rsa_pk_t *pk)
{
    static mont_public_cache_t cache;

    if(cache.pk == pk) {
        return &cache;
    }

    memset((uint8_t *)&cache, 0, sizeof(cache));
    cache.pk = pk;

    bn_decode(cache.n, BN_MAX_DIGITS, pk->modulus, RSA_MAX_MODULUS_LEN);
    bn_decode(cache.e, BN_MAX_DIGITS, pk->exponent, RSA_MAX_MODULUS_LEN);
    cache.ndigits = bn_digits(cache.n, BN_MAX_DIGITS);
    cache.edigits = bn_digits(cache.e, BN_MAX_DIGITS);
    cache.n0inv = montgomery_n0inv(cache.n[0]);
    montgomery_precompute_constants(cache.r, cache.rr, cache.n, cache.ndigits);
    BN_ASSIGN_DIGIT(cache.one, 1, cache.ndigits);

    return &cache;
}

static void classic_mod_exp(bn_t *a, bn_t *b, bn_t *c, uint32_t cdigits, bn_t *n, uint32_t ndigits)
{
    bn_t result[BN_MAX_DIGITS], base[BN_MAX_DIGITS];
    int bit;

    cdigits = bn_digits(c, cdigits);
    if(cdigits == 0) {
        BN_ASSIGN_DIGIT(a, 1, ndigits);
        return;
    }

    if(bn_cmp(b, n, ndigits) >= 0) {
        bn_mod(base, b, ndigits, n, ndigits);
    } else {
        bn_assign(base, b, ndigits);
    }

    bit = (int)highest_bit_index(c, cdigits);
    bn_assign(result, base, ndigits);
    for(bit--; bit>=0; bit--) {
        bn_mod_mul(result, result, result, n, ndigits);
        if(exponent_bit(c, (uint32_t)bit)) {
            bn_mod_mul(result, result, base, n, ndigits);
        }
    }

    bn_assign(a, result, ndigits);

    // Clear potentially sensitive information
    memset((uint8_t *)result, 0, sizeof(result));
    memset((uint8_t *)base, 0, sizeof(base));
}

static void montgomery_mul(bn_t *a, bn_t *b, bn_t *c, bn_t *n, uint32_t digits, bn_t n0inv)
{
    bn_t t[BN_MAX_DIGITS+1], m, borrow;
    dbn_t acc, carry;
    uint32_t i, j;

    bn_assign_zero(t, digits + 1);
    for(i=0; i<digits; i++) {
        bn_t bi = b[i];

        m = (bn_t)(((dbn_t)t[0] + (dbn_t)b[i] * c[0]) * n0inv);
        carry = 0;

        for(j=0; j<digits; j++) {
            acc = (dbn_t)bi * c[j] + t[j] + carry;
            carry = acc >> BN_DIGIT_BITS;
            acc = (dbn_t)m * n[j] + (bn_t)acc;
            carry += acc >> BN_DIGIT_BITS;
            if(j > 0) {
                t[j-1] = (bn_t)acc;
            }
        }

        acc = (dbn_t)t[digits] + carry;
        t[digits-1] = (bn_t)acc;
        t[digits] = (bn_t)(acc >> BN_DIGIT_BITS);
    }

    bn_assign(a, t, digits);
    if(t[digits] || bn_cmp(a, n, digits) >= 0) {
        borrow = bn_sub(a, a, n, digits);
        (void)borrow;
    }

    // Clear potentially sensitive information
    memset((uint8_t *)t, 0, sizeof(t));
}

static bn_t montgomery_n0inv(bn_t n0)
{
    bn_t x;

    x = 1;
    for(uint32_t i=0; i<5; i++) {
        x *= 2 - n0 * x;
    }

    return (bn_t)(0 - x);
}

static uint32_t exponent_bit(bn_t *a, uint32_t bit)
{
    return (a[bit / BN_DIGIT_BITS] >> (bit % BN_DIGIT_BITS)) & 1U;
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
