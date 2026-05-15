#include <string.h>
#include "mul.h"

#define NTT_MAX_LEN 512
#define NTT_PRIMES 3

typedef __uint128_t u128_t;

static const uint32_t ntt_mods[NTT_PRIMES] = {998244353U, 1004535809U, 469762049U};
static const uint32_t ntt_roots[NTT_PRIMES] = {3U, 3U, 3U};

static uint32_t mod_pow(uint32_t a, uint32_t e, uint32_t mod);
static uint32_t mod_inv(uint32_t a, uint32_t mod);
static void ntt(uint32_t *a, uint32_t n, int invert, uint32_t mod, uint32_t root);
static uint32_t next_power_of_two(uint32_t n);
static u128_t crt3(uint32_t r0, uint32_t r1, uint32_t r2);

void mul_bn_mul(bn_t *a, bn_t *b, bn_t *c, uint32_t digits)
{
    uint32_t len;
    uint32_t residues[NTT_PRIMES][2*BN_MAX_DIGITS];

    len = next_power_of_two(2 * digits);

    for(uint32_t p=0; p<NTT_PRIMES; p++) {
        uint32_t fa[NTT_MAX_LEN], fb[NTT_MAX_LEN];
        uint32_t mod = ntt_mods[p];

        memset((uint8_t *)fa, 0, sizeof(fa));
        memset((uint8_t *)fb, 0, sizeof(fb));
        for(uint32_t i=0; i<digits; i++) {
            fa[i] = b[i] % mod;
            fb[i] = c[i] % mod;
        }

        ntt(fa, len, 0, mod, ntt_roots[p]);
        ntt(fb, len, 0, mod, ntt_roots[p]);
        for(uint32_t i=0; i<len; i++) {
            fa[i] = (uint32_t)(((dbn_t)fa[i] * fb[i]) % mod);
        }
        ntt(fa, len, 1, mod, ntt_roots[p]);

        for(uint32_t i=0; i<2*digits; i++) {
            residues[p][i] = fa[i];
        }

        // Clear potentially sensitive information
        memset((uint8_t *)fa, 0, sizeof(fa));
        memset((uint8_t *)fb, 0, sizeof(fb));
    }

    {
        u128_t carry = 0;
        for(uint32_t i=0; i<2*digits; i++) {
            u128_t coeff = crt3(residues[0][i], residues[1][i], residues[2][i]) + carry;
            a[i] = (bn_t)coeff;
            carry = coeff >> BN_DIGIT_BITS;
        }
    }

    // Clear potentially sensitive information
    memset((uint8_t *)residues, 0, sizeof(residues));
}

static u128_t crt3(uint32_t r0, uint32_t r1, uint32_t r2)
{
    static uint32_t inv_m0_m1 = 0;
    static uint32_t inv_m0m1_m2 = 0;
    const uint64_t m0 = 998244353ULL;
    const uint64_t m1 = 1004535809ULL;
    const uint64_t m2 = 469762049ULL;
    uint64_t t1, t2, x_mod_m2;

    if(inv_m0_m1 == 0) {
        inv_m0_m1 = mod_inv((uint32_t)(m0 % m1), (uint32_t)m1);
        inv_m0m1_m2 = mod_inv((uint32_t)((m0 % m2) * (m1 % m2) % m2), (uint32_t)m2);
    }

    t1 = ((r1 + m1 - (uint32_t)(r0 % m1)) * (uint64_t)inv_m0_m1) % m1;
    x_mod_m2 = (r0 + (m0 % m2) * (t1 % m2)) % m2;
    t2 = ((r2 + m2 - x_mod_m2) * (uint64_t)inv_m0m1_m2) % m2;

    return (u128_t)r0 + (u128_t)m0 * t1 + (u128_t)m0 * m1 * t2;
}

static void ntt(uint32_t *a, uint32_t n, int invert, uint32_t mod, uint32_t root)
{
    for(uint32_t i=1, j=0; i<n; i++) {
        uint32_t bit = n >> 1;
        for(; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if(i < j) {
            uint32_t tmp = a[i];
            a[i] = a[j];
            a[j] = tmp;
        }
    }

    for(uint32_t len=2; len<=n; len<<=1) {
        uint32_t wlen = mod_pow(root, (mod - 1) / len, mod);
        if(invert) {
            wlen = mod_inv(wlen, mod);
        }

        for(uint32_t i=0; i<n; i+=len) {
            uint32_t w = 1;
            for(uint32_t j=0; j<len/2; j++) {
                uint32_t u = a[i+j];
                uint32_t v = (uint32_t)(((dbn_t)a[i+j+len/2] * w) % mod);
                a[i+j] = u + v < mod ? u + v : u + v - mod;
                a[i+j+len/2] = u >= v ? u - v : u + mod - v;
                w = (uint32_t)(((dbn_t)w * wlen) % mod);
            }
        }
    }

    if(invert) {
        uint32_t n_inv = mod_inv(n, mod);
        for(uint32_t i=0; i<n; i++) {
            a[i] = (uint32_t)(((dbn_t)a[i] * n_inv) % mod);
        }
    }
}

static uint32_t mod_pow(uint32_t a, uint32_t e, uint32_t mod)
{
    dbn_t result = 1;
    dbn_t base = a;

    while(e != 0) {
        if(e & 1U) {
            result = (result * base) % mod;
        }
        base = (base * base) % mod;
        e >>= 1;
    }

    return (uint32_t)result;
}

static uint32_t mod_inv(uint32_t a, uint32_t mod)
{
    return mod_pow(a, mod - 2, mod);
}

static uint32_t next_power_of_two(uint32_t n)
{
    uint32_t value = 1;
    while(value < n) {
        value <<= 1;
    }
    return value;
}
