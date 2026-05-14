#include <pthread.h>
#include <string.h>
#include "modmul.h"
#include "../mul/mul.h"

#define BARRETT_CACHE_SIZE 4

typedef struct {
    int valid;
    uint32_t digits;
    bn_t n[BN_MAX_DIGITS];
    bn_t mu[BN_MAX_DIGITS+2];
} barrett_cache_t;

static barrett_cache_t g_barrett_cache[BARRETT_CACHE_SIZE];
static uint32_t g_barrett_cache_next = 0;
static pthread_mutex_t g_barrett_cache_lock = PTHREAD_MUTEX_INITIALIZER;

static void barrett_get_mu(bn_t *mu, bn_t *n, uint32_t digits);
static void barrett_reduce(bn_t *a, bn_t *x, bn_t *n, uint32_t digits);
static void base_mod_mul(bn_t *a, bn_t *b, bn_t *c, bn_t *d, uint32_t digits);

void modmul_bn_mod_mul(bn_t *a, bn_t *b, bn_t *c, bn_t *d, uint32_t digits)
{
    bn_t t[2*BN_MAX_DIGITS];

    if(digits == 0 || d[0] == 0) {
        base_mod_mul(a, b, c, d, digits);
        return;
    }

    mul_bn_mul(t, b, c, digits);
    barrett_reduce(a, t, d, digits);

    // Clear potentially sensitive information
    memset((uint8_t *)t, 0, sizeof(t));
}

static void base_mod_mul(bn_t *a, bn_t *b, bn_t *c, bn_t *d, uint32_t digits)
{
    bn_t t[2*BN_MAX_DIGITS];

    mul_bn_mul(t, b, c, digits);
    bn_mod(a, t, 2*digits, d, digits);

    // Clear potentially sensitive information
    memset((uint8_t *)t, 0, sizeof(t));
}

static void barrett_reduce(bn_t *a, bn_t *x, bn_t *n, uint32_t digits)
{
    bn_t mu[BN_MAX_DIGITS+2], q1[BN_MAX_DIGITS+2], q2[2*BN_MAX_DIGITS+4];
    bn_t q3[BN_MAX_DIGITS+2], r[BN_MAX_DIGITS+2], r2_full[2*BN_MAX_DIGITS+4];
    bn_t n_ext[BN_MAX_DIGITS+2];
    uint32_t r_digits = digits + 2;
    uint32_t q1_digits = digits + 1;
    uint32_t q3_digits = digits + 1;

    barrett_get_mu(mu, n, digits);

    bn_assign_zero(q1, BN_MAX_DIGITS+2);
    bn_assign(q1, x + digits - 1, q1_digits);
    mul_bn_mul(q2, q1, mu, q1_digits);

    bn_assign_zero(q3, BN_MAX_DIGITS+2);
    bn_assign(q3, q2 + digits + 1, q3_digits);

    bn_assign_zero(r, BN_MAX_DIGITS+2);
    bn_assign(r, x, digits + 1);

    bn_assign_zero(n_ext, BN_MAX_DIGITS+2);
    bn_assign(n_ext, n, digits);

    bn_assign_zero(r2_full, 2*(digits+1));
    mul_bn_mul(r2_full, q3, n_ext, digits + 1);
    r2_full[digits+1] = 0;

    if(bn_cmp(r, r2_full, r_digits) >= 0) {
        bn_sub(r, r, r2_full, r_digits);
    } else {
        bn_t base_ext[BN_MAX_DIGITS+2];
        bn_assign_zero(base_ext, BN_MAX_DIGITS+2);
        base_ext[digits+1] = 1;
        bn_sub(r, r2_full, r, r_digits);
        bn_sub(r, base_ext, r, r_digits);
        memset((uint8_t *)base_ext, 0, sizeof(base_ext));
    }

    uint32_t corrections = 0;
    while(bn_cmp(r, n_ext, r_digits) >= 0 && corrections < 3) {
        bn_sub(r, r, n_ext, r_digits);
        corrections++;
    }

    if(bn_cmp(r, n_ext, r_digits) >= 0) {
        bn_mod(a, x, 2*digits, n, digits);
        memset((uint8_t *)mu, 0, sizeof(mu));
        memset((uint8_t *)q1, 0, sizeof(q1));
        memset((uint8_t *)q2, 0, sizeof(q2));
        memset((uint8_t *)q3, 0, sizeof(q3));
        memset((uint8_t *)r, 0, sizeof(r));
        memset((uint8_t *)r2_full, 0, sizeof(r2_full));
        memset((uint8_t *)n_ext, 0, sizeof(n_ext));
        return;
    }

    bn_assign(a, r, digits);

    // Clear potentially sensitive information
    memset((uint8_t *)mu, 0, sizeof(mu));
    memset((uint8_t *)q1, 0, sizeof(q1));
    memset((uint8_t *)q2, 0, sizeof(q2));
    memset((uint8_t *)q3, 0, sizeof(q3));
    memset((uint8_t *)r, 0, sizeof(r));
    memset((uint8_t *)r2_full, 0, sizeof(r2_full));
    memset((uint8_t *)n_ext, 0, sizeof(n_ext));
}

static void barrett_get_mu(bn_t *mu, bn_t *n, uint32_t digits)
{
    barrett_cache_t *entry;
    bn_t dividend[2*BN_MAX_DIGITS+3], quotient[2*BN_MAX_DIGITS+3], rem[BN_MAX_DIGITS+2];

    pthread_mutex_lock(&g_barrett_cache_lock);
    for(uint32_t i=0; i<BARRETT_CACHE_SIZE; i++) {
        entry = &g_barrett_cache[i];
        if(entry->valid && entry->digits == digits && bn_cmp(entry->n, n, digits) == 0) {
            bn_assign(mu, entry->mu, digits + 1);
            pthread_mutex_unlock(&g_barrett_cache_lock);
            return;
        }
    }

    bn_assign_zero(dividend, 2*BN_MAX_DIGITS+3);
    bn_assign_zero(quotient, 2*BN_MAX_DIGITS+3);
    bn_assign_zero(rem, BN_MAX_DIGITS+2);
    dividend[2*digits] = 1;
    bn_div(quotient, rem, dividend, 2*digits + 1, n, digits);
    bn_assign(mu, quotient, digits + 1);

    entry = &g_barrett_cache[g_barrett_cache_next++ % BARRETT_CACHE_SIZE];
    entry->valid = 1;
    entry->digits = digits;
    bn_assign_zero(entry->n, BN_MAX_DIGITS);
    bn_assign_zero(entry->mu, BN_MAX_DIGITS+2);
    bn_assign(entry->n, n, digits);
    bn_assign(entry->mu, mu, digits + 1);

    pthread_mutex_unlock(&g_barrett_cache_lock);

    // Clear potentially sensitive information
    memset((uint8_t *)dividend, 0, sizeof(dividend));
    memset((uint8_t *)quotient, 0, sizeof(quotient));
    memset((uint8_t *)rem, 0, sizeof(rem));
}
