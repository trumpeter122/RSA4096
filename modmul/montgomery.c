#include <pthread.h>
#include <string.h>
#include "modmul.h"
#include "../mul/mul.h"

#define MONT_CACHE_SIZE 4

typedef struct {
    int valid;
    uint32_t digits;
    bn_t n[BN_MAX_DIGITS];
    bn_t r2[BN_MAX_DIGITS];
    bn_t n0inv;
} mont_cache_t;

static mont_cache_t g_mont_cache[MONT_CACHE_SIZE];
static uint32_t g_mont_cache_next = 0;
static pthread_mutex_t g_mont_cache_lock = PTHREAD_MUTEX_INITIALIZER;

static bn_t montgomery_n0inv(bn_t n0);
static void montgomery_r2(bn_t *r2, bn_t *n, uint32_t digits);
static void montgomery_get_constants(bn_t *r2, bn_t *n0inv, bn_t *n, uint32_t digits);
static void montgomery_mul(bn_t *a, bn_t *b, bn_t *c, bn_t *n, uint32_t digits, bn_t n0inv);
static void reduce_operand(bn_t *a, bn_t *b, bn_t *n, uint32_t digits);

void modmul_bn_mod_mul(bn_t *a, bn_t *b, bn_t *c, bn_t *d, uint32_t digits)
{
    bn_t bb[BN_MAX_DIGITS], cc[BN_MAX_DIGITS], one[BN_MAX_DIGITS];
    bn_t r2[BN_MAX_DIGITS], bm[BN_MAX_DIGITS], cm[BN_MAX_DIGITS], tm[BN_MAX_DIGITS];
    bn_t n0inv;

    reduce_operand(bb, b, d, digits);
    reduce_operand(cc, c, d, digits);

    montgomery_get_constants(r2, &n0inv, d, digits);
    montgomery_mul(bm, bb, r2, d, digits, n0inv);
    montgomery_mul(cm, cc, r2, d, digits, n0inv);
    montgomery_mul(tm, bm, cm, d, digits, n0inv);

    BN_ASSIGN_DIGIT(one, 1, digits);
    montgomery_mul(a, tm, one, d, digits, n0inv);

    // Clear potentially sensitive information
    memset((uint8_t *)bb, 0, sizeof(bb));
    memset((uint8_t *)cc, 0, sizeof(cc));
    memset((uint8_t *)one, 0, sizeof(one));
    memset((uint8_t *)r2, 0, sizeof(r2));
    memset((uint8_t *)bm, 0, sizeof(bm));
    memset((uint8_t *)cm, 0, sizeof(cm));
    memset((uint8_t *)tm, 0, sizeof(tm));
}

static void reduce_operand(bn_t *a, bn_t *b, bn_t *n, uint32_t digits)
{
    if(bn_cmp(b, n, digits) >= 0) {
        bn_mod(a, b, digits, n, digits);
    } else {
        bn_assign(a, b, digits);
    }
}

static void montgomery_r2(bn_t *r2, bn_t *n, uint32_t digits)
{
    bn_t dividend[2*BN_MAX_DIGITS+1], quotient[2*BN_MAX_DIGITS+1];

    bn_assign_zero(dividend, 2*BN_MAX_DIGITS+1);
    dividend[2*digits] = 1;
    bn_div(quotient, r2, dividend, 2*digits + 1, n, digits);

    // Clear potentially sensitive information
    memset((uint8_t *)dividend, 0, sizeof(dividend));
    memset((uint8_t *)quotient, 0, sizeof(quotient));
}

static void montgomery_get_constants(bn_t *r2, bn_t *n0inv, bn_t *n, uint32_t digits)
{
    mont_cache_t *entry;
    uint32_t i;

    pthread_mutex_lock(&g_mont_cache_lock);
    for(i=0; i<MONT_CACHE_SIZE; i++) {
        entry = &g_mont_cache[i];
        if(entry->valid && entry->digits == digits && bn_cmp(entry->n, n, digits) == 0) {
            *n0inv = entry->n0inv;
            bn_assign(r2, entry->r2, digits);
            pthread_mutex_unlock(&g_mont_cache_lock);
            return;
        }
    }

    *n0inv = montgomery_n0inv(n[0]);
    montgomery_r2(r2, n, digits);

    entry = &g_mont_cache[g_mont_cache_next++ % MONT_CACHE_SIZE];
    entry->valid = 1;
    entry->digits = digits;
    entry->n0inv = *n0inv;
    bn_assign_zero(entry->n, BN_MAX_DIGITS);
    bn_assign_zero(entry->r2, BN_MAX_DIGITS);
    bn_assign(entry->n, n, digits);
    bn_assign(entry->r2, r2, digits);
    pthread_mutex_unlock(&g_mont_cache_lock);
}

static void montgomery_mul(bn_t *a, bn_t *b, bn_t *c, bn_t *n, uint32_t digits, bn_t n0inv)
{
    bn_t t[BN_MAX_DIGITS+1], m, borrow;
    unsigned __int128 acc, carry;
    uint32_t i, j;

    bn_assign_zero(t, digits + 1);
    for(i=0; i<digits; i++) {
        m = (bn_t)(((dbn_t)t[0] + (dbn_t)b[i] * c[0]) * n0inv);
        carry = 0;

        for(j=0; j<digits; j++) {
            acc = (unsigned __int128)t[j] + (unsigned __int128)b[i] * c[j] +
                  (unsigned __int128)m * n[j] + carry;
            if(j > 0) {
                t[j-1] = (bn_t)acc;
            }
            carry = acc >> BN_DIGIT_BITS;
        }

        acc = (unsigned __int128)t[digits] + carry;
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
    uint32_t i;

    x = 1;
    for(i=0; i<5; i++) {
        x *= 2 - n0 * x;
    }

    return (bn_t)(0 - x);
}
