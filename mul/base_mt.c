#include <pthread.h>
#include <string.h>
#include "mul.h"

#define MUL_MT_MIN_DIGITS 128

typedef struct {
    bn_t *b;
    bn_t *c;
    bn_t local[2*BN_MAX_DIGITS];
    uint32_t b_start;
    uint32_t b_end;
    uint32_t cdigits;
    uint32_t digits;
} mul_worker_t;

static bn_t bn_add_digit_mul(bn_t *a, bn_t *b, bn_t c, bn_t *d, uint32_t digits);
static void add_product(bn_t *a, bn_t *b, uint32_t digits);
static void base_mul(bn_t *a, bn_t *b, bn_t *c, uint32_t digits);
static void mul_range(mul_worker_t *worker);
static void *mul_range_thread(void *arg);

void mul_bn_mul(bn_t *a, bn_t *b, bn_t *c, uint32_t digits)
{
    mul_worker_t low, high;
    pthread_t thread;
    uint32_t bdigits, split;
    int threaded;

    if(digits < MUL_MT_MIN_DIGITS) {
        base_mul(a, b, c, digits);
        return;
    }

    bdigits = bn_digits(b, digits);
    split = (bdigits + 1) / 2;

    low.b = b;
    low.c = c;
    low.b_start = 0;
    low.b_end = split;
    low.cdigits = bn_digits(c, digits);
    low.digits = digits;

    high.b = b;
    high.c = c;
    high.b_start = split;
    high.b_end = bdigits;
    high.cdigits = low.cdigits;
    high.digits = digits;

    bn_assign_zero(low.local, 2*digits);
    bn_assign_zero(high.local, 2*digits);

    threaded = (pthread_create(&thread, NULL, mul_range_thread, &high) == 0);
    mul_range(&low);
    if(threaded) {
        pthread_join(thread, NULL);
    } else {
        mul_range(&high);
    }

    bn_assign(a, low.local, 2*digits);
    add_product(a, high.local, 2*digits);

    // Clear potentially sensitive information
    memset((uint8_t *)&low, 0, sizeof(low));
    memset((uint8_t *)&high, 0, sizeof(high));
}

static void base_mul(bn_t *a, bn_t *b, bn_t *c, uint32_t digits)
{
    bn_t t[2*BN_MAX_DIGITS];
    uint32_t bdigits, cdigits, i;

    bn_assign_zero(t, 2*digits);
    bdigits = bn_digits(b, digits);
    cdigits = bn_digits(c, digits);

    for(i=0; i<bdigits; i++) {
        t[i+cdigits] += bn_add_digit_mul(&t[i], &t[i], b[i], c, cdigits);
    }

    bn_assign(a, t, 2*digits);

    // Clear potentially sensitive information
    memset((uint8_t *)t, 0, sizeof(t));
}

static void *mul_range_thread(void *arg)
{
    mul_range((mul_worker_t *)arg);
    return NULL;
}

static void mul_range(mul_worker_t *worker)
{
    uint32_t i;

    for(i=worker->b_start; i<worker->b_end; i++) {
        worker->local[i+worker->cdigits] += bn_add_digit_mul(&worker->local[i],
                                                             &worker->local[i],
                                                             worker->b[i],
                                                             worker->c,
                                                             worker->cdigits);
    }
}

static void add_product(bn_t *a, bn_t *b, uint32_t digits)
{
    dbn_t sum;
    bn_t carry;
    uint32_t i;

    carry = 0;
    for(i=0; i<digits; i++) {
        sum = (dbn_t)a[i] + b[i] + carry;
        a[i] = (bn_t)sum;
        carry = (bn_t)(sum >> BN_DIGIT_BITS);
    }
}

static bn_t bn_add_digit_mul(bn_t *a, bn_t *b, bn_t c, bn_t *d, uint32_t digits)
{
    dbn_t result;
    bn_t carry, rh, rl;
    uint32_t i;

    if(c == 0)
        return 0;

    carry = 0;
    for(i=0; i<digits; i++) {
        result = (dbn_t)c * d[i];
        rl = result & BN_MAX_DIGIT;
        rh = (result >> BN_DIGIT_BITS) & BN_MAX_DIGIT;
        if((a[i] = b[i] + carry) < carry) {
            carry = 1;
        } else {
            carry = 0;
        }
        if((a[i] += rl) < rl) {
            carry++;
        }
        carry += rh;
    }

    return carry;
}
