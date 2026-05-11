#include <pthread.h>
#include <string.h>
#include "modmul.h"
#include "../mul/mul.h"

#define MODMUL_MT_MIN_DIGITS 128

typedef struct {
    bn_t *out;
    bn_t *in;
    bn_t *rr;
    bn_t *n;
    uint32_t digits;
    bn_t n0inv;
} montgomery_worker_t;

static bn_t montgomery_n0inv(bn_t n0);
static void base_mod_mul(bn_t *a, bn_t *b, bn_t *c, bn_t *d, uint32_t digits);
static void montgomery_reduce(bn_t *a, bn_t *t, bn_t *n, uint32_t digits, bn_t n0inv);
static void montgomery_mul(bn_t *a, bn_t *b, bn_t *c, bn_t *n, uint32_t digits, bn_t n0inv);
static void reduce_operand(bn_t *a, bn_t *b, bn_t *n, uint32_t digits);
static void *to_montgomery_thread(void *arg);

void modmul_bn_mod_mul(bn_t *a, bn_t *b, bn_t *c, bn_t *d, uint32_t digits)
{
    bn_t bb[BN_MAX_DIGITS], cc[BN_MAX_DIGITS], r[BN_MAX_DIGITS], rr[BN_MAX_DIGITS];
    bn_t tmp[2*BN_MAX_DIGITS], bm[BN_MAX_DIGITS], cm[BN_MAX_DIGITS], tm[BN_MAX_DIGITS];
    bn_t n0inv;
    montgomery_worker_t worker;
    pthread_t thread;
    int threaded;

    if(digits < MODMUL_MT_MIN_DIGITS || d[0] == 0 || ((d[0] & 1) == 0)) {
        base_mod_mul(a, b, c, d, digits);
        return;
    }

    reduce_operand(bb, b, d, digits);
    reduce_operand(cc, c, d, digits);

    bn_assign_zero(tmp, 2*digits);
    tmp[digits] = 1;
    bn_mod(r, tmp, digits+1, d, digits);

    mul_bn_mul(tmp, r, r, digits);
    bn_mod(rr, tmp, 2*digits, d, digits);

    n0inv = montgomery_n0inv(d[0]);

    worker.out = cm;
    worker.in = cc;
    worker.rr = rr;
    worker.n = d;
    worker.digits = digits;
    worker.n0inv = n0inv;

    threaded = (pthread_create(&thread, NULL, to_montgomery_thread, &worker) == 0);
    montgomery_mul(bm, bb, rr, d, digits, n0inv);
    if(threaded) {
        pthread_join(thread, NULL);
    } else {
        montgomery_mul(cm, cc, rr, d, digits, n0inv);
    }

    montgomery_mul(tm, bm, cm, d, digits, n0inv);

    BN_ASSIGN_DIGIT(r, 1, digits);
    montgomery_mul(a, tm, r, d, digits, n0inv);

    // Clear potentially sensitive information
    memset((uint8_t *)bb, 0, sizeof(bb));
    memset((uint8_t *)cc, 0, sizeof(cc));
    memset((uint8_t *)r, 0, sizeof(r));
    memset((uint8_t *)rr, 0, sizeof(rr));
    memset((uint8_t *)tmp, 0, sizeof(tmp));
    memset((uint8_t *)bm, 0, sizeof(bm));
    memset((uint8_t *)cm, 0, sizeof(cm));
    memset((uint8_t *)tm, 0, sizeof(tm));
    memset((uint8_t *)&worker, 0, sizeof(worker));
}

static void *to_montgomery_thread(void *arg)
{
    montgomery_worker_t *worker = (montgomery_worker_t *)arg;

    montgomery_mul(worker->out, worker->in, worker->rr, worker->n, worker->digits, worker->n0inv);
    return NULL;
}

static void base_mod_mul(bn_t *a, bn_t *b, bn_t *c, bn_t *d, uint32_t digits)
{
    bn_t t[2*BN_MAX_DIGITS];

    mul_bn_mul(t, b, c, digits);
    bn_mod(a, t, 2*digits, d, digits);

    // Clear potentially sensitive information
    memset((uint8_t *)t, 0, sizeof(t));
}

static void reduce_operand(bn_t *a, bn_t *b, bn_t *n, uint32_t digits)
{
    if(bn_cmp(b, n, digits) >= 0) {
        bn_mod(a, b, digits, n, digits);
    } else {
        bn_assign(a, b, digits);
    }
}

static void montgomery_mul(bn_t *a, bn_t *b, bn_t *c, bn_t *n, uint32_t digits, bn_t n0inv)
{
    bn_t t[2*BN_MAX_DIGITS+1];

    bn_assign_zero(t, 2*digits+1);
    mul_bn_mul(t, b, c, digits);
    montgomery_reduce(a, t, n, digits, n0inv);

    // Clear potentially sensitive information
    memset((uint8_t *)t, 0, sizeof(t));
}

static void montgomery_reduce(bn_t *a, bn_t *t, bn_t *n, uint32_t digits, bn_t n0inv)
{
    dbn_t uv;
    uint32_t i, j, k;

    for(i=0; i<digits; i++) {
        bn_t m = t[i] * n0inv;
        bn_t carry = 0;

        for(j=0; j<digits; j++) {
            uv = (dbn_t)t[i+j] + (dbn_t)m * n[j] + carry;
            t[i+j] = (bn_t)uv;
            carry = (bn_t)(uv >> BN_DIGIT_BITS);
        }

        k = i + digits;
        while(carry != 0) {
            uv = (dbn_t)t[k] + carry;
            t[k] = (bn_t)uv;
            carry = (bn_t)(uv >> BN_DIGIT_BITS);
            k++;
        }
    }

    bn_assign(a, t+digits, digits);
    if(t[2*digits] != 0 || bn_cmp(a, n, digits) >= 0) {
        bn_sub(a, a, n, digits);
    }
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
