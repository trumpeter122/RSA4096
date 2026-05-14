#include <string.h>
#include "mul.h"

/*
 * The coefficient-space Toom-3 implementation below is kept as an experimental
 * implementation scaffold. RSA4096 operands are small enough that the safe
 * fallback is preferable until a full signed-magnitude interpolation path is
 * added.
 */
#define TOOM_MIN_DIGITS (BN_MAX_DIGITS + 1)

typedef __int128 i128_t;

static void school_mul(bn_t *a, bn_t *b, bn_t *c, uint32_t digits);
static void toom3_mul(bn_t *a, bn_t *b, bn_t *c, uint32_t digits);
static void eval_chunk(int64_t *out, bn_t *x, uint32_t start, uint32_t chunk,
                       uint32_t digits, int point);
static void mul_eval(i128_t *out, int64_t *x, int64_t *y, uint32_t len);
static void add_shifted(i128_t *out, uint32_t out_len, i128_t *x, uint32_t len,
                        uint32_t shift);
static void normalize_product(bn_t *a, i128_t *x, uint32_t digits);

void mul_bn_mul(bn_t *a, bn_t *b, bn_t *c, uint32_t digits)
{
    if(digits < TOOM_MIN_DIGITS) {
        school_mul(a, b, c, digits);
        return;
    }

    toom3_mul(a, b, c, digits);
}

static void toom3_mul(bn_t *a, bn_t *b, bn_t *c, uint32_t digits)
{
    uint32_t chunk = (digits + 2) / 3;
    uint32_t vlen = 2 * chunk;
    uint32_t out_len = 2 * digits;

    int64_t b0[BN_MAX_DIGITS], b1[BN_MAX_DIGITS], bm1[BN_MAX_DIGITS], b2[BN_MAX_DIGITS], binf[BN_MAX_DIGITS];
    int64_t c0[BN_MAX_DIGITS], c1[BN_MAX_DIGITS], cm1[BN_MAX_DIGITS], c2[BN_MAX_DIGITS], cinf[BN_MAX_DIGITS];
    i128_t w0[2*BN_MAX_DIGITS], w1[2*BN_MAX_DIGITS], wm1[2*BN_MAX_DIGITS], w2[2*BN_MAX_DIGITS], winf[2*BN_MAX_DIGITS];
    i128_t r0[2*BN_MAX_DIGITS], r1[2*BN_MAX_DIGITS], r2[2*BN_MAX_DIGITS], r3[2*BN_MAX_DIGITS], r4[2*BN_MAX_DIGITS];
    i128_t acc[2*BN_MAX_DIGITS+6];

    memset((uint8_t *)w0, 0, sizeof(w0));
    memset((uint8_t *)w1, 0, sizeof(w1));
    memset((uint8_t *)wm1, 0, sizeof(wm1));
    memset((uint8_t *)w2, 0, sizeof(w2));
    memset((uint8_t *)winf, 0, sizeof(winf));
    memset((uint8_t *)acc, 0, sizeof(acc));

    eval_chunk(b0, b, 0, chunk, digits, 0);
    eval_chunk(c0, c, 0, chunk, digits, 0);
    eval_chunk(b1, b, 0, chunk, digits, 1);
    eval_chunk(c1, c, 0, chunk, digits, 1);
    eval_chunk(bm1, b, 0, chunk, digits, -1);
    eval_chunk(cm1, c, 0, chunk, digits, -1);
    eval_chunk(b2, b, 0, chunk, digits, 2);
    eval_chunk(c2, c, 0, chunk, digits, 2);
    eval_chunk(binf, b, 2*chunk, chunk, digits, 0);
    eval_chunk(cinf, c, 2*chunk, chunk, digits, 0);

    mul_eval(w0, b0, c0, chunk);
    mul_eval(w1, b1, c1, chunk);
    mul_eval(wm1, bm1, cm1, chunk);
    mul_eval(w2, b2, c2, chunk);
    mul_eval(winf, binf, cinf, chunk);

    for(uint32_t i=0; i<vlen; i++) {
        r0[i] = w0[i];
        r4[i] = winf[i];
        r3[i] = (w2[i] - w1[i]) / 3;
        r1[i] = (w1[i] - wm1[i]) / 2;
        r2[i] = wm1[i] - w0[i];
        r3[i] = (r2[i] - r3[i]) / 2 + 2 * winf[i];
        r2[i] = r2[i] + r1[i] - winf[i];
        r1[i] = r1[i] - r3[i];
    }

    add_shifted(acc, out_len, r0, vlen, 0);
    add_shifted(acc, out_len, r1, vlen, chunk);
    add_shifted(acc, out_len, r2, vlen, 2*chunk);
    add_shifted(acc, out_len, r3, vlen, 3*chunk);
    add_shifted(acc, out_len, r4, vlen, 4*chunk);
    normalize_product(a, acc, digits);

    // Clear potentially sensitive information
    memset((uint8_t *)b0, 0, sizeof(b0));
    memset((uint8_t *)b1, 0, sizeof(b1));
    memset((uint8_t *)bm1, 0, sizeof(bm1));
    memset((uint8_t *)b2, 0, sizeof(b2));
    memset((uint8_t *)binf, 0, sizeof(binf));
    memset((uint8_t *)c0, 0, sizeof(c0));
    memset((uint8_t *)c1, 0, sizeof(c1));
    memset((uint8_t *)cm1, 0, sizeof(cm1));
    memset((uint8_t *)c2, 0, sizeof(c2));
    memset((uint8_t *)cinf, 0, sizeof(cinf));
    memset((uint8_t *)w0, 0, sizeof(w0));
    memset((uint8_t *)w1, 0, sizeof(w1));
    memset((uint8_t *)wm1, 0, sizeof(wm1));
    memset((uint8_t *)w2, 0, sizeof(w2));
    memset((uint8_t *)winf, 0, sizeof(winf));
    memset((uint8_t *)r0, 0, sizeof(r0));
    memset((uint8_t *)r1, 0, sizeof(r1));
    memset((uint8_t *)r2, 0, sizeof(r2));
    memset((uint8_t *)r3, 0, sizeof(r3));
    memset((uint8_t *)r4, 0, sizeof(r4));
    memset((uint8_t *)acc, 0, sizeof(acc));
}

static void eval_chunk(int64_t *out, bn_t *x, uint32_t start, uint32_t chunk,
                       uint32_t digits, int point)
{
    for(uint32_t i=0; i<chunk; i++) {
        uint32_t i0 = start + i;
        uint32_t i1 = start + chunk + i;
        uint32_t i2 = start + 2*chunk + i;
        int64_t x0 = (i0 < digits) ? x[i0] : 0;
        int64_t x1 = (i1 < digits) ? x[i1] : 0;
        int64_t x2 = (i2 < digits) ? x[i2] : 0;

        if(point == 0) {
            out[i] = x0;
        } else if(point == 1) {
            out[i] = x0 + x1 + x2;
        } else if(point == -1) {
            out[i] = x0 - x1 + x2;
        } else {
            out[i] = x0 + 2*x1 + 4*x2;
        }
    }
}

static void mul_eval(i128_t *out, int64_t *x, int64_t *y, uint32_t len)
{
    for(uint32_t i=0; i<2*len; i++) {
        out[i] = 0;
    }

    for(uint32_t i=0; i<len; i++) {
        for(uint32_t j=0; j<len; j++) {
            out[i+j] += (i128_t)x[i] * y[j];
        }
    }
}

static void add_shifted(i128_t *out, uint32_t out_len, i128_t *x, uint32_t len,
                        uint32_t shift)
{
    for(uint32_t i=0; i<len && i+shift<out_len; i++) {
        out[i+shift] += x[i];
    }
}

static void normalize_product(bn_t *a, i128_t *x, uint32_t digits)
{
    i128_t carry = 0;
    i128_t base = ((i128_t)1) << BN_DIGIT_BITS;

    for(uint32_t i=0; i<2*digits; i++) {
        i128_t v = x[i] + carry;
        bn_t low = (bn_t)v;
        a[i] = low;
        carry = (v - (i128_t)low) / base;
    }
}

static void school_mul(bn_t *a, bn_t *b, bn_t *c, uint32_t digits)
{
    dbn_t uv;
    uint32_t i, j, k;

    bn_assign_zero(a, 2*digits);
    for(i=0; i<digits; i++) {
        bn_t carry = 0;
        for(j=0; j<digits; j++) {
            uv = (dbn_t)a[i+j] + (dbn_t)b[i] * c[j] + carry;
            a[i+j] = (bn_t)uv;
            carry = (bn_t)(uv >> BN_DIGIT_BITS);
        }
        k = i + digits;
        while(carry != 0) {
            uv = (dbn_t)a[k] + carry;
            a[k] = (bn_t)uv;
            carry = (bn_t)(uv >> BN_DIGIT_BITS);
            k++;
        }
    }
}
