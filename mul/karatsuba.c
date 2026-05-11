#include <string.h>
#include "mul.h"

#define KARATSUBA_THRESHOLD 16
#define KARATSUBA_MAX_DIGITS 256

static void add_into(bn_t *a, uint32_t alen, bn_t *b, uint32_t blen, uint32_t offset);
static void sub_from(bn_t *a, uint32_t alen, bn_t *b, uint32_t blen);
static void school_mul(bn_t *a, bn_t *b, bn_t *c, uint32_t digits);
static void karatsuba_mul(bn_t *a, bn_t *b, bn_t *c, uint32_t digits);
static uint32_t next_power_of_two(uint32_t digits);

void mul_bn_mul(bn_t *a, bn_t *b, bn_t *c, uint32_t digits)
{
    bn_t bb[KARATSUBA_MAX_DIGITS], cc[KARATSUBA_MAX_DIGITS];
    bn_t tt[2*KARATSUBA_MAX_DIGITS];
    uint32_t padded;

    padded = next_power_of_two(digits);
    if(padded > KARATSUBA_MAX_DIGITS) {
        school_mul(a, b, c, digits);
        return;
    }

    bn_assign_zero(bb, padded);
    bn_assign_zero(cc, padded);
    bn_assign_zero(tt, 2*padded);
    bn_assign(bb, b, digits);
    bn_assign(cc, c, digits);

    karatsuba_mul(tt, bb, cc, padded);
    bn_assign(a, tt, 2*digits);

    // Clear potentially sensitive information
    memset((uint8_t *)bb, 0, sizeof(bb));
    memset((uint8_t *)cc, 0, sizeof(cc));
    memset((uint8_t *)tt, 0, sizeof(tt));
}

static void karatsuba_mul(bn_t *a, bn_t *b, bn_t *c, uint32_t digits)
{
    uint32_t half;

    bn_assign_zero(a, 2*digits);
    if(digits <= KARATSUBA_THRESHOLD) {
        school_mul(a, b, c, digits);
        return;
    }

    half = digits / 2;

    bn_t z0[2*half], z2[2*half];
    bn_t bsum[half+1], csum[half+1], z1[2*(half+1)];

    bn_assign_zero(z0, 2*half);
    bn_assign_zero(z2, 2*half);
    bn_assign_zero(bsum, half+1);
    bn_assign_zero(csum, half+1);
    bn_assign_zero(z1, 2*(half+1));

    karatsuba_mul(z0, b, c, half);
    karatsuba_mul(z2, b+half, c+half, half);

    bsum[half] = bn_add(bsum, b, b+half, half);
    csum[half] = bn_add(csum, c, c+half, half);
    school_mul(z1, bsum, csum, half+1);

    sub_from(z1, 2*(half+1), z0, 2*half);
    sub_from(z1, 2*(half+1), z2, 2*half);

    add_into(a, 2*digits, z0, 2*half, 0);
    add_into(a, 2*digits, z1, 2*(half+1), half);
    add_into(a, 2*digits, z2, 2*half, 2*half);

    // Clear potentially sensitive information
    memset((uint8_t *)z0, 0, sizeof(z0));
    memset((uint8_t *)z2, 0, sizeof(z2));
    memset((uint8_t *)bsum, 0, sizeof(bsum));
    memset((uint8_t *)csum, 0, sizeof(csum));
    memset((uint8_t *)z1, 0, sizeof(z1));
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

static void add_into(bn_t *a, uint32_t alen, bn_t *b, uint32_t blen, uint32_t offset)
{
    dbn_t sum;
    bn_t carry;
    uint32_t i;

    carry = 0;
    for(i=0; i<blen && i+offset<alen; i++) {
        sum = (dbn_t)a[i+offset] + b[i] + carry;
        a[i+offset] = (bn_t)sum;
        carry = (bn_t)(sum >> BN_DIGIT_BITS);
    }
    while(carry != 0 && i+offset<alen) {
        sum = (dbn_t)a[i+offset] + carry;
        a[i+offset] = (bn_t)sum;
        carry = (bn_t)(sum >> BN_DIGIT_BITS);
        i++;
    }
}

static void sub_from(bn_t *a, uint32_t alen, bn_t *b, uint32_t blen)
{
    dbn_t diff;
    bn_t borrow;
    uint32_t i;

    borrow = 0;
    for(i=0; i<alen; i++) {
        bn_t bi = (i < blen) ? b[i] : 0;
        diff = (dbn_t)a[i] - bi - borrow;
        a[i] = (bn_t)diff;
        borrow = (diff >> 63) & 1;
    }
}

static uint32_t next_power_of_two(uint32_t digits)
{
    uint32_t n;

    n = 1;
    while(n < digits) {
        n <<= 1;
    }

    return n;
}
