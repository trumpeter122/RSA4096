/*****************************************************************************
Filename    : bignum.c
Author      : 
Date        : 
Description : 整理数据
*****************************************************************************/
#include <string.h>
#include "bignum.h"

#define KARATSUBA_THRESHOLD 16
#define KARATSUBA_MAX_DIGITS 256

static bn_t bn_sub_digit_mul(bn_t *a, bn_t *b, bn_t c, bn_t *d, uint32_t digits);
static uint32_t bn_digit_bits(bn_t a);
static void karatsuba_add_into(bn_t *a, uint32_t alen, bn_t *b, uint32_t blen, uint32_t offset);
static void karatsuba_sub_from(bn_t *a, uint32_t alen, bn_t *b, uint32_t blen);
static void school_mul(bn_t *a, bn_t *b, bn_t *c, uint32_t digits);
static void karatsuba_mul(bn_t *a, bn_t *b, bn_t *c, uint32_t digits);
static uint32_t next_power_of_two(uint32_t digits);
static bn_t montgomery_n0inv(bn_t n0);
static void montgomery_r2(bn_t *r2, bn_t *n, uint32_t digits);
static void montgomery_mul(bn_t *a, bn_t *b, bn_t *c, bn_t *n, uint32_t digits, bn_t n0inv);
static void reduce_operand(bn_t *a, bn_t *b, bn_t *n, uint32_t digits);

void bn_decode(bn_t *bn, uint32_t digits, uint8_t *hexarr, uint32_t size)
{
    bn_t t;
    int j;
    uint32_t i, u;
    for(i=0,j=size-1; i<digits && j>=0; i++) {
        t = 0;
        for(u=0; j>=0 && u<BN_DIGIT_BITS; j--, u+=8) {
            t |= ((bn_t)hexarr[j]) << u;
        }
        bn[i] = t;
    }

    for(; i<digits; i++) {
        bn[i] = 0;
    }
}

void bn_encode(uint8_t *hexarr, uint32_t size, bn_t *bn, uint32_t digits)
{
    bn_t t;
    int j;
    uint32_t i, u;

    for(i=0,j=size-1; i<digits && j>=0; i++) {
        t = bn[i];
        for(u=0; j>=0 && u<BN_DIGIT_BITS; j--, u+=8) {
            hexarr[j] = (uint8_t)(t >> u);
        }
    }

    for(; j>=0; j--) {
        hexarr[j] = 0;
    }
}

void bn_assign(bn_t *a, bn_t *b, uint32_t digits)
{
    uint32_t i;
    for(i=0; i<digits; i++) {
        a[i] = b[i];
    }
}

void bn_assign_zero(bn_t *a, uint32_t digits)
{
    uint32_t i;
    for(i=0; i<digits; i++) {
        a[i] = 0;
    }
}

bn_t bn_add(bn_t *a, bn_t *b, bn_t *c, uint32_t digits)
{
    bn_t ai, carry;
    uint32_t i;

    carry = 0;
    for(i=0; i<digits; i++) {
        if((ai = b[i] + carry) < carry) {
            ai = c[i];
        } else if((ai += c[i]) < c[i]) {
            carry = 1;
        } else {
            carry = 0;
        }
        a[i] = ai;
    }

    return carry;
}

bn_t bn_sub(bn_t *a, bn_t *b, bn_t *c, uint32_t digits)
{
    bn_t ai, borrow;
    uint32_t i;

    borrow = 0;
    for(i=0; i<digits; i++) {
        if((ai = b[i] - borrow) > (BN_MAX_DIGIT - borrow)) {
            ai = BN_MAX_DIGIT - c[i];
        } else if((ai -= c[i]) > (BN_MAX_DIGIT - c[i])) {
            borrow = 1;
        } else {
            borrow = 0;
        }
        a[i] = ai;
    }

    return borrow;
}

void bn_mul(bn_t *a, bn_t *b, bn_t *c, uint32_t digits)
{
    bn_t bb[KARATSUBA_MAX_DIGITS], cc[KARATSUBA_MAX_DIGITS];
    bn_t tt[2*KARATSUBA_MAX_DIGITS];
    uint32_t padded;

    padded = next_power_of_two(digits);
    if(padded > KARATSUBA_MAX_DIGITS) {
        padded = KARATSUBA_MAX_DIGITS;
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

void bn_div(bn_t *a, bn_t *b, bn_t *c, uint32_t cdigits, bn_t *d, uint32_t ddigits)
{
    dbn_t tmp;
    bn_t ai, t, cc[2*BN_MAX_DIGITS+1], dd[BN_MAX_DIGITS];
    int i;
    uint32_t dddigits, shift;

    dddigits = bn_digits(d, ddigits);
    if(dddigits == 0)
        return;

    shift = BN_DIGIT_BITS - bn_digit_bits(d[dddigits-1]);
    bn_assign_zero(cc, dddigits);
    cc[cdigits] = bn_shift_l(cc, c, shift, cdigits);
    bn_shift_l(dd, d, shift, dddigits);
    t = dd[dddigits-1];

    bn_assign_zero(a, cdigits);
    i = cdigits - dddigits;
    for(; i>=0; i--) {
        if(t == BN_MAX_DIGIT) {
            ai = cc[i+dddigits];
        } else {
            tmp = cc[i+dddigits-1];
            tmp += (dbn_t)cc[i+dddigits] << BN_DIGIT_BITS;
            ai = tmp / (t + 1);
        }

        cc[i+dddigits] -= bn_sub_digit_mul(&cc[i], &cc[i], ai, dd, dddigits);
        // printf("cc[%d]: %08X\n", i, cc[i+dddigits]);
        while(cc[i+dddigits] || (bn_cmp(&cc[i], dd, dddigits) >= 0)) {
            ai++;
            cc[i+dddigits] -= bn_sub(&cc[i], &cc[i], dd, dddigits);
        }
        a[i] = ai;
        // printf("ai[%d]: %08X\n", i, ai);
    }

    bn_assign_zero(b, ddigits);
    bn_shift_r(b, cc, shift, dddigits);

    // Clear potentially sensitive information
    memset((uint8_t *)cc, 0, sizeof(cc));
    memset((uint8_t *)dd, 0, sizeof(dd));
}

bn_t bn_shift_l(bn_t *a, bn_t *b, uint32_t c, uint32_t digits)
{
    bn_t bi, carry;
    uint32_t i, t;

    if(c >= BN_DIGIT_BITS)
        return 0;

    t = BN_DIGIT_BITS - c;
    carry = 0;
    for(i=0; i<digits; i++) {
        bi = b[i];
        a[i] = (bi << c) | carry;
        carry = c ? (bi >> t) : 0;
    }

    return carry;
}

bn_t bn_shift_r(bn_t *a, bn_t *b, uint32_t c, uint32_t digits)
{
    bn_t bi, carry;
    int i;
    uint32_t t;

    if(c >= BN_DIGIT_BITS)
        return 0;

    t = BN_DIGIT_BITS - c;
    carry = 0;
    i = digits - 1;
    for(; i>=0; i--) {
        bi = b[i];
        a[i] = (bi >> c) | carry;
        carry = c ? (bi << t) : 0;
    }

    return carry;
}

void bn_mod(bn_t *a, bn_t *b, uint32_t bdigits, bn_t *c, uint32_t cdigits)
{
    bn_t t[2*BN_MAX_DIGITS] = {0};

    bn_div(t, a, b, bdigits, c, cdigits);

    // Clear potentially sensitive information
    memset((uint8_t *)t, 0, sizeof(t));
}

void bn_mod_mul(bn_t *a, bn_t *b, bn_t *c, bn_t *d, uint32_t digits)
{
    bn_t bb[BN_MAX_DIGITS], cc[BN_MAX_DIGITS], one[BN_MAX_DIGITS];
    bn_t r2[BN_MAX_DIGITS], bm[BN_MAX_DIGITS], cm[BN_MAX_DIGITS], tm[BN_MAX_DIGITS];
    bn_t n0inv;

    if(digits == 0 || d[0] == 0 || ((d[0] & 1) == 0)) {
        bn_t t[2*BN_MAX_DIGITS];

        bn_mul(t, b, c, digits);
        bn_mod(a, t, 2*digits, d, digits);

        // Clear potentially sensitive information
        memset((uint8_t *)t, 0, sizeof(t));
        return;
    }

    reduce_operand(bb, b, d, digits);
    reduce_operand(cc, c, d, digits);

    montgomery_r2(r2, d, digits);
    n0inv = montgomery_n0inv(d[0]);
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



void bn_mod_exp(bn_t *a, bn_t *b, bn_t *c, uint32_t cdigits, bn_t *d, uint32_t ddigits)
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

int bn_cmp(bn_t *a, bn_t *b, uint32_t digits)
{
    int i;
    for(i=digits-1; i>=0; i--) {
        if(a[i] > b[i])     return 1;
        if(a[i] < b[i])     return -1;
    }

    return 0;
}

uint32_t bn_digits(bn_t *a, uint32_t digits)
{
    int i;
    for(i=digits-1; i>=0; i--) {
        if(a[i])    break;
    }

    return (i + 1);
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

    karatsuba_sub_from(z1, 2*(half+1), z0, 2*half);
    karatsuba_sub_from(z1, 2*(half+1), z2, 2*half);

    karatsuba_add_into(a, 2*digits, z0, 2*half, 0);
    karatsuba_add_into(a, 2*digits, z1, 2*(half+1), half);
    karatsuba_add_into(a, 2*digits, z2, 2*half, 2*half);

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

static void karatsuba_add_into(bn_t *a, uint32_t alen, bn_t *b, uint32_t blen, uint32_t offset)
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

static void karatsuba_sub_from(bn_t *a, uint32_t alen, bn_t *b, uint32_t blen)
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

static bn_t bn_sub_digit_mul(bn_t *a, bn_t *b, bn_t c, bn_t *d, uint32_t digits)
{
    dbn_t result;
    bn_t borrow, rh, rl;
    uint32_t i;

    if(c == 0)
        return 0;

    borrow = 0;
    for(i=0; i<digits; i++) {
        result = (dbn_t)c * d[i];
        rl = result & BN_MAX_DIGIT;
        rh = (result >> BN_DIGIT_BITS) & BN_MAX_DIGIT;
        if((a[i] = b[i] - borrow) > (BN_MAX_DIGIT - borrow)) {
            borrow = 1;
        } else {
            borrow = 0;
        }
        if((a[i] -= rl) > (BN_MAX_DIGIT - rl)) {
            borrow++;
        }
        borrow += rh;
    }

    return borrow;
}

static uint32_t bn_digit_bits(bn_t a)
{
    uint32_t i;
    for(i=0; i<BN_DIGIT_BITS; i++) {
        if(a == 0)  break;
        a >>= 1;
    }

    return i;
}
