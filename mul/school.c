#include "mul.h"

void mul_bn_mul(bn_t *a, bn_t *b, bn_t *c, uint32_t digits)
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
