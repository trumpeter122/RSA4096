#ifndef __MUL_H__
#define __MUL_H__

#include "../bignum.h"

void mul_bn_mul(bn_t *a, bn_t *b, bn_t *c, uint32_t digits);

#endif  // __MUL_H__
