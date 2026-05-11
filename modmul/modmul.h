#ifndef __MODMUL_H__
#define __MODMUL_H__

#include "../bignum.h"

void modmul_bn_mod_mul(bn_t *a, bn_t *b, bn_t *c, bn_t *d, uint32_t digits);

#endif  // __MODMUL_H__
