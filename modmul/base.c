#include <string.h>
#include "modmul.h"

void modmul_bn_mod_mul(bn_t *a, bn_t *b, bn_t *c, bn_t *d, uint32_t digits)
{
    bn_t t[2*BN_MAX_DIGITS];

    bn_mul(t, b, c, digits);
    bn_mod(a, t, 2*digits, d, digits);

    // Clear potentially sensitive information
    memset((uint8_t *)t, 0, sizeof(t));
}
