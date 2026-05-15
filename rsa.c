/*****************************************************************************
Filename    : rsa.c
Author      : 
Date        : 
Description :
*****************************************************************************/
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#include "rsa.h"
#include "bignum.h"

#define MONT_WINDOW_BITS 6
#define MONT_WINDOW_TABLE (1U << (MONT_WINDOW_BITS - 1))
#define MONT_WINDOW_MIN_BITS 256

typedef struct {
    rsa_sk_t *sk;
    bn_t p[BN_MAX_DIGITS], q[BN_MAX_DIGITS], dp[BN_MAX_DIGITS], dq[BN_MAX_DIGITS];
    bn_t qinv[BN_MAX_DIGITS], rp[BN_MAX_DIGITS], rrp[BN_MAX_DIGITS], rq[BN_MAX_DIGITS], rrq[BN_MAX_DIGITS];
    bn_t p_n0inv, q_n0inv;
    uint32_t pdigits, qdigits, dpdigits, dqdigits;
} mont_private_cache_t;

typedef struct {
    rsa_pk_t *pk;
    bn_t n[BN_MAX_DIGITS], e[BN_MAX_DIGITS], r[BN_MAX_DIGITS], rr[BN_MAX_DIGITS], one[BN_MAX_DIGITS];
    bn_t n0inv;
    uint32_t ndigits, edigits;
} mont_public_cache_t;

static int private_block_operation(uint8_t *out, uint32_t *out_len, uint8_t *in, uint32_t in_len, rsa_sk_t *sk);
static int public_block_operation(uint8_t *out, uint32_t *out_len, uint8_t *in, uint32_t in_len, rsa_pk_t *pk);
static int rsa_private_compute_mont_crt(bn_t *out, bn_t *in, rsa_sk_t *sk, bn_t *n, uint32_t ndigits);
static int rsa_public_compute_mont(bn_t *out, bn_t *in, rsa_pk_t *pk);
static mont_private_cache_t *get_private_cache(rsa_sk_t *sk);
static mont_public_cache_t *get_public_cache(rsa_pk_t *pk);
static bn_t montgomery_n0inv(bn_t n0);
static uint32_t exponent_bit(bn_t *a, uint32_t bit);
static uint32_t highest_bit_index(bn_t *a, uint32_t digits);
static void montgomery_mul(bn_t *a, bn_t *b, bn_t *c, bn_t *n, uint32_t digits, bn_t n0inv);
static void montgomery_mul_64(bn_t *a, bn_t *b, bn_t *c, bn_t *n, bn_t n0inv);
static void montgomery_mul_128(bn_t *a, bn_t *b, bn_t *c, bn_t *n, bn_t n0inv);
static void montgomery_mul_var(bn_t *a, bn_t *b, bn_t *c, bn_t *n, uint32_t digits, bn_t n0inv);
static void montgomery_precompute_constants(bn_t *r, bn_t *rr, bn_t *n, uint32_t digits);
static void montgomery_mod_mul_precomp(bn_t *a, bn_t *b, bn_t *c, bn_t *n,
                                       uint32_t digits, bn_t n0inv, bn_t *rr);
static void montgomery_mod_exp_precomp(bn_t *a, bn_t *b, bn_t *c, uint32_t cdigits,
                                       bn_t *n, uint32_t ndigits, bn_t n0inv, bn_t *r, bn_t *rr);
static void montgomery_mod_exp_binary(bn_t *a, bn_t *base_m, bn_t *c, uint32_t cdigits,
                                      bn_t *n, uint32_t ndigits, bn_t n0inv);
static void montgomery_mod_exp_window(bn_t *a, bn_t *base_m, bn_t *c, uint32_t cdigits,
                                      bn_t *n, uint32_t ndigits, bn_t n0inv, bn_t *r);
static void montgomery_mod_exp_65537(bn_t *a, bn_t *b, bn_t *n, uint32_t ndigits,
                                     bn_t n0inv, bn_t *rr, bn_t *one);
static void classic_mod_exp(bn_t *a, bn_t *b, bn_t *c, uint32_t cdigits, bn_t *n, uint32_t ndigits);

int rsa_private_encrypt_any_len(uint8_t *out, uint32_t *out_len, uint8_t *in, uint32_t in_len, rsa_sk_t *sk){
	int status=0;
	uint32_t len=0;
	uint8_t *tmp_o=out;
	*out_len=0;
	for(uint32_t i=0;i<in_len && status==0;i+=(RSA_MAX_MODULUS_LEN-11)){
		if((in_len-i)>(RSA_MAX_MODULUS_LEN-11)){
			status=rsa_private_encrypt(tmp_o,&len,in+i,RSA_MAX_MODULUS_LEN-11,sk);
		}
		else{
			status=rsa_private_encrypt(tmp_o,&len,in+i,in_len-i,sk);
			*out_len+=len;
			break;
		}
		tmp_o=tmp_o+len;
		*out_len+=len;
	}
	tmp_o=NULL;
	return status;
}

int rsa_public_encrypt_any_len(uint8_t *out, uint32_t *out_len, uint8_t *in, uint32_t in_len, rsa_pk_t *pk){
	int status=0;
	uint32_t len=0;
	uint8_t *tmp_o=out;
	*out_len=0;
	for(uint32_t i=0;i<in_len && status==0;i+=(RSA_MAX_MODULUS_LEN-11)){
		if((in_len-i)>(RSA_MAX_MODULUS_LEN-11)){
			status=rsa_public_encrypt(tmp_o,&len,in+i,RSA_MAX_MODULUS_LEN-11,pk);
			tmp_o=tmp_o+len;
			*out_len+=len;
		}
		else{
			status=rsa_public_encrypt(tmp_o,&len,in+i,in_len-i,pk);
			*out_len+=len;
			break;
		}
	}
	tmp_o=NULL;
	return status;
}


int rsa_private_decrypt_any_len(uint8_t *out, uint32_t *out_len, uint8_t *in, uint32_t in_len, rsa_sk_t *sk){
	int status=0;
	uint32_t len=0;
	uint8_t *tmp_o=out;
	*out_len=0;
	for(uint32_t i=0;i<in_len && status==0;i+=RSA_MAX_MODULUS_LEN){
		if((in_len-i)>RSA_MAX_MODULUS_LEN){
			status=rsa_private_decrypt(tmp_o,&len,in+i,RSA_MAX_MODULUS_LEN,sk);
			tmp_o=tmp_o+len;
			*out_len+=len;
		}
		else{
			status=rsa_private_decrypt(tmp_o,&len,in+i,in_len-i,sk);
			*out_len+=len;
			break;
		}
	}
	tmp_o=NULL;
	return status;
}

int rsa_private_encrypt(uint8_t *out, uint32_t *out_len, uint8_t *in, uint32_t in_len, rsa_sk_t *sk)
{
    int status;
    uint8_t pkcs_block[RSA_MAX_MODULUS_LEN];
    uint32_t i, modulus_len;

    modulus_len = (sk->bits + 7) / 8;
    if(in_len + 11 > modulus_len)
        return ERR_WRONG_LEN;

    pkcs_block[0] = 0;
    pkcs_block[1] = 1;
    for(i=2; i<modulus_len-in_len-1; i++) {
        pkcs_block[i] = 0xFF;
    }

    pkcs_block[i++] = 0;

    memcpy((uint8_t *)&pkcs_block[i], (uint8_t *)in, in_len);

    status = private_block_operation(out, out_len, pkcs_block, modulus_len, sk);

    // Clear potentially sensitive information
    memset((uint8_t *)pkcs_block, 0, sizeof(pkcs_block));

    return status;
}

int rsa_private_decrypt(uint8_t *out, uint32_t *out_len, uint8_t *in, uint32_t in_len, rsa_sk_t *sk)
{
    int status;
    uint8_t pkcs_block[RSA_MAX_MODULUS_LEN];
    uint32_t i, modulus_len, pkcs_block_len;

    modulus_len = (sk->bits + 7) / 8;
    if(in_len > modulus_len)
        return ERR_WRONG_LEN;

    status = private_block_operation(pkcs_block, &pkcs_block_len, in, in_len, sk);
    if(status != 0)
        return status;

    if(pkcs_block_len != modulus_len)
        return ERR_WRONG_LEN;

    if((pkcs_block[0] != 0) || (pkcs_block[1] != 2))
        return ERR_WRONG_DATA;

    for(i=2; i<modulus_len-1; i++) {
        if(pkcs_block[i] == 0)  break;
    }

    i++;
    if(i >= modulus_len)
        return ERR_WRONG_DATA;
    *out_len = modulus_len - i;
    if(*out_len + 11 > modulus_len)
        return ERR_WRONG_DATA;
    memcpy((uint8_t *)out, (uint8_t *)&pkcs_block[i], *out_len);
    // Clear potentially sensitive information
    memset((uint8_t *)pkcs_block, 0, sizeof(pkcs_block));

    return status;
}

static int private_block_operation(uint8_t *out, uint32_t *out_len, uint8_t *in, uint32_t in_len, rsa_sk_t *sk)
{
    uint32_t ndigits;
    bn_t c[BN_MAX_DIGITS], m[BN_MAX_DIGITS], n[BN_MAX_DIGITS];

    bn_decode(c, BN_MAX_DIGITS, in, in_len);
    bn_decode(n, BN_MAX_DIGITS, sk->modulus, RSA_MAX_MODULUS_LEN);

    ndigits = bn_digits(n, BN_MAX_DIGITS);

    if(bn_cmp(c, n, ndigits) >= 0)
        return ERR_WRONG_DATA;

    rsa_private_compute_mont_crt(m, c, sk, n, ndigits);

    *out_len = (sk->bits + 7) / 8;
    bn_encode(out, *out_len, m, ndigits);

    // Clear potentially sensitive information
    memset((uint8_t *)c, 0, sizeof(c));
    memset((uint8_t *)m, 0, sizeof(m));

    return 0;
}

// Public encryption
void generate_rand(uint8_t *block, uint32_t block_len)
{
    uint32_t i;
	srand ((unsigned)time(NULL));   // real rand message
    for(i=0; i<block_len; i++) {
        block[i] = rand();
		while(block[i]==0)
			block[i]=rand();
    }
}


// int rsa_public_decrypt_any_len(uint8_t *out, uint32_t *out_len, uint8_t *in, uint32_t in_len, rsa_pk_t *pk){
// 	int status=0;
// 	int len=0;
// 	uint8_t *tmp_o=out;
// 	*out_len=0;
// 	for(int i=0;i<in_len && status==0;i+=RSA_MAX_MODULUS_LEN){
// 		if((in_len-i)>RSA_MAX_MODULUS_LEN){
// 			status=rsa_public_decrypt(tmp_o,&len,in,RSA_MAX_MODULUS_LEN,pk);
// 		}
// 		else{
// 			status=rsa_public_decrypt(tmp_o,&len,in,in_len-i,pk);
// 			break;
// 		}
// 		tmp_o=tmp_o+len;
// 		*out_len+=len;
// 	}
// 	tmp_o=NULL;
// 	free(tmp_o);
// 	// *out_len=len;
// 	return status;
// }

int rsa_public_encrypt(uint8_t *out, uint32_t *out_len, uint8_t *in, uint32_t in_len, rsa_pk_t *pk)
{
    int status;
    uint8_t byte, pkcs_block[RSA_MAX_MODULUS_LEN];
    uint32_t i, modulus_len;

    modulus_len = (pk->bits + 7) / 8;
    if(in_len + 11 > modulus_len) {//padding len
        return ERR_WRONG_LEN;
    }

    pkcs_block[0] = 0;
    pkcs_block[1] = 2;
    for(i=2; i<modulus_len-in_len-1; i++) {
        do {
            generate_rand(&byte, 1);
        } while(byte == 0);
        pkcs_block[i] = byte;
    }
    pkcs_block[i++] = 0;

    memcpy((uint8_t *)&pkcs_block[i], (uint8_t *)in, in_len);
    status = public_block_operation(out, out_len, pkcs_block, modulus_len, pk);
    // Clear potentially sensitive information
    byte = 0;
    memset((uint8_t *)pkcs_block, 0, sizeof(pkcs_block));

    return status;
}

// int rsa_public_decrypt(uint8_t *out, uint32_t *out_len, uint8_t *in, uint32_t in_len, rsa_pk_t *pk) {
//     int status;
//     uint8_t pkcs_block[RSA_MAX_MODULUS_LEN];
//     uint32_t i, modulus_len, pkcs_block_len;

//     modulus_len = (pk->bits + 7) / 8;
//     if (in_len > modulus_len)
//         return ERR_WRONG_LEN;

//     status = public_block_operation(pkcs_block, &pkcs_block_len, in, in_len, pk);
//     if (status != 0)
//         return status;

//     if (pkcs_block_len != modulus_len)
//         return ERR_WRONG_LEN;

//     if ((pkcs_block[0] != 0) || (pkcs_block[1] != 1))
//         return ERR_WRONG_DATA;

//     for (i = 2; i < modulus_len - 1; i++) {
//         if (pkcs_block[i] != 0xFF) break;
//     }

//     if (pkcs_block[i++] != 0)
//         return ERR_WRONG_DATA;

//     *out_len = modulus_len - i;
//     if (*out_len + 11 > modulus_len)
//         return ERR_WRONG_DATA;

//     memcpy((uint8_t *) out, (uint8_t *) &pkcs_block[i], *out_len);

//     // Clear potentially sensitive information
//     memset((uint8_t *) pkcs_block, 0, sizeof(pkcs_block));

//     return status;
// }

static int public_block_operation(uint8_t *out, uint32_t *out_len, uint8_t *in, uint32_t in_len, rsa_pk_t *pk)
{
    uint32_t ndigits;
    bn_t c[BN_MAX_DIGITS], m[BN_MAX_DIGITS], n[BN_MAX_DIGITS];

    bn_decode(m, BN_MAX_DIGITS, in, in_len);
    bn_decode(n, BN_MAX_DIGITS, pk->modulus, RSA_MAX_MODULUS_LEN);

    ndigits = bn_digits(n, BN_MAX_DIGITS);

    if(bn_cmp(m, n, ndigits) >= 0) {
        return ERR_WRONG_DATA;
    }

    rsa_public_compute_mont(c, m, pk);

    *out_len = (pk->bits + 7) / 8;
    bn_encode(out, *out_len, c, ndigits);

    // Clear potentially sensitive information
    memset((uint8_t *)c, 0, sizeof(c));
    memset((uint8_t *)m, 0, sizeof(m));

    return 0;
}

static int rsa_private_compute_mont_crt(bn_t *out, bn_t *in, rsa_sk_t *sk, bn_t *n, uint32_t ndigits)
{
    mont_private_cache_t *cache;
    bn_t cp[BN_MAX_DIGITS], cq[BN_MAX_DIGITS];
    bn_t m1[BN_MAX_DIGITS], m2[BN_MAX_DIGITS], diff[BN_MAX_DIGITS], h[BN_MAX_DIGITS];
    bn_t qh[2*BN_MAX_DIGITS], m2n[BN_MAX_DIGITS];

    cache = get_private_cache(sk);

    bn_mod(cp, in, ndigits, cache->p, cache->pdigits);
    bn_mod(cq, in, ndigits, cache->q, cache->qdigits);
    montgomery_mod_exp_precomp(m1, cp, cache->dp, cache->dpdigits, cache->p, cache->pdigits,
                               cache->p_n0inv, cache->rp, cache->rrp);
    montgomery_mod_exp_precomp(m2, cq, cache->dq, cache->dqdigits, cache->q, cache->qdigits,
                               cache->q_n0inv, cache->rq, cache->rrq);

    if(bn_cmp(m1, m2, cache->pdigits) >= 0) {
        bn_sub(diff, m1, m2, cache->pdigits);
    } else {
        bn_sub(diff, m2, m1, cache->pdigits);
        bn_sub(diff, cache->p, diff, cache->pdigits);
    }

    montgomery_mod_mul_precomp(h, diff, cache->qinv, cache->p, cache->pdigits,
                               cache->p_n0inv, cache->rrp);
    bn_assign_zero(qh, 2*BN_MAX_DIGITS);
    bn_mul(qh, cache->q, h, cache->pdigits);

    bn_assign_zero(m2n, ndigits);
    bn_assign(m2n, m2, cache->qdigits);
    bn_add(out, qh, m2n, ndigits);
    if(bn_cmp(out, n, ndigits) >= 0) {
        bn_sub(out, out, n, ndigits);
    }

    // Clear potentially sensitive information
    memset((uint8_t *)cp, 0, sizeof(cp));
    memset((uint8_t *)cq, 0, sizeof(cq));
    memset((uint8_t *)m1, 0, sizeof(m1));
    memset((uint8_t *)m2, 0, sizeof(m2));
    memset((uint8_t *)diff, 0, sizeof(diff));
    memset((uint8_t *)h, 0, sizeof(h));
    memset((uint8_t *)qh, 0, sizeof(qh));
    memset((uint8_t *)m2n, 0, sizeof(m2n));

    return 0;
}

static int rsa_public_compute_mont(bn_t *out, bn_t *in, rsa_pk_t *pk)
{
    mont_public_cache_t *cache;

    cache = get_public_cache(pk);
    if(cache->edigits == 1 && cache->e[0] == 65537) {
        montgomery_mod_exp_65537(out, in, cache->n, cache->ndigits,
                                 cache->n0inv, cache->rr, cache->one);
    } else {
        montgomery_mod_exp_precomp(out, in, cache->e, cache->edigits, cache->n, cache->ndigits,
                                   cache->n0inv, cache->r, cache->rr);
    }

    return 0;
}

static mont_private_cache_t *get_private_cache(rsa_sk_t *sk)
{
    static mont_private_cache_t cache;

    if(cache.sk == sk) {
        return &cache;
    }

    memset((uint8_t *)&cache, 0, sizeof(cache));
    cache.sk = sk;

    bn_decode(cache.p, BN_MAX_DIGITS, sk->prime1, RSA_MAX_PRIME_LEN);
    bn_decode(cache.q, BN_MAX_DIGITS, sk->prime2, RSA_MAX_PRIME_LEN);
    bn_decode(cache.dp, BN_MAX_DIGITS, sk->prime_exponent1, RSA_MAX_PRIME_LEN);
    bn_decode(cache.dq, BN_MAX_DIGITS, sk->prime_exponent2, RSA_MAX_PRIME_LEN);
    bn_decode(cache.qinv, BN_MAX_DIGITS, sk->coefficient, RSA_MAX_PRIME_LEN);

    cache.pdigits = bn_digits(cache.p, BN_MAX_DIGITS);
    cache.qdigits = bn_digits(cache.q, BN_MAX_DIGITS);
    cache.dpdigits = bn_digits(cache.dp, cache.pdigits);
    cache.dqdigits = bn_digits(cache.dq, cache.qdigits);
    cache.p_n0inv = montgomery_n0inv(cache.p[0]);
    cache.q_n0inv = montgomery_n0inv(cache.q[0]);
    montgomery_precompute_constants(cache.rp, cache.rrp, cache.p, cache.pdigits);
    montgomery_precompute_constants(cache.rq, cache.rrq, cache.q, cache.qdigits);

    return &cache;
}

static mont_public_cache_t *get_public_cache(rsa_pk_t *pk)
{
    static mont_public_cache_t cache;

    if(cache.pk == pk) {
        return &cache;
    }

    memset((uint8_t *)&cache, 0, sizeof(cache));
    cache.pk = pk;

    bn_decode(cache.n, BN_MAX_DIGITS, pk->modulus, RSA_MAX_MODULUS_LEN);
    bn_decode(cache.e, BN_MAX_DIGITS, pk->exponent, RSA_MAX_MODULUS_LEN);
    cache.ndigits = bn_digits(cache.n, BN_MAX_DIGITS);
    cache.edigits = bn_digits(cache.e, BN_MAX_DIGITS);
    cache.n0inv = montgomery_n0inv(cache.n[0]);
    montgomery_precompute_constants(cache.r, cache.rr, cache.n, cache.ndigits);
    BN_ASSIGN_DIGIT(cache.one, 1, cache.ndigits);

    return &cache;
}

static void montgomery_mod_exp_precomp(bn_t *a, bn_t *b, bn_t *c, uint32_t cdigits,
                                       bn_t *n, uint32_t ndigits, bn_t n0inv, bn_t *r, bn_t *rr)
{
    bn_t base[BN_MAX_DIGITS], base_m[BN_MAX_DIGITS], one[BN_MAX_DIGITS];
    uint32_t bits;

    if(ndigits == 0 || n[0] == 0 || ((n[0] & 1) == 0)) {
        classic_mod_exp(a, b, c, cdigits, n, ndigits);
        return;
    }

    cdigits = bn_digits(c, cdigits);
    if(cdigits == 0) {
        BN_ASSIGN_DIGIT(a, 1, ndigits);
        return;
    }

    if(bn_cmp(b, n, ndigits) >= 0) {
        bn_mod(base, b, ndigits, n, ndigits);
    } else {
        bn_assign(base, b, ndigits);
    }

    montgomery_mul(base_m, base, rr, n, ndigits, n0inv);
    bits = highest_bit_index(c, cdigits) + 1;
    if(bits >= MONT_WINDOW_MIN_BITS) {
        montgomery_mod_exp_window(a, base_m, c, cdigits, n, ndigits, n0inv, r);
    } else {
        montgomery_mod_exp_binary(a, base_m, c, cdigits, n, ndigits, n0inv);
    }

    BN_ASSIGN_DIGIT(one, 1, ndigits);
    montgomery_mul(a, a, one, n, ndigits, n0inv);

    // Clear potentially sensitive information
    memset((uint8_t *)base, 0, sizeof(base));
    memset((uint8_t *)base_m, 0, sizeof(base_m));
    memset((uint8_t *)one, 0, sizeof(one));
}

static void montgomery_mod_exp_binary(bn_t *a, bn_t *base_m, bn_t *c, uint32_t cdigits,
                                      bn_t *n, uint32_t ndigits, bn_t n0inv)
{
    bn_t result[BN_MAX_DIGITS];
    int bit;

    bit = (int)highest_bit_index(c, cdigits);
    bn_assign(result, base_m, ndigits);
    for(bit--; bit>=0; bit--) {
        montgomery_mul(result, result, result, n, ndigits, n0inv);
        if(exponent_bit(c, (uint32_t)bit)) {
            montgomery_mul(result, result, base_m, n, ndigits, n0inv);
        }
    }

    bn_assign(a, result, ndigits);

    // Clear potentially sensitive information
    memset((uint8_t *)result, 0, sizeof(result));
}

static void montgomery_mod_exp_window(bn_t *a, bn_t *base_m, bn_t *c, uint32_t cdigits,
                                      bn_t *n, uint32_t ndigits, bn_t n0inv, bn_t *r)
{
    bn_t table[MONT_WINDOW_TABLE][BN_MAX_DIGITS], base2[BN_MAX_DIGITS], result[BN_MAX_DIGITS];
    int bit;

    bn_assign(table[0], base_m, ndigits);
    montgomery_mul(base2, base_m, base_m, n, ndigits, n0inv);
    for(uint32_t i=1; i<MONT_WINDOW_TABLE; i++) {
        montgomery_mul(table[i], table[i-1], base2, n, ndigits, n0inv);
    }

    bn_assign(result, r, ndigits);
    bit = (int)highest_bit_index(c, cdigits);
    while(bit >= 0) {
        if(!exponent_bit(c, (uint32_t)bit)) {
            montgomery_mul(result, result, result, n, ndigits, n0inv);
            bit--;
        } else {
            int low = bit - MONT_WINDOW_BITS + 1;
            uint32_t value = 0;

            if(low < 0) {
                low = 0;
            }
            while(!exponent_bit(c, (uint32_t)low)) {
                low++;
            }

            for(int j=bit; j>=low; j--) {
                value = (value << 1) | exponent_bit(c, (uint32_t)j);
                montgomery_mul(result, result, result, n, ndigits, n0inv);
            }
            montgomery_mul(result, result, table[value >> 1], n, ndigits, n0inv);
            bit = low - 1;
        }
    }

    bn_assign(a, result, ndigits);

    // Clear potentially sensitive information
    memset((uint8_t *)table, 0, sizeof(table));
    memset((uint8_t *)base2, 0, sizeof(base2));
    memset((uint8_t *)result, 0, sizeof(result));
}

static void montgomery_mod_exp_65537(bn_t *a, bn_t *b, bn_t *n, uint32_t ndigits,
                                     bn_t n0inv, bn_t *rr, bn_t *one)
{
    bn_t base[BN_MAX_DIGITS], base_m[BN_MAX_DIGITS], result[BN_MAX_DIGITS];

    if(bn_cmp(b, n, ndigits) >= 0) {
        bn_mod(base, b, ndigits, n, ndigits);
    } else {
        bn_assign(base, b, ndigits);
    }

    montgomery_mul(base_m, base, rr, n, ndigits, n0inv);
    bn_assign(result, base_m, ndigits);
    for(uint32_t i=0; i<16; i++) {
        montgomery_mul(result, result, result, n, ndigits, n0inv);
    }
    montgomery_mul(result, result, base_m, n, ndigits, n0inv);
    montgomery_mul(a, result, one, n, ndigits, n0inv);

    // Clear potentially sensitive information
    memset((uint8_t *)base, 0, sizeof(base));
    memset((uint8_t *)base_m, 0, sizeof(base_m));
    memset((uint8_t *)result, 0, sizeof(result));
}

static void montgomery_mod_mul_precomp(bn_t *a, bn_t *b, bn_t *c, bn_t *n,
                                       uint32_t digits, bn_t n0inv, bn_t *rr)
{
    bn_t bm[BN_MAX_DIGITS], cm[BN_MAX_DIGITS], tm[BN_MAX_DIGITS], one[BN_MAX_DIGITS];

    montgomery_mul(bm, b, rr, n, digits, n0inv);
    montgomery_mul(cm, c, rr, n, digits, n0inv);
    montgomery_mul(tm, bm, cm, n, digits, n0inv);
    BN_ASSIGN_DIGIT(one, 1, digits);
    montgomery_mul(a, tm, one, n, digits, n0inv);

    // Clear potentially sensitive information
    memset((uint8_t *)bm, 0, sizeof(bm));
    memset((uint8_t *)cm, 0, sizeof(cm));
    memset((uint8_t *)tm, 0, sizeof(tm));
    memset((uint8_t *)one, 0, sizeof(one));
}

static void montgomery_precompute_constants(bn_t *r, bn_t *rr, bn_t *n, uint32_t digits)
{
    bn_t tmp[2*BN_MAX_DIGITS];

    bn_assign_zero(tmp, 2*digits);
    tmp[digits] = 1;
    bn_mod(r, tmp, digits + 1, n, digits);

    bn_mul(tmp, r, r, digits);
    bn_mod(rr, tmp, 2*digits, n, digits);

    // Clear potentially sensitive information
    memset((uint8_t *)tmp, 0, sizeof(tmp));
}

static void classic_mod_exp(bn_t *a, bn_t *b, bn_t *c, uint32_t cdigits, bn_t *n, uint32_t ndigits)
{
    bn_t result[BN_MAX_DIGITS], base[BN_MAX_DIGITS];
    int bit;

    cdigits = bn_digits(c, cdigits);
    if(cdigits == 0) {
        BN_ASSIGN_DIGIT(a, 1, ndigits);
        return;
    }

    if(bn_cmp(b, n, ndigits) >= 0) {
        bn_mod(base, b, ndigits, n, ndigits);
    } else {
        bn_assign(base, b, ndigits);
    }

    bit = (int)highest_bit_index(c, cdigits);
    bn_assign(result, base, ndigits);
    for(bit--; bit>=0; bit--) {
        bn_mod_mul(result, result, result, n, ndigits);
        if(exponent_bit(c, (uint32_t)bit)) {
            bn_mod_mul(result, result, base, n, ndigits);
        }
    }

    bn_assign(a, result, ndigits);

    // Clear potentially sensitive information
    memset((uint8_t *)result, 0, sizeof(result));
    memset((uint8_t *)base, 0, sizeof(base));
}

#define DEFINE_MONTGOMERY_MUL_FIXED(name, fixed_digits)                                      \
static void name(bn_t *a, bn_t *b, bn_t *c, bn_t *n, bn_t n0inv)                             \
{                                                                                             \
    bn_t t[(fixed_digits)+1], m;                                                              \
    dbn_t acc, carry;                                                                         \
    uint32_t i, j;                                                                            \
                                                                                              \
    memset((uint8_t *)t, 0, sizeof(t));                                                       \
    for(i=0; i<(fixed_digits); i++) {                                                         \
        bn_t bi = b[i];                                                                       \
                                                                                              \
        acc = (dbn_t)bi * c[0] + t[0];                                                        \
        carry = acc >> BN_DIGIT_BITS;                                                         \
        m = (bn_t)acc * n0inv;                                                                \
        acc = (dbn_t)m * n[0] + (bn_t)acc;                                                    \
        carry += acc >> BN_DIGIT_BITS;                                                        \
                                                                                              \
        for(j=1; j<(fixed_digits); j++) {                                                     \
            acc = (dbn_t)bi * c[j] + t[j] + carry;                                            \
            carry = acc >> BN_DIGIT_BITS;                                                     \
            acc = (dbn_t)m * n[j] + (bn_t)acc;                                                \
            t[j-1] = (bn_t)acc;                                                               \
            carry += acc >> BN_DIGIT_BITS;                                                    \
        }                                                                                     \
                                                                                              \
        acc = (dbn_t)t[(fixed_digits)] + carry;                                               \
        t[(fixed_digits)-1] = (bn_t)acc;                                                      \
        t[(fixed_digits)] = (bn_t)(acc >> BN_DIGIT_BITS);                                     \
    }                                                                                         \
                                                                                              \
    bn_assign(a, t, (fixed_digits));                                                          \
    if(t[(fixed_digits)] || bn_cmp(a, n, (fixed_digits)) >= 0) {                              \
        bn_sub(a, a, n, (fixed_digits));                                                      \
    }                                                                                         \
                                                                                              \
    memset((uint8_t *)t, 0, sizeof(t));                                                       \
}

DEFINE_MONTGOMERY_MUL_FIXED(montgomery_mul_64, 64)
DEFINE_MONTGOMERY_MUL_FIXED(montgomery_mul_128, 128)

static void montgomery_mul(bn_t *a, bn_t *b, bn_t *c, bn_t *n, uint32_t digits, bn_t n0inv)
{
    if(digits == 64) {
        montgomery_mul_64(a, b, c, n, n0inv);
    } else if(digits == 128) {
        montgomery_mul_128(a, b, c, n, n0inv);
    } else {
        montgomery_mul_var(a, b, c, n, digits, n0inv);
    }
}

static void montgomery_mul_var(bn_t *a, bn_t *b, bn_t *c, bn_t *n, uint32_t digits, bn_t n0inv)
{
    bn_t t[BN_MAX_DIGITS+1], m;
    dbn_t acc, carry;
    uint32_t i, j;

    bn_assign_zero(t, digits + 1);
    for(i=0; i<digits; i++) {
        bn_t bi = b[i];

        acc = (dbn_t)bi * c[0] + t[0];
        carry = acc >> BN_DIGIT_BITS;
        m = (bn_t)acc * n0inv;
        acc = (dbn_t)m * n[0] + (bn_t)acc;
        carry += acc >> BN_DIGIT_BITS;

        for(j=1; j<digits; j++) {
            acc = (dbn_t)bi * c[j] + t[j] + carry;
            carry = acc >> BN_DIGIT_BITS;
            acc = (dbn_t)m * n[j] + (bn_t)acc;
            t[j-1] = (bn_t)acc;
            carry += acc >> BN_DIGIT_BITS;
        }

        acc = (dbn_t)t[digits] + carry;
        t[digits-1] = (bn_t)acc;
        t[digits] = (bn_t)(acc >> BN_DIGIT_BITS);
    }

    bn_assign(a, t, digits);
    if(t[digits] || bn_cmp(a, n, digits) >= 0) {
        bn_sub(a, a, n, digits);
    }

    // Clear potentially sensitive information
    memset((uint8_t *)t, 0, sizeof(t));
}

static bn_t montgomery_n0inv(bn_t n0)
{
    bn_t x;

    x = 1;
    for(uint32_t i=0; i<5; i++) {
        x *= 2 - n0 * x;
    }

    return (bn_t)(0 - x);
}

static uint32_t exponent_bit(bn_t *a, uint32_t bit)
{
    return (a[bit / BN_DIGIT_BITS] >> (bit % BN_DIGIT_BITS)) & 1U;
}

static uint32_t highest_bit_index(bn_t *a, uint32_t digits)
{
    bn_t top;
    uint32_t bit;

    top = a[digits-1];
    bit = (digits - 1) * BN_DIGIT_BITS;
    while(top >>= 1) {
        bit++;
    }

    return bit;
}
