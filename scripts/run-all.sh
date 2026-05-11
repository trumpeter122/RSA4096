#!/usr/bin/env sh
set -eu

MUL_IMPLS="base karatsuba"
MODMUL_IMPLS="base montgomery"
RSA_IMPLS="base square_multiply crt"

mkdir -p out

for mul in $MUL_IMPLS; do
    for modmul in $MODMUL_IMPLS; do
        for rsa in $RSA_IMPLS; do
            name="mul-${mul}_modmul-${modmul}_rsa-${rsa}"
            echo "==> $name"
            make clean >/dev/null
            make MUL="$mul" MODMUL="$modmul" RSA="$rsa" >/dev/null
            ./release/main >"out/${name}.log" 2>&1
        done
    done
done
