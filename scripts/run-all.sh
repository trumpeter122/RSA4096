#!/usr/bin/env sh
set -eu

MUL_IMPLS="base base_mt karatsuba karatsuba_mt toom_cook ntt"
MODMUL_IMPLS="base montgomery montgomery_mt barrett"
RSA_IMPLS="base square_multiply crt crt_mt mont_crt"

mkdir -p out

for mul in $MUL_IMPLS; do
  for modmul in $MODMUL_IMPLS; do
    for rsa in $RSA_IMPLS; do
      name="mul-${mul}_modmul-${modmul}_rsa-${rsa}"
      log="out/${name}.log"
      prof="out/${name}.prof"

      if [ -f "$log" ] && [ -f "$prof" ]; then
        echo "==> $name exists, skipping"
        continue
      fi

      echo "==> $name"

      make clean >/dev/null
      make MUL="$mul" MODMUL="$modmul" RSA="$rsa" >/dev/null

      make clean >/dev/null
      make MUL="$mul" MODMUL="$modmul" RSA="$rsa" 'RELCFLAGS=-O3 -DNDEBUG -pg' >/dev/null

      ./release/main >"$log" 2>&1
      gprof ./release/main gmon.out >"$prof"
      rm gmon.out
    done
  done
done
