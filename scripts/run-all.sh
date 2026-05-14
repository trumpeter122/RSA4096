#!/usr/bin/env sh
set -eu

TIMEOUT_MULTIPLIER=2
FASTEST_RUN_SECONDS=0

MUL_IMPLS="base karatsuba toom_cook ntt"
MODMUL_IMPLS="base montgomery barrett"
RSA_IMPLS="base square_multiply crt mont_crt"

run_benchmark() {
  log="$1"
  prof="$2"
  started_at="$(date +%s)"

  if [ "$FASTEST_RUN_SECONDS" -gt 0 ]; then
    timeout_seconds=$((FASTEST_RUN_SECONDS * TIMEOUT_MULTIPLIER))
    if timeout "$timeout_seconds" ./release/main >"$log" 2>&1; then
      :
    else
      rm -f "$log" "$prof" gmon.out
      return 0
    fi
  else
    ./release/main >"$log" 2>&1
  fi

  finished_at="$(date +%s)"
  elapsed_seconds=$((finished_at - started_at))
  if [ "$elapsed_seconds" -lt 1 ]; then
    elapsed_seconds=1
  fi
  if [ "$FASTEST_RUN_SECONDS" -eq 0 ] || [ "$elapsed_seconds" -lt "$FASTEST_RUN_SECONDS" ]; then
    FASTEST_RUN_SECONDS="$elapsed_seconds"
  fi

  gprof ./release/main gmon.out >"$prof"
  rm -f gmon.out
}

mkdir -p out

MODMUL_IMPLS=$(echo "$MODMUL_IMPLS" | tr ' ' '\n' | shuf)
MUL_IMPLS=$(echo "$MUL_IMPLS" | tr ' ' '\n' | shuf)
RSA_IMPLS=$(echo "$RSA_IMPLS" | tr ' ' '\n' | shuf)

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

      run_benchmark "$log" "$prof"
    done
  done
done
