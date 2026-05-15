#!/usr/bin/env sh
set -eu

TIMEOUT_MULTIPLIER=2
FASTEST_RUN_SECONDS=0

MUL_IMPLS="base karatsuba toom_cook ntt"
MODMUL_IMPLS="base montgomery barrett"
RSA_IMPLS="base square_multiply crt mont_crt"

MUL_IMPLS=$(echo "$MUL_IMPLS" | tr ' ' '\n' | shuf)
MODMUL_IMPLS=$(echo "$MODMUL_IMPLS" | tr ' ' '\n' | shuf)
RSA_IMPLS=$(echo "$RSA_IMPLS" | tr ' ' '\n' | shuf)

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

run_one() {
  mul="$1"
  modmul="$2"
  rsa="$3"

  name="mul-${mul}_modmul-${modmul}_rsa-${rsa}"
  log="out/${name}.log"
  prof="out/${name}.prof"

  current_round=$((current_round + 1))

  echo "==> $current_round/$total_rounds | Current record: ${FASTEST_RUN_SECONDS}s"

  if [ -f "$log" ] && [ -f "$prof" ]; then
    echo "$name exists, skipping"
    return
  fi

  echo "$name running"

  make clean >/dev/null
  make MUL="$mul" MODMUL="$modmul" RSA="$rsa" >/dev/null

  make clean >/dev/null
  make MUL="$mul" MODMUL="$modmul" RSA="$rsa" 'RELCFLAGS=-O3 -DNDEBUG -pg' >/dev/null

  run_benchmark "$log" "$prof"
}

mkdir -p out

total_rounds=$(($(echo "$MUL_IMPLS" | wc -l) * $(echo "$MODMUL_IMPLS" | wc -l) * $(echo "$RSA_IMPLS" | wc -l)))
current_round=0

run_one base base base

for mul in $MUL_IMPLS; do
  for modmul in $MODMUL_IMPLS; do
    for rsa in $RSA_IMPLS; do
      if [ "$mul" != "base" ] || [ "$modmul" != "base" ] || [ "$rsa" != "base" ]; then
        run_one "$mul" "$modmul" "$rsa"
      fi
    done
  done
done
