# RSA4096

RSA4096 implementation in C with a modular optimization structure. The project
keeps the original top-level RSA interface while allowing multiplication,
modular multiplication, and RSA computation strategies to be selected at build
time.

## Project Structure

```text
.
├── main.c                  # Test / benchmark entry point
├── keys.h                  # RSA test keys
├── rsa.h, rsa.c            # Public RSA API, padding, block handling
├── bignum.h, bignum.c      # Big number helpers and module dispatch
├── Makefile                # Build configuration and module selection
├── mul/                    # Multiplication implementations
│   ├── base.c
│   ├── ...
│   └── mul.h
├── modmul/                 # Modular multiplication implementations
│   ├── base.c
│   ├── ...
│   └── modmul.h
├── rsa/                    # RSA computation implementations
│   ├── base.c
│   ├── ...
│   └── compute.h
└── scripts/                # Benchmark and analysis helpers
    ├── run-all.sh
    └── analyze.py
```

The modular boundaries are:

- `mul`: implements `mul_bn_mul`, used by `bn_mul`.
- `modmul`: implements `modmul_bn_mod_mul`, used by `bn_mod_mul`.
- `rsa`: implements modular exponentiation and RSA private/public compute
  routines declared in `rsa/compute.h`.

Each module has a `base` implementation for the original behavior plus
optimized variants. `main.c` and `keys.h` are intentionally left unchanged.

## Build

Default release build:

```sh
make
```

The default optimized combination is the one with the best performance:

```text
MUL=base_mt MODMUL=montgomery_mt RSA=mont_crt
```

Build the original/base combination:

```sh
make clean
make MUL=base MODMUL=base RSA=base
```

Build with a specific combination:

```sh
make clean
make MUL=karatsuba_mt MODMUL=montgomery_mt RSA=crt_mt
```

Debug build:

```sh
make clean
make debug
gdb ./debug/main
```

Run the release binary:

```sh
./release/main > out.log 2>&1
```

## Available Implementations

`MUL` choices:

- `base`
- `base_mt`
- `karatsuba`
- `karatsuba_mt`

`MODMUL` choices:

- `base`
- `montgomery`
- `montgomery_mt`

`RSA` choices:

- `base`
- `square_multiply`
- `crt`
- `crt_mt`
- `mont_crt`

## Scripts

Run all module combinations and collect logs/profiles under `out/`:

```sh
scripts/run-all.sh
```

Analyze the generated results:

```sh
scripts/analyze.py out
```

Optionally write CSV output:

```sh
scripts/analyze.py out --csv results.csv
```

## Notes

- The all-base build is intended to match the original implementation.
- Optimizations are isolated inside `mul/`, `modmul/`, and `rsa/`.
