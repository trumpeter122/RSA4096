#
# Compiler flags
#
CC     = gcc
CFLAGS = -Wall -Wextra -I.
THREADFLAGS = -pthread

MUL    ?= base_mt
MODMUL ?= montgomery_mt
RSA    ?= mont_crt

MUL_IMPLS    = base base_mt karatsuba karatsuba_mt toom_cook ntt
MODMUL_IMPLS = base montgomery montgomery_mt barrett
RSA_IMPLS    = base square_multiply crt crt_mt mont_crt

#
# Project files
#

SRCS = main.c rsa.c bignum.c mul/$(MUL).c modmul/$(MODMUL).c rsa/$(RSA).c
OBJS = $(SRCS:.c=.o)
EXE  = main

#
# Debug build settings
#
DBGDIR = debug
DBGEXE = $(DBGDIR)/$(EXE)
DBGOBJS = $(addprefix $(DBGDIR)/, $(OBJS))
DBGCFLAGS = -g -O0 -DDEBUG

#
# Release build settings
#
RELDIR = release
RELEXE = $(RELDIR)/$(EXE)
RELOBJS = $(addprefix $(RELDIR)/, $(OBJS))
RELCFLAGS = -O3 -DNDEBUG

.PHONY: all check-modules clean debug prep release remake

# Default build
all: check-modules prep release

check-modules:
	@test -f mul/$(MUL).c || (echo "Unknown MUL implementation: $(MUL)"; exit 1)
	@test -f modmul/$(MODMUL).c || (echo "Unknown MODMUL implementation: $(MODMUL)"; exit 1)
	@test -f rsa/$(RSA).c || (echo "Unknown RSA implementation: $(RSA)"; exit 1)

#
# Debug rules
#
debug: check-modules prep $(DBGEXE)

$(DBGEXE): $(DBGOBJS)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(DBGCFLAGS) $(THREADFLAGS) -o $(DBGEXE) $^

$(DBGDIR)/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) -c $(CFLAGS) $(DBGCFLAGS) $(THREADFLAGS) -o $@ $<

#
# Release rules
#
release: check-modules prep $(RELEXE)

$(RELEXE): $(RELOBJS)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(RELCFLAGS) $(THREADFLAGS) -o $(RELEXE) $^

$(RELDIR)/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) -c $(CFLAGS) $(RELCFLAGS) $(THREADFLAGS) -o $@ $<

#
# Other rules
#
prep:
	@mkdir -p $(DBGDIR) $(RELDIR)

remake: clean all

clean:
	rm -f $(RELEXE) $(DBGEXE) \
		$(RELDIR)/*.o $(RELDIR)/mul/*.o $(RELDIR)/modmul/*.o $(RELDIR)/rsa/*.o \
		$(DBGDIR)/*.o $(DBGDIR)/mul/*.o $(DBGDIR)/modmul/*.o $(DBGDIR)/rsa/*.o
