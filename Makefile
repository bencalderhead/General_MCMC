CC = gcc
CFLAGS = -ggdb -march=native -pipe -std=c99 -pedantic -Wall -Wextra -Wconversion
CPPFLAGS = -I.

LDLIBS = -lm

RM = rm -f

TARGET = libgmcmc.so

OBJS = src/error.o \
       src/rng/rng.o src/rng/mt64.o src/rng/dsfmt.o \
       src/priors/prior.o src/priors/uniform.o src/priors/normal.o src/priors/lognormal.o src/priors/gamma.o

TEST_OBJS = test/rng/rng.o test/rng/mt64.o test/rng/dsfmt.o \
            test/priors/prior.o test/priors/uniform.o test/priors/normal.o test/priors/lognormal.o src/priors/gamma.o
TESTS = rng mt64 dsfmt prior gamma normal lognormal uniform

.PHONY: all test clean

all: $(TARGET) $(TESTS)

test: $(TESTS)

clean:
	$(RM) $(TARGET) $(TESTS) $(OBJS) $(TEST_OBJS)

# Target library (default target)
$(TARGET): CFLAGS += -fPIC
$(TARGET): $(OBJS)

# Tests
$(TESTS): LDFLAGS += -L.
$(TESTS): LDLIBS += -lgmcmc -lcunit -lcurses
$(TESTS):
	$(CC) $(^) -o $(@) $(LDFLAGS) $(LDLIBS)

rng: test/rng/rng.o
mt64: test/rng/mt64.o
dsfmt: test/rng/dsfmt.o
prior: test/priors/prior.o
uniform: test/priors/uniform.o
normal: test/priors/normal.o
lognormal: test/priors/lognormal.o
gamma: test/priors/gamma.o

%.so:
	$(CC) $(^) -shared -o $(@) $(LDFLAGS) $(LDLIBS)
