LIBPICNIC=libpicnic3-l1.a

SOURCES=$(wildcard *.c) ../../../mupq/common/fips202.c ../../../mupq/common/keccakf1600.c

#CFLAGS=-O3 -Wall -Wconversion -Wextra -Wpedantic -Wvla -Werror -Wmissing-prototypes -Wredundant-decls -std=c99 -fvisibility=default
CFLAGS=-Og -g -Wall -Wextra -Wvla -Werror -Wredundant-decls -std=c99 -I../../../mupq/common

all: $(LIBPICNIC)

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

$(LIBPICNIC): $(SOURCES:.c=.o) $(SOURCES:.s=.o)
	ar rcs $@ $^

kats: $(LIBPICNIC) tests/kats_test.c
	$(CC) $(CFLAGS) -o kats.out tests/kats_test.c ${LIBPICNIC}

bench: $(OBJECTS) tests/bench.c tests/bench_timing.c tests/bench_utils.c
	$(CC) $(CFLAGS) -D_GNU_SOURCE -o bench.out tests/bench.c tests/bench_timing.c tests/bench_utils.c $(LIBPICNIC)
clean:
	rm -f $(wildcard *.o) $(wildcard sha3/*.o) $(LIBPICNIC) $(EXECUTABLE_APITEST) $(EXECUTABLE_EXAMPLE) $(EXECUTABLE_NISTKAT)
