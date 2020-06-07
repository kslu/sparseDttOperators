CC=gcc
CFLAGS=-O2

all: speed_tikhonov4x4 speed_tikhonov8x8 speed_tikhonov32 speed_tikhonov128 speed_ideallowpass4x4 speed_ideallowpass8x8 speed_ideallowpass32 speed_ideallowpass128

speed_tikhonov4x4: grfilter.o speed_tikhonov4x4.o
	$(CC) -o $@ $?

speed_tikhonov4x4.o: src/speed_tikhonov4x4.c
	$(CC) -c $<

speed_tikhonov8x8: grfilter.o speed_tikhonov8x8.o
	$(CC) -o $@ $?

speed_tikhonov8x8.o: src/speed_tikhonov8x8.c
	$(CC) -c $<

speed_tikhonov32: grfilter.o speed_tikhonov32.o
	$(CC) -o $@ $?

speed_tikhonov32.o: src/speed_tikhonov32.c
	$(CC) -c $<

speed_tikhonov128: grfilter.o speed_tikhonov128.o
	$(CC) -o $@ $?

speed_tikhonov128.o: src/speed_tikhonov128.c
	$(CC) -c $<

speed_ideallowpass4x4: grfilter.o speed_ideallowpass4x4.o
	$(CC) -o $@ $?

speed_ideallowpass4x4.o: src/speed_ideallowpass4x4.c
	$(CC) -c $<

speed_ideallowpass8x8: grfilter.o speed_ideallowpass8x8.o
	$(CC) -o $@ $?

speed_ideallowpass8x8.o: src/speed_ideallowpass8x8.c
	$(CC) -c $<

speed_ideallowpass32: grfilter.o speed_ideallowpass32.o
	$(CC) -o $@ $?

speed_ideallowpass32.o: src/speed_ideallowpass32.c
	$(CC) -c $<

speed_ideallowpass128: grfilter.o speed_ideallowpass128.o
	$(CC) -o $@ $?

speed_ideallowpass128.o: src/speed_ideallowpass128.c
	$(CC) -c $<

grfilter.o: src/grfilter.c
	@echo "grfilter.o"
	$(CC) -c $<

clean:
		rm -rf *.o speed_tikhonov4x4 speed_tikhonov8x8 speed_tikhonov32 speed_tikhonov128 speed_ideallowpass4x4 speed_ideallowpass8x8 speed_ideallowpass32 speed_ideallowpass128
