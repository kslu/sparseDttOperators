CC=gcc
CFLAGS=-O2

all: speed_tikhonov4x4 speed_ideallowpass32

speed_tikhonov4x4: grfilter.o speed_tikhonov4x4.o
	$(CC) -o $@ $?

speed_tikhonov4x4.o: src/speed_tikhonov4x4.c
	$(CC) -c $<

speed_ideallowpass32: grfilter.o speed_ideallowpass32.o
	$(CC) -o $@ $?

speed_ideallowpass32.o: src/speed_ideallowpass32.c
	$(CC) -c $<

grfilter.o: src/grfilter.c
	@echo "grfilter.o"
	$(CC) -c $<

clean:
		rm -rf *.o speed_tikhonov4x4 speed_ideallowpass32
