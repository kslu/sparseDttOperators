CC=gcc
CFLAGS=-O2
NFLAG=-Wno-nullability-completeness

all: speed_tikhonov16x16 speed_tikhonov64 speed_ideallowpass16x16 speed_ideallowpass64 speed_exp16x16 speed_exp64 

speed_tikhonov4x4: grfilter.o dct.o speed_tikhonov4x4.o
	$(CC) -o $@ $? -lm

speed_tikhonov4x4.o: src/speed_tikhonov4x4.c
	$(CC) $(NFLAG) -c $<

speed_tikhonov8x8: grfilter.o dct.o speed_tikhonov8x8.o
	$(CC) -o $@ $? -lm

speed_tikhonov8x8.o: src/speed_tikhonov8x8.c
	$(CC) $(NFLAG) -c $<

speed_tikhonov16x16: grfilter.o dct.o speed_tikhonov16x16.o
	$(CC) -o $@ $? -lm

speed_tikhonov16x16.o: src/speed_tikhonov16x16.c
	$(CC) $(NFLAG) -c $<

speed_tikhonov32: grfilter.o dct.o speed_tikhonov32.o
	$(CC) -o $@ $? -lm

speed_tikhonov32.o: src/speed_tikhonov32.c
	$(CC) $(NFLAG) -c $<

speed_tikhonov64: grfilter.o dct.o speed_tikhonov64.o
	$(CC) -o $@ $? -lm

speed_tikhonov64.o: src/speed_tikhonov64.c
	$(CC) $(NFLAG) -c $<

speed_tikhonov128: grfilter.o dct.o speed_tikhonov128.o
	$(CC) -o $@ $? -lm

speed_tikhonov128.o: src/speed_tikhonov128.c
	$(CC) $(NFLAG) -c $<

speed_diffusion4x4: grfilter.o dct.o speed_diffusion4x4.o
	$(CC) -o $@ $? -lm

speed_diffusion4x4.o: src/speed_diffusion4x4.c
	$(CC) $(NFLAG) -c $<

speed_diffusion8x8: grfilter.o dct.o speed_diffusion8x8.o
	$(CC) -o $@ $? -lm

speed_diffusion8x8.o: src/speed_diffusion8x8.c
	$(CC) $(NFLAG) -c $<

speed_diffusion16x16: grfilter.o dct.o speed_diffusion16x16.o
	$(CC) -o $@ $? -lm

speed_diffusion16x16.o: src/speed_diffusion16x16.c
	$(CC) $(NFLAG) -c $<

speed_diffusion32: grfilter.o dct.o speed_diffusion32.o
	$(CC) -o $@ $? -lm

speed_diffusion32.o: src/speed_diffusion32.c
	$(CC) $(NFLAG) -c $<

speed_diffusion64: grfilter.o dct.o speed_diffusion64.o
	$(CC) -o $@ $? -lm

speed_diffusion64.o: src/speed_diffusion64.c
	$(CC) $(NFLAG) -c $<

speed_diffusion128: grfilter.o dct.o speed_diffusion128.o
	$(CC) -o $@ $? -lm

speed_diffusion128.o: src/speed_diffusion128.c
	$(CC) $(NFLAG) -c $<

speed_ideallowpass4x4: grfilter.o dct.o speed_ideallowpass4x4.o
	$(CC) -o $@ $? -lm

speed_ideallowpass4x4.o: src/speed_ideallowpass4x4.c
	$(CC) $(NFLAG) -c $<

speed_ideallowpass8x8: grfilter.o dct.o speed_ideallowpass8x8.o
	$(CC) -o $@ $? -lm

speed_ideallowpass8x8.o: src/speed_ideallowpass8x8.c
	$(CC) $(NFLAG) -c $<

speed_ideallowpass16x16: grfilter.o dct.o speed_ideallowpass16x16.o
	$(CC) -o $@ $? -lm

speed_ideallowpass16x16.o: src/speed_ideallowpass16x16.c
	$(CC) $(NFLAG) -c $<

speed_ideallowpass32: grfilter.o dct.o speed_ideallowpass32.o
	$(CC) -o $@ $? -lm

speed_ideallowpass32.o: src/speed_ideallowpass32.c
	$(CC) $(NFLAG) -c $<

speed_ideallowpass64: grfilter.o dct.o speed_ideallowpass64.o
	$(CC) -o $@ $? -lm

speed_ideallowpass64.o: src/speed_ideallowpass64.c
	$(CC) $(NFLAG) -c $<

speed_ideallowpass128: grfilter.o dct.o speed_ideallowpass128.o
	$(CC) -o $@ $? -lm

speed_ideallowpass128.o: src/speed_ideallowpass128.c
	$(CC) $(NFLAG) -c $<

speed_exp4x4: grfilter.o dct.o speed_exp4x4.o
	$(CC) -o $@ $? -lm

speed_exp4x4.o: src/speed_exp4x4.c
	$(CC) $(NFLAG) -c $<

speed_exp8x8: grfilter.o dct.o speed_exp8x8.o
	$(CC) -o $@ $? -lm

speed_exp8x8.o: src/speed_exp8x8.c
	$(CC) $(NFLAG) -c $<

speed_exp16x16: grfilter.o dct.o speed_exp16x16.o
	$(CC) -o $@ $? -lm

speed_exp16x16.o: src/speed_exp16x16.c
	$(CC) $(NFLAG) -c $<

speed_exp32: grfilter.o dct.o speed_exp32.o
	$(CC) -o $@ $? -lm

speed_exp32.o: src/speed_exp32.c
	$(CC) $(NFLAG) -c $<

speed_exp64: grfilter.o dct.o speed_exp64.o
	$(CC) -o $@ $? -lm

speed_exp64.o: src/speed_exp64.c
	$(CC) $(NFLAG) -c $<

speed_exp128: grfilter.o dct.o speed_exp128.o
	$(CC) -o $@ $? -lm

speed_exp128.o: src/speed_exp128.c
	$(CC) $(NFLAG) -c $<

grfilter.o: src/grfilter.c
	@echo "grfilter.o"
	$(CC) $(NFLAG) -c $<

dct.o: src/dct.c
	@echo "dct.o"
	$(CC) $(NFLAG) -c $<

test_dct: dct.o test_dct.o
	$(CC) -o $@ $? -lm

test_dct.o: test/test_dct.c
	$(CC) $(NFLAG) -c $<

test_arma: grfilter.o dct.o test_arma.o
	$(CC) -o $@ $? -lm

test_arma.o: test/test_arma.c
	$(CC) $(NFLAG) -c $<

test_cheby: grfilter.o dct.o test_cheby.o
	$(CC) -o $@ $? -lm

test_cheby.o: test/test_cheby.c
	$(CC) $(NFLAG) -c $<

clean:
		rm -rf *.o speed_tikhonov4x4 speed_tikhonov8x8 speed_tikhonov16x16 speed_tikhonov32 speed_tikhonov64 speed_tikhonov128 speed_ideallowpass4x4 speed_ideallowpass8x8 speed_ideallowpass16x16 speed_ideallowpass32 speed_ideallowpass64 speed_ideallowpass128 speed_exp4x4 speed_exp8x8 speed_exp16x16 speed_exp32 speed_exp64 speed_exp128 speed_diffusion4x4 speed_diffusion8x8 speed_diffusion16x16 speed_diffusion32 speed_diffusion64 speed_diffusion128 test_dct test_cheby test_arma
