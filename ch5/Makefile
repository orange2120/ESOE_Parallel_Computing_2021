CC = mpicc
CFLAGS = -Wall --std=c11 -lm -O3
LDFLAGS = -lgmp 

hps:
	$(CC) $(CFLAGS) $(LDFLAGS) 5_11_harmonic_progression_sum.c -o $@

hps_gmp:
	$(CC) $(CFLAGS) $(LDFLAGS) 5_11_harmonic_progression_sum_gmp.c -o $@

clean:
	rm hps hps_gmp