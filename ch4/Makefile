

CC = mpicc
CFLAGS = -Wall --std=c11 -lm


prime:
	$(CC) $(CFLAGS) 4_8_prime.c -o $@

prime_gap:
	$(CC) $(CFLAGS) 4_9_lgap_prime.c -o $@

pi_calc:
	$(CC) $(CFLAGS) 4_12_pi_calc.c -o $@

clean:
	rm prime prime_gap pi_calc
