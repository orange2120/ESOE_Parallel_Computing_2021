life:
	mpicc -Wall --std=c11 -lm 6_13_Life.c -o $@

life_seq:
	mpicc -Wall --std=c11 -lm 6_13_Life_seq.c -o $@

gen:
	g++ -Wall gen.cpp -o $@

clean:
	rm life life_seq gen