all:
	odin build . -out:mp_vec_sum.exe

opti:
	odin build . -out:mp_vec_sum_opti.exe -o:speed -no-bounds-check

clean:
	rm -f ./mp_vec_sum.exe \
	./mp_vec_sum_opti.exe

run:
	./mp_vec_sum.exe

run_opti:
	./mp_vec_sum_opti.exe
