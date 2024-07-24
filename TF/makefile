cc = g++

main.exe: ./BIN/main.o ./BIN/dati.o ./BIN/TF.o ./BIN/solutore.o ./BIN/bisezione.o ./BIN/integrali.o ./BIN/trapezi.o
	${cc} -o $@ $^
./BIN/main.o: main.cpp ./LIB/TF.h
	${cc} -o $@ -c $<
./BIN/solutore.o: solutore.cpp ./LIB/solutore.h
	${cc} -o $@ -c $<
./BIN/bisezione.o: bisezione.cpp ./LIB/solutore.h ./LIB/bisezione.h
	${cc} -o $@ -c $<
./BIN/integrali.o: integrali.cpp ./LIB/integrali.h
	${cc} -o $@ -c $<
./BIN/trapezi.o: trapezi.cpp ./LIB/integrali.h ./LIB/trapezi.h
	${cc} -o $@ -c $<
./BIN/dati.o: dati.cpp ./LIB/dati.h ./LIB/trapezi.h ./LIB/bisezione.h
	${CC} -o $@ -c $<
./BIN/TF.o: TF.cpp ./LIB/dati.h 
	${CC} -o $@ -c $<

clean:
	rm ./BIN/*.o