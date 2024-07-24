cc = g++  -Wextra -g -fsanitize=address -fsanitize-recover=address


DEPS = BIN/main.o BIN/KS.o BIN/lapack.o BIN/solutore.o  BIN/bisezione.o BIN/integrali.o BIN/trapezi.o BIN/dati.o
LDFLAGS = -llapack

main.exe: ${DEPS}
	${cc} -o $@ $^ ${LDFLAGS}
BIN/main.o: main.cpp
	${cc} -o $@ -c $<
BIN/lapack.o: lapack.cpp LIB/lapack.h
	${cc} -o $@ -c $<
BIN/solutore.o: solutore.cpp LIB/solutore.h
	${cc} -o $@ -c $<
BIN/bisezione.o: bisezione.cpp LIB/bisezione.h
	${cc} -o $@ -c $<
BIN/integrali.o: integrali.cpp LIB/integrali.h
	${cc} -o $@ -c $<
BIN/trapezi.o: trapezi.cpp LIB/trapezi.h
	${cc} -o $@ -c $<
BIN/KS.o: KS.cpp LIB/KS.h
	${cc} -o $@ -c $<
BIN/dati.o: dati.cpp LIB/dati.h
	${cc} -o $@ -c $<

clean:
	rm BIN/*.o