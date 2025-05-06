libqsld.a: linalg/vector.d linalg/matrix.d quantum/qc.d
	dmd -lib -O -of=libqsld.a linalg/vector.d linalg/matrix.d quantum/qc.d

main: main.d libqsld.a
	dmd -O -L-L. -L="-lqsld" -of=main main.d

.PHONY: all

all: libqsld.a main

clean: 
	rm main main.o libqsld.a
