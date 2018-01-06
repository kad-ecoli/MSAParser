CC=g++
CFLAGS=-O3
LDFLAGS=-static

prog=trimMSA rmRedundantSeq realignMSA fasta2aln fastaCov

all: ${prog}

trimMSA: trimMSA.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

rmRedundantSeq: rmRedundantSeq.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

realignMSA: realignMSA.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

fasta2aln: fasta2aln.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

fastaCov: fastaCov.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

install: ${prog}
	cp ${prog} ../bin

clean:
	rm ${prog}
