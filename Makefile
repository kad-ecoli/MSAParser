CC=g++
CFLAGS=-O3
LDFLAGS=-static

prog=trimMSA rmRedundantSeq realignMSA fasta2aln fastaCov fastNf calNf cleanFastaHeader AlnAaProb RemoveNonQueryPosition fasta2pfam calNf_ly

all: ${prog}

trimMSA: trimMSA.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

rmRedundantSeq: rmRedundantSeq.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

realignMSA: realignMSA.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

fasta2aln: fasta2aln.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

fasta2pfam: fasta2pfam.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

fastaCov: fastaCov.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

calNf: calNf.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

fastNf: fastNf.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

cleanFastaHeader: cleanFastaHeader.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

AlnAaProb: AlnAaProb.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

RemoveNonQueryPosition: RemoveNonQueryPosition.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

calNf_ly: calNf_ly.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

install: ${prog}
	cp ${prog} ../bin

clean:
	rm ${prog}
