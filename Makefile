CC=g++
CFLAGS=-O3
LDFLAGS=-static

prog=trimMSA rmRedundantSeq realignMSA fasta2aln fastaCov fastNf calNf cleanFastaHeader cleanFastaBody AlnAaProb RemoveNonQueryPosition fasta2pfam calNf_ly a3m2msa unaligna3m fasta2crc64 fastaOneLine trimBlastN afa2blasttab fixAlnX

all: ${prog}

trimMSA: trimMSA.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

fixAlnX: fixAlnX.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

afa2blasttab: afa2blasttab.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

trimBlastN: trimBlastN.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

fastaNA: fastaNA.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

fastaOneLine: fastaOneLine.cpp
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

cleanFastaBody: cleanFastaBody.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

AlnAaProb: AlnAaProb.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

RemoveNonQueryPosition: RemoveNonQueryPosition.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

calNf_ly: calNf_ly.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

a3m2msa: a3m2msa.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

unaligna3m: unaligna3m.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

fasta2crc64: fasta2crc64.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

install: ${prog}
	cp ${prog} ../bin

clean:
	rm ${prog}
