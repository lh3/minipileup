CC ?=			gcc
CXX ?=		g++
CFLAGS ?=		-g -Wall -O3
CXXFLAGS ?=	$(CFLAGS)
CPPFLAGS=
INCLUDES=
OBJS ?=		razf.o bgzf.o hts.o sam.o faidx.o bedidx.o
PROG ?=		minipileup
LIBS ?=		-lpthread -lz -lm

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address -ldl
endif

.SUFFIXES:.c .cpp .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

.cpp.o:
		$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

minipileup:$(OBJS) pileup.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr *.o a.out $(PROG) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c *.cpp)

# DO NOT DELETE

bedidx.o: ksort.h kseq.h khash.h
bgzf.o: bgzf.h
faidx.o: faidx.h khash.h razf.h
hts.o: bgzf.h hts.h kseq.h khash.h ksort.h
pileup.o: sam.h bgzf.h hts.h faidx.h ksort.h ketopt.h
razf.o: razf.h
sam.o: sam.h bgzf.h hts.h khash.h kseq.h kstring.h
