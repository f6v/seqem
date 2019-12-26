FLAGS = -fno-exceptions -c -w -O3 -g -Wall -DLINUX -D_LIBC_REENTRANT -D_THREAD_SAFE -D_REENTRANT
COMPILER = g++

.SUFFIXES: .cpp 

.cpp.o:
	$(COMPILER) $(FLAGS) $*.cpp

SRCS=	main.cpp		\
	analyze.cpp		\
	ccds.cpp		\
	combinatorics.cpp	\
	sort.cpp		\
	invmat.cpp		\
	binsearch.cpp		\
	read_data.cpp		\
	read_param.cpp		\
	sorter_t.cpp		\
	searcher_t.cpp		\
	ccds_i.cpp

CSRCS=	dcdflib.c

OBJS=$(SRCS:.cpp=.o)

OBJSC=$(CSRCS:.c=.o)

all: seqem

seqem:	$(OBJS) $(OBJSC)
		$(COMPILER) -o seqem $(OBJS) $(OBJSC) -lpthread -lm

clean:
		rm $(OBJS) $(OBJSC)
