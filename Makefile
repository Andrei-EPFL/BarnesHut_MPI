OPTIM+=-O1 -march=native
CXX=g++
CC=g++
LD=${CXX}
CXXFLAGS+=-Wall -Wextra -Werror -pedantic -std=c++11 $(OPTIM)
LDFLAGS+=-lm $(CXXFLAGS)
OBJS=main.o update_node.o dynamics.o

all: main

main: $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)

clean:
	rm -f main *.o *~

