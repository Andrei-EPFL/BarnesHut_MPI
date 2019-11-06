OPTIM+=-O1 -march=native
CXX=mpicxx
LD=${CXX}
CXXFLAGS+=-Wall -Wextra -Werror -pedantic -std=c++11 $(OPTIM)
LDFLAGS+=-lm $(CXXFLAGS)

OBJS=main.o update_node.o dynamics.o node_func.o

all: main

main: $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)

clean:
	rm -f main *.o *~
