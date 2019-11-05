OPTIM+=-O1 -march=native
CXX=mpicxx
LD=${CXX}
CXXFLAGS+=-Wall -Wextra -Werror  -std=c++11 $(OPTIM)
LDFLAGS+=-lm $(CXXFLAGS)

OBJS=main.o update_node.o dynamics.o

all: main

main: $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)

clean:
	rm -f main *.o *~
