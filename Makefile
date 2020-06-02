ALL_EXE = test
all: $(ALL_EXE) 

NE10LOC=/home/pedro/github/Ne10

CC=gcc
CFLAGS=-I. -I$(NE10LOC)/inc -Wall

DEPS = Makefile 

test: test.o cacode.o
# $@: target file name, $^: all prerequisites
	$(CC) -g -o $@ $^ -L$(NE10LOC)/build/modules/ -lNE10 -lm

clean:
	rm -rf *.o $(ALL_EXE) 

#pattern rules
%.o: %.c $(DEPS)
# $<: first prerequisite
	$(CC) -g -c -o $@ $< $(CFLAGS)


