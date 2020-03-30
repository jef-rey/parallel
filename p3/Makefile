CC = gcc
CFLAGS = -g -Wall -W -Werror -pedantic -std=c99
LDLIBS = -lm
TARGET = tmin
all: tmin ptrap

tmin:
	    $(CC) $(CFLAGS) -o $(TARGET) $(TARGET).c $(LDLIBS)

ptrap:
	mpicc -g -Wall -std=c99 -o ptrap ptrap.c $(LDLIBS)

clean:
	    -rm -f tmin ptrap 
