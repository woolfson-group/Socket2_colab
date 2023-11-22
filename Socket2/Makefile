# simple Makefile for SOCKET

TOP             = .
WORKDIR         = $(TOP)

INC             = /usr/include
WORKINC         = $(WORKDIR)

#INCLUDES        = -I$(INC) -I$(WORKINC)

LIBS            = -lm
CC		= cc
CFLAGS		= $(INCLUDES) $(LIBS)
#WARNINGS	= -woff all

OBJECTS	= socket.o daisies.o duplicat.o geometry.o knobs.o order.o read.o register.o statchar.o

HEADERS = socket.h preproc.h global.h prototyp.h aminoa1.h

socket : $(OBJECTS) $(HEADERS)
	$(CC) -o socket2 $(OBJECTS) $(CFLAGS) $(WARNINGS)

socket.o : socket.c $(HEADERS) aminoa2.h

daisies.o : daisies.c $(HEADERS)

duplicat.o : duplicat.c $(HEADERS)

geometry.o : geometry.c $(HEADERS)

knobs.o : knobs.c $(HEADERS)

order.o : order.c $(HEADERS)

read.o : read.c $(HEADERS)

register.o : register.c $(HEADERS)

statchar.o : statchar.c $(HEADERS)

.PHONY : clean
clean :
	rm socket2 $(OBJECTS)
