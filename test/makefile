CFLAGS = -Wall -lm

main: main.c scan.c scan.h
	gcc -o $@ main.c scan.c  $(CFLAGS)
	
clean:
	-rm -f *.o *~ *.gch main
