CC	= gcc   
CFLAGS    = -fopenmp
LIBS	  = -lm

par:	par.c
	$(CC) $(CFLAGS) -o par par.c $(LIBS)
par2:	par2.c
	$(CC) $(CFLAGS) -o par2 par2.c $(LIBS)
qr_par:	qr_par.c
	$(CC) $(CFLAGS) -o qr_par qr_par.c $(LIBS)
qr_par2:	qr_par2.c
	$(CC) $(CFLAGS) -o qr_par2 qr_par2.c $(LIBS)
test:	test.c
	$(CC) $(CFLAGS) -o test test.c $(LIBS)

clean: 
	rm -f *.o par par2 test
