CFLAGS= -Wall -Werror -fpic -c
MFLAGS= -lsundials_ida -lsundials_nvecserial -lm

%.o:%.c
	gcc -I../include/ $(CFLAGS) $< -o $@

test_freq: helper.o test_freq.c
	gcc -std=gnu99 -I../include/ test_freq.c helper.o -L. -Wl,-rpath=. $(MFLAGS) -o $@

test_sig: helper.o test_sigmoid.c
	gcc -std=gnu99 -I../include/ test_sigmoid.c helper.o -L. -Wl,-rpath=. $(MFLAGS) -o $@
	
test_wien: helper.o test_wiener.c
	gcc -std=gnu99 -I../include/ test_wiener.c helper.o -L. -Wl,-rpath=. $(MFLAGS) -o $@
