MACHINEFLAGS = 
CFLAGS = $(MACHINEFLAGS)
LFLAGS = -static -lm $(MACHINEFLAGS)
OBJLIST = dots.o abin.o readPDBrecs.o geom3d.o utility.o select.o \
	 parse.o atomprops.o stdconntable.o autobondrot.o hybrid_36_c.o

.c.o:
	cc -c $*.c $(CFLAGS)

probe: probe.o $(OBJLIST)
	cc -o $@ probe.o $(OBJLIST) $(LFLAGS)

clean:
	@rm -f *.o *.ckp

install: probe
	mv probe /local/bin

installtest: probe
	mv probe /local/bin/probeTest

# DO NOT DELETE THIS LINE -- make depend uses it
abin.o:	./abin.h	./geom3d.h	./utility.h	abin.c
atomprops.o:	./atomprops.h	atomprops.c
autobondrot.o:	./abin.h	./autobondrot.h	./geom3d.h	./readPDBrecs.h \
		./utility.h	autobondrot.c
dots.o:	./dots.h	./geom3d.h	dots.c
geom3d.o:	./geom3d.h	geom3d.c
parse.o:	./parse.h	./utility.h	parse.c
probe.o:	./abin.h	./atomprops.h	./autobondrot.h	./dots.h	./geom3d.h \
		./parse.h	./probe.h	./readPDBrecs.h	./select.h \
		./stdconntable.h	./utility.h	probe.c
readPDBrecs.o:	./geom3d.h	./readPDBrecs.h	./utility.h	./hybrid_36_c.h readPDBrecs.c
select.o:	./abin.h	./atomprops.h	./geom3d.h	./parse.h	./select.h \
		./stdconntable.h	./utility.h	select.c
stdconntable.o:	./stdconntable.h	stdconntable.c
utility.o:	utility.c
hybrid_36_c.o:  ./hybrid_36_c.h	hybrid_36_c.c
# DO NOT DELETE THIS 2nd LINE -- make depend uses it
