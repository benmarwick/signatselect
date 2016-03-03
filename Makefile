CD = ocamldoc
CO = ocamlopt
CC = ocamlc

DEP = str.cmxa unix.cmxa nums.cmxa hash2.cmx aux.cmx io.cmx

CKEYS = -o
LIBKEYS = -c
DOCKEYS = -html -colorize-code -d

GSLDIR = -I ~/.opam/system/lib/gsl/
GSLO = bigarray.cmxa gsl.cmxa
GSLC = bigarray.cma gsl.cma


hash2.cmo hash2.cmx: hash2.ml hash2.mli
	$(CC) $(LIBKEYS) hash2.mli
	$(CC) $(LIBKEYS) hash2.ml
	$(CO) $(LIBKEYS) hash2.mli
	$(CO) $(LIBKEYS) hash2.ml

aux.cmo aux.cmx: aux.ml aux.mli hash2.cmo hash2.cmx
	$(CC) $(LIBKEYS) unix.cma nums.cma hash2.cmo aux.mli
	$(CC) $(LIBKEYS) unix.cma nums.cma hash2.cmo aux.ml
	$(CO) $(LIBKEYS) unix.cmxa nums.cmxa hash2.cmx aux.mli
	$(CO) $(LIBKEYS) unix.cmxa nums.cmxa hash2.cmx aux.ml

io.cmo io.cmx: io.mli io.ml
	$(CC) $(LIBKEYS)  io.mli
	$(CC) $(LIBKEYS) io.ml
	$(CO) $(LIBKEYS)  io.mli
	$(CO) $(LIBKEYS) io.ml

tsinfer : tsinfer.ml
	$(CO) $(CKEYS) $@ $(GSLDIR) $(GSLO) $(DEP) $?

all : 
	make hash2.cmo
	make aux.cmo
	make io.cmo
	make tsinfer

clean :
	rm -f *.cmi
	rm -f *.o
	rm -f *~

clean_all :
	rm -f *.cmi
	rm -f *.cmo
	rm -f *.o
	rm -f *.cmx
	rm -f *~
