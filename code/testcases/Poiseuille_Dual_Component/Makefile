test:
	rm -f *.pdf *.h5
	mv ../../constants.h ../../constants.h.store
	cp constants.h ../../constants.h
	cd ../../; make cleandata; make
	mv ../../constants.h.store ../../constants.h
	mv ../../*.h5 .