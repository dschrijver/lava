COMPILER = gcc
OPT_FLAGS = -std=c17 -O3
DEBUG_FLAGS = -Wall -Wextra
# Static linking with hdf5, link with compression library too. Order matters, hdf5 must be linked first!
LIBRARIES = -I/usr/local/hdf5/include -L/usr/local/hdf5/lib -l:libhdf5.a -lm -lz -lraylib
# LIBRARIES = -I/usr/local/hdf5/include -L/usr/local/hdf5/lib -Wl,-rpath /usr/local/hdf5/lib -l:libhdf5.so -lm -lz

obstacle: clean obstacle.out
	./obstacle.out
	python3 obstacle_read.py

poiseuille: clean poiseuille.out
	./poiseuille.out
	python3 poiseuille_read.py

poiseuille-sc: clean poiseuille-sc.out
	./poiseuille-sc.out

rayleigh: clean rayleigh.out
	./rayleigh.out

stefan: clean stefan.out
	./stefan.out

channel: clean channel.out
	./channel.out

%.out: %.c
	$(COMPILER) $< -o $@ $(OPT_FLAGS) $(DEBUG_FLAGS) $(LIBRARIES)

clean:
	rm -f *.out

cleandata:
	rm -f *.h5
