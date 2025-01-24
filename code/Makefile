COMPILER = gcc
OPT_FLAGS = -std=c17 -O3
DEBUG_FLAGS = -Wall -Wextra -Wno-unused-variable
LIBRARIES = -I/usr/local/hdf5/include -L/usr/local/hdf5/lib -l:libhdf5.a -lm -lz
SOURCE_FILES = $(wildcard main.c src/*.c)

all: clean lava
	./lava

lava:
	$(COMPILER) $(SOURCE_FILES) -o $@ $(OPT_FLAGS) $(DEBUG_FLAGS) $(LIBRARIES)

clean:
	rm -f main.out

cleandata:
	rm -f *.h5