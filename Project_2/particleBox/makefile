
CC = gcc # Compiler 
CFLAGS = -Wall -Wextra -O2 -g -Iinclude # Compiler flags: -Wall (enable all warnings), -Wextra (enable extra warnings), -O2 (optimization level 2), -g (include debug information), -Iinclude (add include directory)
LDFLAGS = -lm # Linker flags: -lm (link math library)

SRC = src/config.c src/particle.c src/sectors.c src/forces.c src/placement.c src/spatialHashing.c src/quicksort.c src/simulation.c  src/main.c
OBJ = $(SRC:.c=.o)

BIN = runld

all: $(BIN)

$(BIN): $(OBJ)
	$(CC) $(OBJ) -o $@ $(LDFLAGS)

test: $(OBJ)
	$(CC) $(OBJ) -o testLangevin $(LDFLAGS)


clean:
	rm -f $(OBJ) $(BIN)

