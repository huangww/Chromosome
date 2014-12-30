CC = gcc-4.9
CFLAGS = -Wall -std=c99 -fopenmp
TARGET = a.out
SRC = pi.c

.PHONY : default all run

default: $(TARGET)
all: default


$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $@

run:
	./$(TARGET)

