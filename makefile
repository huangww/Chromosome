LIBS = -llapack -lm
CC = gcc
CFLAGS = -Wall -std=c99 -O3
SRCDIR = code
BUILDDIR = build
TARGET = $(BUILDDIR)/a.out

VPATH = code

.PHONY : default all run clean movie plot debug

default: $(TARGET)
all: default

SRC = $(foreach sdir, $(SRCDIR),  $(wildcard $(sdir)/*.c))
HEADERS = $(foreach sdir, $(SRCDIR),  $(wildcard $(sdir)/*.h))
OBJECTS = $(patsubst %.c, $(BUILDDIR)/%.o, $(notdir $(SRC)))

$(BUILDDIR)/%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -Wall $(LIBS) -o $@

run:
	./$(TARGET)

clean :
	-rm -rf $(BUILDDIR)/*

movie:
	vpython script/movie.py

plot:
	python script/movie.py

debug: $(TARGET)
	sudo gdb $(TARGET)
