LIBS = -llapack -lm -lblas
CC = gcc-4.9
CFLAGS = -Wall -std=c99 -g
SRCDIR = code
BUILDDIR = build
TARGET = $(BUILDDIR)/a.out

VPATH = code
# VPATH = used/plainC

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
	$(CC) $(OBJECTS) $(CFLAGS) $(LIBS) -o $@

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
