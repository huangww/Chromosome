LIBS = -llapack -lm -lblas
CC = clang++
CFLAGS = -Wall -std=c++11 -O3
SRCDIR = code
BUILDDIR = build
TARGET = $(BUILDDIR)/a.out

VPATH = code

.PHONY : default all run clean movie plot debug

default: $(TARGET)
all: default
debug: CFLAGS += -g
debug: default

SRC = $(foreach sdir, $(SRCDIR),  $(wildcard $(sdir)/*.cpp))
HEADERS = $(foreach sdir, $(SRCDIR),  $(wildcard $(sdir)/*.hpp))
OBJECTS = $(patsubst %.cpp, $(BUILDDIR)/%.o, $(notdir $(SRC)))

$(BUILDDIR)/%.o: %.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) $(CFLAGS) $(LIBS) -o $@

run:
	./$(TARGET)

clean :
	-rm -rf $(BUILDDIR)/*

movie:
	python script/movie.py
	# python script/sfd-movie.py

plot:
	python script/plotfig.py

