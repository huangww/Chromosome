CC = clang++
SRCDIR = code
BUILDDIR = build
TARGET = $(BUILDDIR)/a.out
CXXFLAGS += -Wall -std=c++11 -O3
CXXFLAGS += -I/usr/local/include
LDFLAGS += -L/usr/local/lib
LDLIBS = -llapack -lm -lgsl -lgslcblas 

VPATH = code

.PHONY : default all run clean movie plot debug

default: $(TARGET)
all: default
debug: CXXFLAGS += -g
debug: default

SRC = $(foreach sdir, $(SRCDIR),  $(wildcard $(sdir)/*.cpp))
HEADERS = $(foreach sdir, $(SRCDIR),  $(wildcard $(sdir)/*.hpp))
OBJECTS = $(patsubst %.cpp, $(BUILDDIR)/%.o, $(notdir $(SRC)))

$(BUILDDIR)/%.o: %.cpp $(HEADERS)
	$(CC) $(CXXFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@ $(LDLIBS)

run:
	./$(TARGET) input.in

clean :
	-rm -rf $(BUILDDIR)/*

movie:
	python script/movie.py
	# python script/sfd-movie.py

plot:
	python script/plotfig.py
