LIBS = -llapack -lm -lblas
CC = clang++
CFLAGS = -Wall -std=c++11
SRCDIR = code
BUILDDIR = build
TARGET = $(BUILDDIR)/a.out

VPATH = code

.PHONY : default debug eigen run clean movie plot
default: CFLAGS += -O3 -DPLAIN
default: $(TARGET) 

debug: CFLAGS += -g -DDEBUG 
debug: $(TARGET)

eigen: CFLAGS += -O3 -DEIGEN 
eigen: $(TARGET)



SRC = $(foreach sdir, $(SRCDIR),  $(wildcard $(sdir)/*.cpp))
HEADERS = $(foreach sdir, $(SRCDIR),  $(wildcard $(sdir)/*.hpp))
OBJECTS = $(patsubst %.cpp, $(BUILDDIR)/%.o, $(notdir $(SRC)))

$(BUILDDIR)/%.o: %.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) $(CFLAGS) $(LIBS) -o $@

run:
	./$(TARGET) input.br.in

clean :
	-rm -rf $(BUILDDIR)/*

movie:
	python script/movie.py
	# python script/sfd-movie.py

plot:
	python script/plotfig.py

