CC := g++
SRCDIR := ../../../src

BUILDDIR := ../../../build
TARGET := synthetic

SOURCES := $(SRCDIR)/fmmiolib.cpp $(SRCDIR)/fmmlib.cpp  $(SRCDIR)/velocitylib.cpp $(SRCDIR)/numtools.cpp

OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.cpp=.o))

CFLAGS := -std=c++11 -O3 -Wall
LIB := -larmadillo
INC := -I ../../../include

all: $(TARGET)

$(TARGET): $(OBJECTS) synthetic.o
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)

synthetic.o: synthetic.cpp
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(TARGET) synthetic.o"; $(RM) -r $(TARGET) synthetic.o

.PHONY: clean
