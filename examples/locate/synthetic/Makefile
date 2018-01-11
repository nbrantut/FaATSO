CC := g++
SRCDIR := ../../../src

BUILDDIR := ../../../build
TARGET := synth_loc

LOCBUILDDIR := build
LOCBINDIR := bin

SOURCES := $(SRCDIR)/fmmiolib.cpp $(SRCDIR)/fmmlib.cpp  $(SRCDIR)/velocitylib.cpp $(SRCDIR)/numtools.cpp

OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.cpp=.o))

CFLAGS := -std=c++11 -O3 -Wall
LIB := -larmadillo
INC := -I ../../../include

all: $(LOCBINDIR)/$(TARGET)

$(LOCBINDIR)/$(TARGET): $(OBJECTS) $(LOCBUILDDIR)/$(TARGET).o
	@mkdir -p $(LOCBINDIR)
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $@ $(LIB)

$(LOCBUILDDIR)/%.o: %.cpp
	@mkdir -p $(LOCBUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(LOCBINDIR)/$(TARGET) $(LOCBUILDDIR)/$(TARGET).o"; $(RM) -r $(LOCBINDIR)/$(TARGET) $(LOCBUILDDIR)/$(TARGET).o

.PHONY: clean
