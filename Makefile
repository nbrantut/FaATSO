CC := g++
SRCDIR := src

BUILDDIR := build
TARGET_FMM := bin/fmm
TARGET_TOMO := bin/tomo

LIBSRC_FMM := $(SRCDIR)/fmmiolib.cpp $(SRCDIR)/fmmlib.cpp $(SRCDIR)/raylib.cpp $(SRCDIR)/velocitylib.cpp $(SRCDIR)/numtools.cpp

SOURCES_FMM := $(LIBSRC_FMM) $(SRCDIR)/main_fmm.cpp
OBJECTS_FMM := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES_FMM:.cpp=.o))

SOURCES_TOMO := $(LIBSRC_FMM) $(SRCDIR)/tomolib.cpp $(SRCDIR)/main_tomo.cpp
OBJECTS_TOMO := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES_TOMO:.cpp=.o))

CFLAGS := -std=c++11 -O3 -Wall
LIB := -larmadillo
INC := -I include

all: $(TARGET_FMM) $(TARGET_TOMO)

fmm: $(TARGET_FMM)

tomo: $(TARGET_TOMO)

$(TARGET_FMM): $(OBJECTS_FMM)
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET_FMM)"; $(CC) $^ -o $(TARGET_FMM) 

$(TARGET_TOMO): $(OBJECTS_TOMO)
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET_TOMO) $(LIB)"; $(CC) $^ -o $(TARGET_TOMO) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET_FMM) $(TARGET_TOMO)"; $(RM) -r $(BUILDDIR) $(TARGET_FMM) $(TARGET_TOMO)

.PHONY: clean
