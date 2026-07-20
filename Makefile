# version string for the header: use an exact tag when present, otherwise the commit string
VERSION := $(shell git describe --tags --exact-match HEAD 2>/dev/null || git rev-parse --short HEAD)
CURRENTDATE := $(shell date --iso=seconds)
LOG_FILE := $(shell git log --format=reference > ./release/git_idragra.log)

# Windows OS variables & settings
DEL = rm
CP = cp
EXE = .exe
WIN = 1

# Compiler settings
# -cpp: activates compiler pre processing
# -DGIT_VERSION: sets the macro GIT_VERSION in the code (actually used only in cli_main.f90)
# -g: enables debug with breakpoints 

MINGW64_BINDIR ?= /mingw64/bin
export PATH := $(MINGW64_BINDIR):$(PATH)

CC = gfortran
CPP = gfortran -cpp
IS_RELEASE := false

# -g for gdb, -O0 zero optimization or -Og
### for debug ###
TYPE := debug
GFFLAGS = -cpp -DGIT_VERSION=\"$(VERSION)\" -DCOMP_DATE=\"$(CURRENTDATE)\" -DWIN=$(WIN) -g -Wall  -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan -c 
# nnnooo  GFFLAGS = -cpp -DGIT_VERSION=\"$(VERSION)\" -DCOMP_DATE=\"$(CURRENTDATE)\" -DWIN=$(WIN) -g -Wall  -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=denorm -funsafe-math-optimizations -finit-real=nan -c 
#GFFLAGS = -g -O0 -Wall -Wextra -Wshadow -pedantic -static -c
#GFFLAGS =  -cpp -DMY_VERSION=\"$(VERSION)\" -g -Wall -c

ifeq ($(IS_RELEASE),true)
   ### for release ###
   # -ffree-line-length-512 manage long commands in the code
   TYPE := release
   GFFLAGS = -cpp -DGIT_VERSION=\"$(VERSION)\" -DCOMP_DATE=\"$(CURRENTDATE)\" -DWIN=$(WIN) -ffast-math  -O3 -ffree-line-length-0 -c
endif

LDFLAGS = 

APPNAME = idragra_$(VERSION)_$(TYPE)
APPALIAS = idragra_latest
EXT = .f90
SRCDIR = src
OBJDIR = obj
RELDIR = release

# List of file names, without extention separated by space
# Check the list sequence according to compile order
# interf_bilancio 

# FILES = mod_constants utility mod_parameters mod_common mod_grid mod_meteo mod_xls mod_crop_parameters \
# 		mod_watsources crop_cycle_distr interf curve_number evaporation_iter \
# 		percolation cap_rise transpiration mod_runoff serbatoio1 serbatoio2 TDx_index mod_irrigation \
# 		mod_CSW_balance \
# 		print_debug \
# 		main

FILES = mod_constants mod_utility mod_parameters mod_grid mod_common mod_evapotranspiration \
		mod_meteo \
		mod_crop_phenology mod_crop_soil_water mod_runoff mod_TDx_index mod_irrigation \
		mod_system \
		cli_watsources \
		cli_crop_parameters \
		cli_save_outputs cli_read_parameter\
		cli_simulation_manager \
		print_debug
		

#### User, don't touch the following line ####

# Builds the app
# force removing main.o in order to update program metadata
$(APPNAME): $(patsubst %, $(OBJDIR)/%.o, $(FILES))
	$(DEL) -f /$(OBJDIR)/cli_main.o
	$(CC) -o $(OBJDIR)/cli_main.o -J$(OBJDIR) $(GFFLAGS) $(SRCDIR)/cli_main.f90
	$(CPP) -g -o $(RELDIR)/$@$(EXE) $^ $(OBJDIR)/cli_main.o $(LDFLAGS) -static
	$(CP) $(RELDIR)/$@$(EXE) $(RELDIR)/$(APPALIAS)$(EXE)

# The following sets a function that creates makefile rules
# Additionally:
# -o: set path to output
# -J: set path for *.mod
define make-o-rule
$(OBJDIR)/$1.o: $(SRCDIR)/$1.f90
	$(CC) -o $(OBJDIR)/$1.o -J$(OBJDIR) $(GFFLAGS) $(SRCDIR)/$1.f90
all: $1.o
endef

$(foreach element,$(FILES),$(eval $(call make-o-rule,$(element))))

all:
	$(APPNAME)

##### Set the obj folder empty ####
# Cleans complete project
# call as: make cleanall
.PHONY: cleanall
cleanall:
	$(DEL) $(wildcard ./$(OBJDIR)/*.mod)
	$(DEL) $(wildcard ./$(OBJDIR)/*.o)
	$(DEL) $(wildcard ./$(RELDIR)/*.exe)

.PHONY: cleanmain
cleanmain:
	@echo "hello from cleanmain"
	$(DEL) -f /$(OBJDIR)/main.o
