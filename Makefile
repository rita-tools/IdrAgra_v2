# version string for the header: use an exact tag when present, otherwise the commit string
VERSION := $(shell git describe --tags --exact-match HEAD 2>/dev/null || git rev-parse --short HEAD)
CURRENTDATE := $(shell powershell -NoProfile -Command "Get-Date -Format yyyy-MM-ddTHH:mm:ss")

# Windows OS variables & settings
DEL = powershell -NoProfile -Command
EXE = .exe
WIN = 1

# Compiler settings
# -cpp: activates compiler pre processing
# -DGIT_VERSION: sets the macro GIT_VERSION in the code (actually used only in cli_main.f90)
# -g: enables debug with breakpoints

CC = gfortran
CPP = gfortran -cpp
# -g for gdb, -O0 zero optimization or -Og
### for debug ###
#GFFLAGS = -cpp -DGIT_VERSION=\"$(VERSION)\" -DCOMP_DATE=\"$(CURRENTDATE)\" -DWIN=$(WIN) -g -Wall -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all,no-array-temps -ffpe-trap=zero,overflow,underflow -finit-real=nan -c
# nnnooo  GFFLAGS = -cpp -DGIT_VERSION=\"$(VERSION)\" -DCOMP_DATE=\"$(CURRENTDATE)\" -DWIN=$(WIN) -g -Wall  -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=denorm -funsafe-math-optimizations -finit-real=nan -c
#GFFLAGS = -g -O0 -Wall -Wextra -Wshadow -pedantic -static -c
#GFFLAGS =  -cpp -DMY_VERSION=\"$(VERSION)\" -g -Wall -c
### for release ###
# -ffree-line-length-512 manage long commands in the code
GFFLAGS = -cpp -DGIT_VERSION=\"$(VERSION)\" -DCOMP_DATE=\"$(CURRENTDATE)\" -DWIN=$(WIN) -ffast-math  -O3 -ffree-line-length-0 -c
LDFLAGS =

APPNAME = idragra
EXT = .f90
SRCDIR = src
OBJDIR = obj
RELDIR = release

# List of file names, without extention separated by space
# Check the list sequence according to compile order
# interf_bilancio

FILES = mod_constants mod_utility mod_parameters mod_grid mod_common mod_evapotranspiration \
		mod_meteo \
		mod_crop_phenology mod_crop_soil_water mod_runoff mod_TDx_index mod_irrigation \
		mod_system \
		cli_watsources \
		cli_crop_parameters \
		cli_save_outputs cli_read_parameter\
		cli_simulation_manager \

#### User, don't touch the following line ####

# Builds the app
# force removing main.o in order to update program metadata
$(APPNAME): $(patsubst %, $(OBJDIR)/%.o, $(FILES))
	$(DEL) "if (Test-Path '$(OBJDIR)\cli_main.o') { Remove-Item -Force '$(OBJDIR)\cli_main.o' }"
	$(CC) -o $(OBJDIR)/cli_main.o -J$(OBJDIR) $(GFFLAGS) $(SRCDIR)/cli_main.f90
	$(CPP) -g -o $(RELDIR)/$@$(EXE) $^ $(OBJDIR)/cli_main.o $(LDFLAGS) -static

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
	$(DEL) "if (Test-Path '$(OBJDIR)\*.mod') { Remove-Item -Force '$(OBJDIR)\*.mod' }"
	$(DEL) "if (Test-Path '$(OBJDIR)\*.o') { Remove-Item -Force '$(OBJDIR)\*.o' }"
	$(DEL) "if (Test-Path '$(RELDIR)\*.exe') { Remove-Item -Force '$(RELDIR)\*.exe' }"

.PHONY: cleanmain
cleanmain:
	@echo "hello from cleanmain"
	$(DEL) "if (Test-Path '$(OBJDIR)\cli_main.o') { Remove-Item -Force '$(OBJDIR)\cli_main.o' }"
