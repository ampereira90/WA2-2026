# ACA - Work Assignment 2 - 3D Jacobi
#
#   make [COMPILER=<toolchain>] [OPT="extra flags"] [VARIANTS="poisson7 box27 varcoef"]
#
# Toolchains (load the matching module first, see README.md):
#   gnu            g++            x86 and Arm (default)
#   fujitsu        FCC  trad mode native compiler, A64FX compute nodes only
#   fujitsu-clang  FCC -Nclang    native compiler, A64FX compute nodes only
#   fujitsu-cross  FCCpx          cross compiler on the x86 login nodes
#   intel          icpx           x86 only
#   nvhpc          nvc++          x86 (GPU nodes)
#   clang          clang++        any
#
# Binaries go to bin/<arch>-<toolchain>/ because the same directory is seen
# from the x86 login nodes and from the Arm compute nodes.
#
# The REFERENCE build uses the flags below with OPT empty. Record the exact
# flags you use for every measurement (they are printed by `make info`).

COMPILER ?= gnu
VARIANTS ?= poisson7 box27 varcoef
OPT      ?=

ARCH := $(shell uname -m)
ifeq ($(COMPILER),fujitsu-cross)
  ARCH := aarch64
endif

# ---------------------------------------------------------------- toolchains
ifeq ($(COMPILER),gnu)
  CXX       = g++
  BASEFLAGS = -O2 -g
  ifeq ($(ARCH),x86_64)
    ARCHFLAGS = -march=native
  else
    ARCHFLAGS = -mcpu=native
  endif
endif
ifeq ($(COMPILER),fujitsu)
  CXX       = FCC
  BASEFLAGS = -O2 -g
  ARCHFLAGS = -KA64FX -KSVE
endif
ifeq ($(COMPILER),fujitsu-clang)
  CXX       = FCC
  BASEFLAGS = -Nclang -O2 -g
  ARCHFLAGS = -mcpu=a64fx
endif
ifeq ($(COMPILER),fujitsu-cross)
  CXX       = FCCpx
  BASEFLAGS = -O2 -g
  ARCHFLAGS = -KA64FX -KSVE
endif
ifeq ($(COMPILER),intel)
  CXX       = icpx
  BASEFLAGS = -O2 -g -fp-model=precise
  ARCHFLAGS = -march=core-avx2
endif
ifeq ($(COMPILER),nvhpc)
  CXX       = nvc++
  BASEFLAGS = -O2 -g
  ARCHFLAGS = -tp=zen2
endif
ifeq ($(COMPILER),clang)
  CXX       = clang++
  BASEFLAGS = -O2 -g
  ARCHFLAGS = -march=native
endif
ifndef CXX
  $(error unknown COMPILER=$(COMPILER))
endif

STD      ?= -std=c++11
CXXFLAGS  = $(BASEFLAGS) $(ARCHFLAGS) $(STD) $(OPT)
LDLIBS    = -lm

BINDIR   := bin/$(ARCH)-$(COMPILER)
COMMON    = common/common.cpp
SRC_poisson7 = poisson7/main.cpp poisson7/grid.cpp poisson7/jacobi.cpp
SRC_box27    = box27/main.cpp    box27/grid.cpp    box27/jacobi.cpp box27/stencil.cpp
SRC_varcoef  = varcoef/main.cpp  varcoef/grid.cpp  varcoef/jacobi.cpp

.PHONY: all clean info
all: $(addprefix $(BINDIR)/,$(VARIANTS)) $(BINDIR)/compare

$(BINDIR):
	mkdir -p $@

OBJDIR   := $(BINDIR)/obj

# each translation unit is compiled separately (no LTO in the reference)
$(OBJDIR)/%.o: %.cpp | $(BINDIR)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

.SECONDEXPANSION:
$(BINDIR)/poisson7 $(BINDIR)/box27 $(BINDIR)/varcoef: $$(patsubst %.cpp,$(OBJDIR)/%.o,$$(SRC_$$(@F)) $(COMMON))
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

$(BINDIR)/compare: $(OBJDIR)/tools/compare.o $(OBJDIR)/common/common.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

info:
	@echo "ARCH     = $(ARCH)"
	@echo "COMPILER = $(COMPILER)"
	@echo "CXX      = $(CXX)"
	@echo "CXXFLAGS = $(CXXFLAGS)"
	@$(CXX) --version 2>/dev/null | head -1 || true

clean:
	rm -rf bin
