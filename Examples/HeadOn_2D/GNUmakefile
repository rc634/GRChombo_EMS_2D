# -*- Mode: Makefile -*-

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# Included makefiles need an absolute path to the Chombo installation
# CHOMBO_HOME := Please set the CHOMBO_HOME locally (e.g. export CHOMBO_HOME=... in bash)

ifdef TWOPUNCTURES_SOURCE
ifneq (,$(wildcard $(TWOPUNCTURES_SOURCE)/TwoPunctures.hpp))
  $(info ================== Using TwoPunctures ==================)
  USE_TWOPUNCTURES := TRUE
  # Note this code requires linking with the GNU Scientific Library
  # In addition to the line below, the include paths to the relevant header
  # files might also need to be passed to the preprocessor/compiler
  XTRALIBFLAGS = -lgsl
  cxxcppflags = -DUSE_TWOPUNCTURES
else
  $(error TWOPUNCTURES_SOURCE not set correctly)
endif
endif

DIM = 2

GRCHOMBO_SOURCE = ../../Source

ebase = Main_HeadOn2D

LibNames = AMRTimeDependent AMRTools BoxTools

src_dirs := $(GRCHOMBO_SOURCE)/utils \
            $(GRCHOMBO_SOURCE)/simd  \
            $(GRCHOMBO_SOURCE)/CCZ4  \
            $(GRCHOMBO_SOURCE)/Cartoon  \
            $(GRCHOMBO_SOURCE)/BoxUtils  \
            $(GRCHOMBO_SOURCE)/GRChomboCore  \
            $(GRCHOMBO_SOURCE)/TaggingCriteria  \
            $(GRCHOMBO_SOURCE)/AMRInterpolator  \
            $(GRCHOMBO_SOURCE)/ApparentHorizonFinder  \
            $(GRCHOMBO_SOURCE)/InitialConditions/BlackHoles \
            $(GRCHOMBO_SOURCE)/BlackHoles

ifeq ($(USE_TWOPUNCTURES),TRUE)
  src_dirs += $(TWOPUNCTURES_SOURCE)
endif

include $(CHOMBO_HOME)/mk/Make.test
