# Top-level Makefile for MARLEY
# Dispatches to CMake (preferred) or the hand-written GNU Make recipe (make/build.mk).
#
# Usage:
#   make                   — auto-detect cmake; use it if found
#   make IGNORE_CMAKE=1    — force use of the hand-written make/build.mk recipe
#   make IGNORE_ROOT=1     — pass through to either build system
#   make IGNORE_HEPMC3=1   — pass through to either build system
#   make IGNORE_GSL=1      — use built-in GSL subset (pass through to either build system)
#   make debug             — debug build
#   make clean             — remove build/ entirely
#   make install           — install (GNU Make path only for now)
#   make help              — show this message

CMAKE := $(shell command -v cmake 2>/dev/null)

USE_CMAKE :=
ifndef IGNORE_CMAKE
ifdef CMAKE
USE_CMAKE := 1
endif
endif

CMAKE_FLAGS :=
ifdef IGNORE_ROOT
CMAKE_FLAGS += -DIGNORE_ROOT=ON
endif
ifdef IGNORE_HEPMC3
CMAKE_FLAGS += -DIGNORE_HEPMC3=ON
endif
ifdef IGNORE_GSL
CMAKE_FLAGS += -DIGNORE_GSL=ON
endif

CMAKE_BUILD_FLAGS :=
MAKE_JFLAG := $(firstword $(filter -j%,$(MAKEFLAGS)))
ifneq ($(MAKE_JFLAG),)
ifneq ($(MAKE_JFLAG),-j)
CMAKE_BUILD_FLAGS += --parallel $(patsubst -j%,%,$(MAKE_JFLAG))
else
CMAKE_BUILD_FLAGS += --parallel
endif
endif

PASSTHROUGH_VARS :=
ifdef IGNORE_ROOT
PASSTHROUGH_VARS += IGNORE_ROOT=$(IGNORE_ROOT)
endif
ifdef IGNORE_HEPMC3
PASSTHROUGH_VARS += IGNORE_HEPMC3=$(IGNORE_HEPMC3)
endif
ifdef IGNORE_GSL
PASSTHROUGH_VARS += IGNORE_GSL=$(IGNORE_GSL)
endif
ifdef CXX
PASSTHROUGH_VARS += CXX='$(CXX)'
endif
ifdef CXXFLAGS
PASSTHROUGH_VARS += CXXFLAGS='$(CXXFLAGS)'
endif
ifdef prefix
PASSTHROUGH_VARS += prefix='$(prefix)'
endif
ifdef DESTDIR
PASSTHROUGH_VARS += DESTDIR='$(DESTDIR)'
endif

.PHONY: all clean debug reconfigure test docs doxygen install uninstall help

all:
ifdef USE_CMAKE
	@test -d build/CMakeFiles || cmake -S . -B build $(CMAKE_FLAGS)
	cmake --build build $(CMAKE_BUILD_FLAGS)
else
	@mkdir -p build
	$(MAKE) -C build -f $(CURDIR)/make/build.mk TOP_DIR=$(CURDIR) $(PASSTHROUGH_VARS)
endif

debug:
ifdef USE_CMAKE
	cmake -S . -B build $(CMAKE_FLAGS) -DCMAKE_BUILD_TYPE=Debug
	cmake --build build $(CMAKE_BUILD_FLAGS)
else
	@mkdir -p build
	$(MAKE) -C build -f $(CURDIR)/make/build.mk TOP_DIR=$(CURDIR) debug $(PASSTHROUGH_VARS)
endif

test:
ifdef USE_CMAKE
	@test -d build/CMakeFiles || cmake -S . -B build $(CMAKE_FLAGS)
	cmake --build build --target martest $(CMAKE_BUILD_FLAGS)
	cd $(CURDIR) && MARLEY=$(CURDIR) ctest --test-dir build --output-on-failure
else
	@mkdir -p build
	$(MAKE) -C build -f $(CURDIR)/make/build.mk TOP_DIR=$(CURDIR) test $(PASSTHROUGH_VARS)
endif

docs:
ifdef USE_CMAKE
	@test -d build/CMakeFiles || cmake -S . -B build $(CMAKE_FLAGS)
	cmake --build build --target docs $(CMAKE_BUILD_FLAGS)
else
	@mkdir -p build
	$(MAKE) -C build -f $(CURDIR)/make/build.mk TOP_DIR=$(CURDIR) docs $(PASSTHROUGH_VARS)
endif

doxygen:
ifdef USE_CMAKE
	@test -d build/CMakeFiles || cmake -S . -B build $(CMAKE_FLAGS)
	cmake --build build --target doxygen $(CMAKE_BUILD_FLAGS)
else
	@mkdir -p build
	$(MAKE) -C build -f $(CURDIR)/make/build.mk TOP_DIR=$(CURDIR) doxygen $(PASSTHROUGH_VARS)
endif

install:
ifdef USE_CMAKE
	@test -d build/CMakeFiles || cmake -S . -B build $(CMAKE_FLAGS)
	cmake --install build
else
	@mkdir -p build
	$(MAKE) -C build -f $(CURDIR)/make/build.mk TOP_DIR=$(CURDIR) install $(PASSTHROUGH_VARS)
endif

uninstall:
ifdef USE_CMAKE
	@test -d build/CMakeFiles || (echo "ERROR: No CMake build found in build/. Run 'cmake --install' first." && exit 1)
	cmake --build build --target uninstall
else
	@mkdir -p build
	$(MAKE) -C build -f $(CURDIR)/make/build.mk TOP_DIR=$(CURDIR) uninstall $(PASSTHROUGH_VARS)
endif

reconfigure:
ifdef USE_CMAKE
	$(RM) -rf build/CMakeFiles
	cmake -S . -B build $(CMAKE_FLAGS)
else
	@echo "reconfigure is only supported for the CMake build path."
endif

clean:
	$(RM) -rf build

help:
	@echo "Usage:"
	@echo "  make                   Auto-detect cmake and use it if found"
	@echo "  make IGNORE_CMAKE=1    Force GNU Make path (make/build.mk)"
	@echo "  make IGNORE_ROOT=1     Ignore ROOT"
	@echo "  make IGNORE_HEPMC3=1   Ignore external HepMC3"
	@echo "  make IGNORE_GSL=1      Use built-in GSL subset"
	@echo "  make [debug|test|docs|doxygen|install|uninstall|clean|reconfigure]"
