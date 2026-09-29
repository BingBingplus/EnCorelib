##############################################################################
#  Makefile  --  EndCover for IVP  (EnCorelib)
#
#  Quick start
#  -----------
#      make                       compile ./EnCorelib
#      make run-eg1Volterra       run Eg1 of Table 1
#      make table2                Table 2 of the paper
#      make figure1  figure4      the data behind Figures 1 and 4
#      make check                 self-test
#      make compare               diff the last table2 run against the paper
#      make verify                sample-point check of the covers
#
#  Selecting CAPD
#  --------------
#  The Makefile looks for a CAPD source tree that contains build/libcapd.a,
#  or falls back to capd-config / pkg-config.  Override explicitly with
#
#      make CAPD_MASTER=/path/to/CAPD           (source tree + build/)
#      make CAPD_CONFIG=/path/to/capd-config    (installed CAPD)
#
#  Running an example
#  ------------------
#      make run EG=examples/eg5Lorenz.mk
#      make run EG=examples/eg5Lorenz.mk eps=0.1 T=4 method=0 output_mode=4
#
#  Every field of a .mk file can be overridden on the command line; see
#  examples/_template.txt for the list.
##############################################################################

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wno-sign-compare -Wno-overloaded-virtual

TARGET  := EnCorelib
HEADERS := types.h symparse.h stepAB.h extend.h refine.h endcover.h boundary.h

# ────────────────────────────────────────────────────────────────────────────
#  1. Locate CAPD
# ────────────────────────────────────────────────────────────────────────────
#  THIS IS THE ONLY THING YOU HAVE TO CONFIGURE.  Pick whichever is easiest:
#
#    (a) copy config.mk.example to config.mk and edit the one path in it
#          cp config.mk.example config.mk        # then edit config.mk
#
#    (b) pass it on the command line, once per build
#          make CAPD_MASTER=/absolute/path/to/CAPD
#
#    (c) do nothing, if capd-config is already on your PATH
#
#  config.mk is yours; it is git-ignored and never shipped.
-include config.mk

#  Directories searched when neither config.mk nor the command line says where
#  CAPD is.  A directory qualifies when it contains build/libcapd.a.
_CAPD_CANDIDATES := \
  $(CAPD_MASTER) \
  $(CAPD_DIR) \
  $(HOME)/CAPD \
  $(HOME)/capd \
  /opt/CAPD \
  /usr/local/share/capd

CAPD_MASTER := $(shell \
  for d in $(_CAPD_CANDIDATES); do \
    if [ -f "$$d/build/libcapd.a" ]; then echo "$$d"; break; fi; \
  done)

CAPD_CONFIG ?= $(shell command -v capd-config 2>/dev/null)

ifneq ($(CAPD_MASTER),)
  CAPD_BUILD  := $(CAPD_MASTER)/build
  CAPD_FILIB  := $(CAPD_BUILD)/capdExt/filibsrc
  CAPD_CFLAGS := \
    -I$(CAPD_MASTER)/capdDynSys/include \
    -I$(CAPD_MASTER)/capdAlg/include \
    -I$(CAPD_MASTER)/capdAux/include \
    -I$(CAPD_MASTER)/capdExt/include \
    -I$(CAPD_MASTER)/capdExt/filibsrc \
    -D__USE_FILIB__ -DFILIB_EXTENDED -frounding-math
  CAPD_LIBS   := -L$(CAPD_BUILD) -L$(CAPD_FILIB) -lcapd -lfilib
  $(info CAPD source tree: $(CAPD_MASTER))
else ifneq ($(CAPD_CONFIG),)
  CAPD_CFLAGS := $(shell $(CAPD_CONFIG) --cflags) -frounding-math
  CAPD_LIBS   := $(shell $(CAPD_CONFIG) --libs)
  $(info CAPD via $(CAPD_CONFIG))
else
  $(warning )
  $(warning CAPD not found.  Build CAPD first (see README.md) and then use)
  $(warning     make CAPD_MASTER=/path/to/CAPD)
  $(warning or  make CAPD_CONFIG=/path/to/capd-config)
  $(warning )
  CAPD_CFLAGS :=
  CAPD_LIBS   :=
endif

# ────────────────────────────────────────────────────────────────────────────
#  2. Optional SymEngine (algebraic expansion of the right-hand side)
# ────────────────────────────────────────────────────────────────────────────
SYMENGINE_AUTO := $(shell \
  for d in "$(SYMENGINE_DIR)" /usr/local /usr; do \
    [ -f "$$d/lib/libsymengine.a" ] && echo "$$d" && break; \
  done)

ifneq ($(SYMENGINE_AUTO),)
  SYMENGINE_CFLAGS := -DHAVE_SYMENGINE -I$(SYMENGINE_AUTO)/include
  SYMENGINE_LIBS   := -L$(SYMENGINE_AUTO)/lib -l:libsymengine.a -lgmp -lmpfr -lmpc
  $(info SymEngine: $(SYMENGINE_AUTO)  (optional))
else
  SYMENGINE_CFLAGS :=
  SYMENGINE_LIBS   :=
endif

# ────────────────────────────────────────────────────────────────────────────
#  3. Build
# ────────────────────────────────────────────────────────────────────────────
all: $(TARGET)

$(TARGET): EnCorelib.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(CAPD_CFLAGS) $(SYMENGINE_CFLAGS) \
	    EnCorelib.cpp -o $(TARGET) \
	    $(CAPD_LIBS) $(SYMENGINE_LIBS)
	@echo "Built $(TARGET)"

# Encoretraj: validated trajectory printer used for Figure 1 (top row).
Encoretraj: tools/trajectory.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(CAPD_CFLAGS) $(SYMENGINE_CFLAGS) \
	    tools/trajectory.cpp -o Encoretraj \
	    $(CAPD_LIBS) $(SYMENGINE_LIBS)
	@echo "Built Encoretraj"

# capdref: the reference integrator behind repro/competitors/.
capdref: repro/competitors/capd_compare.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(CAPD_CFLAGS) \
	    repro/competitors/capd_compare.cpp -o capdref \
	    $(CAPD_LIBS)
	@echo "Built capdref"

tools: Encoretraj capdref

clean:
	rm -f $(TARGET) $(TARGET).exe Encoretraj Encoretraj.exe capdref capdref.exe
	rm -f E0.txt E1.txt E_0.txt E_1.txt convex_hull.txt convex_hull_E0.txt out.txt

distclean: clean
	rm -rf build_cmake repro/out

.PHONY: all clean distclean tools

# ────────────────────────────────────────────────────────────────────────────
#  4. Example runner
# ────────────────────────────────────────────────────────────────────────────
#  Defaults (a .mk file included through EG= overrides them; a variable given
#  on the command line overrides both).
n           := 2
var         := x,y
ff          := 2*x-2*x*y,-y+x*y
cen         := 1.0,3.0
wid         := 0.1,0.1
lohi        :=
eps         := 1.0
order       := 20
T           := 1.0
output_mode := 2
mode        := 0
method      := 0
stepB       := 0
tubedegree  := 0
debug       := 0

ifneq ($(EG),)
  include $(EG)
endif

# cen/wid -> "lo1 hi1 lo2 hi2 ..."  (skipped when lohi is given explicitly)
BOX = $(if $(lohi),$(lohi),$(shell awk -v c='$(cen)' -v w='$(wid)' \
        'BEGIN{nc=split(c,C,","); split(w,W,","); \
               for(i=1;i<=nc;i++) printf "%.17g %.17g ",C[i]-W[i],C[i]+W[i]}'))

# `var` and `ff` stay quoted all the way into tools/run.sh, which splits them
# on commas.
RUNARGS = $(output_mode) $(mode) $(method) $(stepB) $(tubedegree) $(n) \
          '$(var)' '$(ff)' $(eps) $(order) $(T) $(debug) $(BOX)

run: $(TARGET)
	@bash tools/run.sh ./$(TARGET) $(RUNARGS)

show-args:
	@SHOW_ONLY=1 bash tools/run.sh ./$(TARGET) $(RUNARGS)

.PHONY: run show-args

# ── Per-example shortcuts ──────────────────────────────────────────────────
run-eg1Volterra:        ; @$(MAKE) --no-print-directory run EG=examples/eg1Volterra.mk
run-eg2VanDerPol:       ; @$(MAKE) --no-print-directory run EG=examples/eg2VanDerPol.mk
run-eg3Quadratic:       ; @$(MAKE) --no-print-directory run EG=examples/eg3Quadratic.mk
run-eg4FitzHughNagumo:  ; @$(MAKE) --no-print-directory run EG=examples/eg4FitzHughNagumo.mk
run-eg5Lorenz:          ; @$(MAKE) --no-print-directory run EG=examples/eg5Lorenz.mk
run-eg6Rossler:         ; @$(MAKE) --no-print-directory run EG=examples/eg6Rossler.mk

# Generic: make run-endcover-eg5Lorenz / make run-boundary-eg5Lorenz
run-endcover-%: $(TARGET) ; @$(MAKE) --no-print-directory run EG=examples/$*.mk mode=0
run-boundary-%: $(TARGET) ; @$(MAKE) --no-print-directory run EG=examples/$*.mk mode=1

.PHONY: run-eg1Volterra run-eg2VanDerPol run-eg3Quadratic \
        run-eg4FitzHughNagumo run-eg5Lorenz run-eg6Rossler

# ────────────────────────────────────────────────────────────────────────────
#  5. The tables and figures of the paper
# ────────────────────────────────────────────────────────────────────────────
table2: $(TARGET)
	@bash repro/table2.sh

table2-%: $(TARGET)
	@bash repro/table2.sh $*

# Cell-by-cell comparison of the last table2 run against the printed table.
compare:
	@$(if $(PYTHON),$(PYTHON),python3) repro/compare.py

figure1: $(TARGET) Encoretraj
	@bash repro/figure1.sh

figure4: $(TARGET)
	@bash repro/figure4.sh

check: $(TARGET)
	@bash repro/check.sh

# Integrate a grid of sample points of B0 with ./Encoretraj and confirm each
# image lies inside some box of the cover.
verify: $(TARGET) Encoretraj
	@bash repro/verify_cover.sh

.PHONY: table2 table2-% compare figure1 figure4 check verify

# ────────────────────────────────────────────────────────────────────────────
#  6. CMake alternative
# ────────────────────────────────────────────────────────────────────────────
cmake-build:
	@mkdir -p build_cmake && cd build_cmake && \
	cmake .. -DCMAKE_BUILD_TYPE=Release && \
	cmake --build . --parallel

.PHONY: cmake-build

# ────────────────────────────────────────────────────────────────────────────
#  7. extras/
# ────────────────────────────────────────────────────────────────────────────
extras/moore_direct_test: extras/moore_direct_test.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(CAPD_CFLAGS) extras/moore_direct_test.cpp \
	    -o $@ $(CAPD_LIBS)

extras/moore_endcover_direct: extras/moore_endcover_direct.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(CAPD_CFLAGS) extras/moore_endcover_direct.cpp \
	    -o $@ $(CAPD_LIBS)

extras-all: extras/moore_direct_test extras/moore_endcover_direct

.PHONY: extras-all
