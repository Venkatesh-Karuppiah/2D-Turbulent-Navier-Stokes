# Turbulent Navier-Stokes Solver - Makefile

# Compiler settings
FC = gfortran
FFLAGS = -O3 -fdefault-real-8 -Wall
LDFLAGS =

# Directories
SRC_DIR = src
SOLVERS_DIR = $(SRC_DIR)/solvers
BUILD_DIR = build
BIN_DIR = bin

# Create output directories
$(shell mkdir -p $(BUILD_DIR) $(BIN_DIR))

# Module dependency order (must be compiled in this order)
MOD_CONSTANTS  = $(BUILD_DIR)/constants.o
MOD_STATE      = $(BUILD_DIR)/solver_state.o
MOD_GEOMETRY   = $(BUILD_DIR)/geometry.o
MOD_VISCOUS    = $(BUILD_DIR)/viscous.o
MOD_SOURCES    = $(BUILD_DIR)/source_terms.o
MOD_BOUNDARY   = $(BUILD_DIR)/boundary.o
MOD_RESIDUAL   = $(BUILD_DIR)/residual.o

# Solver targets
TURB2D   = $(BIN_DIR)/turbnstwod
TURBAX   = $(BIN_DIR)/turbnsaxis

# All targets
all: $(TURB2D) $(TURBAX)

# 2D Solver
$(TURB2D): $(MOD_CONSTANTS) $(MOD_STATE) $(MOD_GEOMETRY) $(MOD_VISCOUS) \
           $(MOD_SOURCES) $(MOD_BOUNDARY) $(MOD_RESIDUAL) \
           $(SOLVERS_DIR)/turbnstwod.f90
	$(FC) $(FFLAGS) -I$(BUILD_DIR) -o $@ $^ $(LDFLAGS)

# Axisymmetric Solver
$(TURBAX): $(MOD_CONSTANTS) $(MOD_STATE) $(MOD_GEOMETRY) $(MOD_VISCOUS) \
           $(MOD_SOURCES) $(MOD_BOUNDARY) $(MOD_RESIDUAL) \
           $(SOLVERS_DIR)/turbnsaxis.f90
	$(FC) $(FFLAGS) -I$(BUILD_DIR) -o $@ $^ $(LDFLAGS)

# Module compilation (dependency order)
$(BUILD_DIR)/constants.o: $(SOLVERS_DIR)/constants.f90
	$(FC) $(FFLAGS) -c -o $@ $<

$(BUILD_DIR)/solver_state.o: $(SOLVERS_DIR)/solver_state.f90 $(BUILD_DIR)/constants.o
	$(FC) $(FFLAGS) -I$(BUILD_DIR) -c -o $@ $<

$(BUILD_DIR)/geometry.o: $(SOLVERS_DIR)/geometry.f90 $(BUILD_DIR)/solver_state.o
	$(FC) $(FFLAGS) -I$(BUILD_DIR) -c -o $@ $<

$(BUILD_DIR)/viscous.o: $(SOLVERS_DIR)/viscous.f90 $(BUILD_DIR)/solver_state.o
	$(FC) $(FFLAGS) -I$(BUILD_DIR) -c -o $@ $<

$(BUILD_DIR)/source_terms.o: $(SOLVERS_DIR)/source_terms.f90 $(BUILD_DIR)/solver_state.o
	$(FC) $(FFLAGS) -I$(BUILD_DIR) -c -o $@ $<

$(BUILD_DIR)/boundary.o: $(SOLVERS_DIR)/boundary.f90 $(BUILD_DIR)/solver_state.o
	$(FC) $(FFLAGS) -I$(BUILD_DIR) -c -o $@ $<

$(BUILD_DIR)/residual.o: $(SOLVERS_DIR)/residual.f90 $(BUILD_DIR)/solver_state.o
	$(FC) $(FFLAGS) -I$(BUILD_DIR) -c -o $@ $<

# Clean build artifacts
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

.PHONY: all clean
