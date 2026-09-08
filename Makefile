.PHONY: all build test clean python-reference

FC = gfortran
FFLAGS ?= -std=f95 -Wall -Wextra -O2
BUILD_DIR ?= build

all: test

$(BUILD_DIR)/brillouin_demo: fortran/brillouin90.f90 programs/brillouin_demo.f90
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -o $@ $^

$(BUILD_DIR)/test_brillouin: fortran/brillouin90.f90 tests_fortran/test_brillouin.f90
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) -fcheck=all -J$(BUILD_DIR) -o $@ $^

build: $(BUILD_DIR)/brillouin_demo

test: $(BUILD_DIR)/test_brillouin
	$(BUILD_DIR)/test_brillouin

python-reference:
	@echo "The former Python implementation is documented in legacy/python_reference."
	@echo "The maintained numerical implementation is the Fortran 90 code."

clean:
	rm -rf $(BUILD_DIR) *.mod
