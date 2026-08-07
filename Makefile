CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic
BIN := bin

.PHONY: all clean test

all: $(BIN)/pdms_elastomer_generator $(BIN)/z_profile \
	$(BIN)/final_snapshot_analyzer $(BIN)/msd_analyzer

$(BIN):
	mkdir -p $(BIN)

$(BIN)/pdms_elastomer_generator: Generator/pdms_elastomer_generator.cpp Generator/pdms_filler_component.hpp | $(BIN)
	$(CXX) $(CXXFLAGS) $< -o $@

$(BIN)/z_profile: Analysis/z_profile.cpp | $(BIN)
	$(CXX) $(CXXFLAGS) $< -o $@

$(BIN)/final_snapshot_analyzer: Analysis/final_snapshot_analyzer.cpp | $(BIN)
	$(CXX) $(CXXFLAGS) $< -o $@

$(BIN)/msd_analyzer: Analysis/msd_analyzer.cpp | $(BIN)
	$(CXX) $(CXXFLAGS) -fopenmp $< -o $@

test: all
	bash tests/smoke_test.sh

clean:
	rm -rf $(BIN)

