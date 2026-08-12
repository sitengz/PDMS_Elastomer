CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic
BIN := bin

.PHONY: all clean test

all: $(BIN)/pdms_elastomer_generator $(BIN)/topology_analyzer \
	$(BIN)/basic_network_analyzer $(BIN)/network_profile_analyzer

$(BIN):
	mkdir -p $(BIN)

$(BIN)/pdms_elastomer_generator: Generator/pdms_elastomer_generator.cpp Generator/pdms_filler_component.hpp | $(BIN)
	$(CXX) $(CXXFLAGS) $< -o $@

$(BIN)/topology_analyzer: Analysis/topology_analyzer.cpp Analysis/network_common.hpp | $(BIN)
	$(CXX) $(CXXFLAGS) $< -o $@

$(BIN)/basic_network_analyzer: Analysis/basic_network_analyzer.cpp Analysis/network_common.hpp | $(BIN)
	$(CXX) $(CXXFLAGS) $< -o $@

$(BIN)/network_profile_analyzer: Analysis/network_profile_analyzer.cpp Analysis/network_common.hpp | $(BIN)
	$(CXX) $(CXXFLAGS) $< -o $@

test: all
	bash tests/smoke_test.sh

clean:
	rm -rf $(BIN)
