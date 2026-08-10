# PDMS Elastomer

This repository generates a coarse-grained PDMS elastomer model with a simple,
component-based interface. The base model contains linear functional strands
and short functional PDMS crosslinkers. Five-bead star moderators and neutral
PDMS filler chains are optional.

## Build

```bash
make
```

Executables are written to `bin/`. To build only the generator:

```bash
make bin/pdms_elastomer_generator
```

## Components and molecule IDs

Molecule IDs are consecutive and grouped in this order:

1. strands;
2. crosslinkers;
3. moderators;
4. filler.

The default model uses 900 linear strands of 128 beads, with both chain ends
functional. Crosslinkers are 32-bead linear PDMS chains with four functional
sites each. Moderators and filler are disabled by default.

## Stoichiometry

`--stoichiometry A:B` specifies the ratio

```text
strand-end functional groups : crosslinker functional groups
```

Only the two functional ends per strand enter this calculation. Moderator
functional groups are extra and do not change the derived crosslinker count.
The generator requires an exact whole-molecule crosslinker result.

At the default `1:1` ratio, with four-functional crosslinkers and no
moderators:

```text
M_crosslinker = 2*M_strand/4 = M_strand/2
```

Therefore the default strand:crosslinker molecule ratio is `2:1`.

## Model-size limit

The complete generated model—including strands, crosslinkers, moderators, and
filler—is limited to 150,000 beads. The generator reports an error instead of
writing an oversized model.

## Text-file input

Settings can also be stored in a reusable text file using `key = value` lines:

```text
strand_length = 128
strand_count = 900
crosslinker_length = 32
functionality = 4
stoichiometry = 1:1
moderator_count = 0
spacing = 7.0
output = data.PDMS_elastomer
```

Lines beginning with `#` are comments. Keys may use underscores or hyphens.
Run the generator with:

```bash
./bin/pdms_elastomer_generator --config examples/01_default/model.conf
```

Command-line values override values from the file, which makes parameter
sweeps straightforward:

```bash
./bin/pdms_elastomer_generator \
    --config examples/01_default/model.conf \
    --strand-count 800 \
    --output data.strands_800
```

The resulting `<case>.info` file records the resolved settings, config-file
path, component counts, molecule-ID ranges, stoichiometry, and bead limit for
future analysis.

## Generate a model

```bash
./bin/pdms_elastomer_generator

./bin/pdms_elastomer_generator \
    --strand-length 96 \
    --strand-count 600 \
    --crosslinker-length 24 \
    --functionality 4 \
    --stoichiometry 1:1

./bin/pdms_elastomer_generator \
    --moderator-count 6 \
    --filler-length 32 \
    --filler-wt 10
```

Each invocation creates a case-named directory containing:

- `data.<case>`: initial LAMMPS data file;
- `in.<case>`: equilibration, crosslinking, and MSD-production input;
- `submit.<case>.sh`: one-node, 96-task Slurm script;
- `<case>.info`: version-3 JSON model metadata.

Omitting `--thickness` produces a cubic, fully periodic bulk system. Supplying
it fixes `Lz`, keeps `Lx = Ly`, uses a nonperiodic z direction, and writes
separate lower and upper wall fixes.

## Generator options

```text
--strand-length N
--strand-count M
--crosslinker-length N
--functionality F
--stoichiometry A:B
--moderator-count M
--config FILE
--filler-length N
--filler-wt X
--filler-seed N
--filler-min-separation X
--crosslink-distribution random|regular
--crosslink-seed N
--mass X
--density X
--target-density X
--bond-length X
--spacing X
--thickness X
--seed N
--output FILE
--help
```

For regular crosslinkers, reactive sites are placed at beads 1, 3, 5, and so
on. Random placement samples distinct sites reproducibly from
`--crosslink-seed`. Filler chains are packed in the central 40% of the box
while avoiding the other components and box boundaries.

## Simulation template

The generated LAMMPS input relaxes the low-density model at 800 K, compresses
and crosslinks it, cools it to 300 K, performs final equilibration, and writes
an independent 1,000,000-step NVT trajectory for MSD analysis.

The generated Slurm file is a template for the Iowa State Nova environment.
Review its modules, partition, memory, wall time, and email before use on a
different cluster.

## Analysis

The current generator writes generic version-3 metadata and the new component
order. Migration of the programs under `Analysis/` is intentionally deferred;
those sources have not been changed in this generator-only refactor.

## Repository layout

- `Generator/pdms_elastomer_generator.cpp`: generic model generator;
- `examples/01_default/`: reproducible current-default sample;
- `Generator/pdms_filler_component.hpp`: neutral PDMS filler builder;
- `Analysis/`: existing post-processing programs, currently unchanged;
- `tests/smoke_test.sh`: generator build, metadata, stoichiometry, and component-order checks.
