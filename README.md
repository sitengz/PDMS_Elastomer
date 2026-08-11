# PDMS Elastomer

This repository generates a coarse-grained PDMS elastomer model with a simple,
component-based interface. Functional strands can be linear chains, rings,
star polymers, or grafted-backbone polymers spanning comb-like and
bottlebrush architectures. The base model also contains short functional PDMS
crosslinkers. Five-bead star moderators and neutral PDMS filler chains are
optional.

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
functional. Ring strands can instead carry a configurable number of regular or
randomly located functional sites. Star strands have functional outer arm ends
and support 3, 4, 6, or 8 arms. Grafted strands place side chains along a
neutral-ended backbone and select a configurable fraction of side-chain ends
as functional. Crosslinkers are 32-bead linear PDMS chains with four functional
sites each. Moderators and filler are disabled by default.

The 2.801 Å bond length and 7.5 Å initial intermolecular spacing are fixed by
the PDMS model rather than treated as generator inputs. The spacing represents
the approximate minimum-energy separation at 800 K.

## Stoichiometry

`--stoichiometry A:B` specifies the ratio

```text
strand functional groups : crosslinker functional groups
```

The strand contribution is `strand_functionality * strand_count`. Moderator
functional groups are extra and do not change the derived crosslinker count.
The generator requires an exact whole-molecule crosslinker result.

At the default `1:1` ratio, with four-functional crosslinkers and no
moderators:

```text
M_crosslinker = F_strand*M_strand/4 = 2*M_strand/4 = M_strand/2
```

Therefore the default strand:crosslinker molecule ratio is `2:1`.

## Model-size guidance

150,000 total beads is the recommended working size for balancing simulation
cost against finite-size effects. It is not a hard generator limit. Larger
models are generated normally, with a warning so the increased memory and
runtime requirements are intentional.

## Text-file input

Settings can also be stored in a reusable text file using `key = value` lines:

```text
strand_length = 128
strand_count = 900
strand_topology = linear
strand_functionality = 2
crosslinker_length = 32
functionality = 4
stoichiometry = 1:1
moderator_count = 0
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
path, component counts, molecule-ID ranges, stoichiometry, and whether the
model exceeds the recommended bead count for future analysis.

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

./bin/pdms_elastomer_generator \
    --strand-topology ring \
    --strand-functionality 4 \
    --strand-reactive-distribution random \
    --strand-reactive-seed 20260810

./bin/pdms_elastomer_generator \
    --strand-topology star \
    --strand-length 32 \
    --strand-arm-count 6 \
    --strand-count 600

./bin/pdms_elastomer_generator \
    --strand-topology grafted \
    --backbone-length 64 \
    --side-chain-length 8 \
    --graft-spacing 12 \
    --graft-functional-fraction 40 \
    --strand-count 1200
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
--strand-topology linear|ring|star|grafted
--strand-functionality F
--strand-arm-count 3|4|6|8
--backbone-length N
--side-chain-length N
--graft-spacing N
--graft-functional-fraction X
--strand-reactive-distribution regular|random
--strand-reactive-seed N
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
--thickness X
--seed N
--output FILE
--help
```

Linear strands retain two functional ends. Ring strands use cyclic bonds,
angles, and dihedrals, and their functionality is set with
`--strand-functionality`. The `regular` distribution spaces ring reactive sites
uniformly, which requires the ring length to be divisible by its functionality.
The `random` distribution samples distinct ring sites reproducibly from
`--strand-reactive-seed`.

For star strands, `--strand-length` is the number of beads in each arm and
`--strand-arm-count` selects one of four architectures. Every outer arm end is
functional, so star functionality equals its arm count. If the length is not
specified, star arms default to 32 beads. Center beads are neutral and
connected as follows:

| Arms | Center beads | Arm distribution across centers |
|---:|---:|:---|
| 3 | 1 | 3 |
| 4 | 1 | 4 |
| 6 | 2 | 3 + 3 |
| 8 | 3 | 3 + 2 + 3 |

The multi-bead cores are linear. For the 6- and 8-arm architectures, every
center bead has total bond coordination four. A complete star contains
`center_beads + arm_count * strand_length` beads.

For grafted strands, `--backbone-length` and `--side-chain-length` replace
`--strand-length`. The unified `grafted` topology covers both comb-like and
bottlebrush structures. `--graft-spacing S` is the number of ungrafted
backbone beads between adjacent grafts, so the graft interval is `S + 1` and
the number of side chains is

```text
ceil(backbone_length / (graft_spacing + 1))
```

Spacing 0 adds one side chain to every backbone bead and is labeled a dense
bottlebrush. Positive spacing produces a comb-like structure; it is a
continuous architectural control rather than a hard comb/bottlebrush boundary.
The total strand size is
`backbone_length + side_chain_count * side_chain_length`.

`--graft-functional-fraction X` requests the percentage of side-chain ends
that are functional. The generator rounds this to a whole number of side
chains, selects those ends evenly along the graft sequence, and records both
the requested and realized fractions in the `.info` file. At least one
side-chain end is functional. Backbone ends are always neutral and are not
included in stoichiometry.

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

The topology-first analyzers read `data.<case>.npt_eq` and the matching
version-3 `<case>.info` file. They reduce linear, ring, star, and grafted
architectures to a common effective-strand graph, classify network defects,
calculate strand conformation and directional periodic-image shortest paths,
and export native Z1+ input with a companion graph mapping.

The distribution analyzer reports reaction and defect profiles along z and
can use `dump.msd.lammpstrj` for origin-layer-resolved dynamics. See
[`Analysis/README.md`](Analysis/README.md) for commands, definitions, output
columns, and the film/Z1+ boundary caveat.

## Repository layout

- `Generator/pdms_elastomer_generator.cpp`: generic model generator;
- `examples/01_default/`: reproducible current-default sample;
- `examples/02_ring_bifunctional/`: ring strands with two regular reactive sites;
- `examples/03_ring_tetrafunctional/`: ring strands with four random reactive sites;
- `examples/04_star_3arm/`: three-arm, one-center star strands;
- `examples/05_star_4arm/`: four-arm, one-center star strands;
- `examples/06_star_6arm/`: six-arm, two-center star strands;
- `examples/07_star_8arm/`: eight-arm, three-center star strands;
- `examples/08_grafted_comb/`: sparse, comb-like grafted strands;
- `examples/09_grafted_bottlebrush/`: dense bottlebrush strands with a graft on every backbone bead;
- `Generator/pdms_filler_component.hpp`: neutral PDMS filler builder;
- `Analysis/`: topology reduction, Z1+ export, z profiles, and layer dynamics;
- `tests/smoke_test.sh`: generator and analyzer regression checks.

## References

1. D. Zhang et al., “Energy renormalization for temperature transferable
   coarse-graining of silicone polymer,” *Physical Chemistry Chemical Physics*
   **26**, 4541–4554 (2024).
   [https://doi.org/10.1039/d3cp05969c](https://doi.org/10.1039/d3cp05969c)
2. M. Safaripour et al., “Predicting ice adhesion of fluid-containing
   elastomers from surface tension/energy and crosslink density,” *Materials &
   Design* **268**, 116498 (2026).
   [https://doi.org/10.1016/j.matdes.2026.116498](https://doi.org/10.1016/j.matdes.2026.116498)
3. A. P. Thompson et al., “LAMMPS - a flexible simulation tool for
   particle-based materials modeling at the atomic, meso, and continuum
   scales,” *Computer Physics Communications* **271**, 108171 (2022).
   [https://doi.org/10.1016/j.cpc.2021.108171](https://doi.org/10.1016/j.cpc.2021.108171)
