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

The independent mechanical-test generator can be built with:

```bash
make bin/tensile_test_generator
```

The independent extended layer-dynamics generator can be built with:

```bash
make bin/layer_dynamics_generator
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

Omitting `--thickness` produces a cubic, fully periodic bulk system. For a
film, `--thickness H` specifies the nominal wall-free material thickness at
300 K. The generator sets

```text
Lz = H + 2*cold_wall_cutoff
```

and calculates `Lx = Ly` from the nominal material volume `Lx*Ly*H`. The full
box uses a nonperiodic z direction and separate lower and upper repulsive wall
fixes. The final equilibrated density profile remains the preferred way to
measure the realized material thickness because the wall is a soft potential.

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
--target-conversion X
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

The generated LAMMPS input applies the same preparation sequence to every
strand topology. Molecular architecture changes the initial bonded graph and
reactive-site locations, but not the standard minimization, curing, cooling,
or final-equilibration stages. The timestep is 5 fs, so 1,000,000 steps
correspond to 5 ns.

### Standard conversion-controlled workflow

The production workflow used when `--target-conversion X` is supplied is:

```text
Initial configuration: data.<case> + <case>.info
  user-selected topology, composition, stoichiometry, density and geometry
                                |
                                v
800 K force field and boundary conditions
  bulk: periodic x/y/z
  film: periodic x/y, fixed z, two repulsive walls
                                |
                                v
Conjugate-gradient energy minimization at fixed box dimensions
  maximum bead displacement = 0.1 A per minimizer iteration
                                |
                                v
Assign Gaussian velocities at 800 K
  reproducible seed; zero net linear and angular momentum
                                |
                                v
800 K NVT equilibration: 1,000,000 steps (5 ns)
                                |
                                +--> data.<case>.rep_800
                                |
                                v
Compress without reactions to 0.5 g/cm3: 1,000,000 steps (5 ns)
  bulk: deform x/y/z
  film: deform x/y while Lz remains fixed
                                |
                                v
Conversion-controlled curing at 800 K and 0.5 g/cm3
  fix bond/create probability = 0.5
  stop at the requested new-bond count or after 5,000,000 steps (25 ns)
                                |
                                v
Remove bond/create and halt fixes; network connectivity is now frozen
                                |
                                v
Compress from 0.5 g/cm3 to the requested target density
  1,000,000 steps (5 ns); bulk x/y/z or film x/y only
                                |
                                v
Post-cure 800 K NVT equilibration: 1,000,000 steps (5 ns)
                                |
                                +--> data.<case>.xlink_800
                                |
                                v
Restore the 300 K PDMS pair, reacted-bond and film-wall parameters
                                |
                                v
Cool from 800 K to 300 K at 1 atm: 1,000,000 steps (5 ns)
  bulk: isotropic NPT
  film: xy-coupled NPT with fixed Lz
                                |
                                +--> data.<case>.300
                                |
                                v
Final 300 K, 1 atm equilibration: 1,000,000 steps (5 ns)
  bulk: isotropic NPT
  film: xy-coupled NPT with fixed Lz
                                |
                                v
Final equilibrated network: data.<case>.npt_eq
  common input for topology, profiles, Z1+, dynamics and tensile tests
```

The fixed part of this protocol contains 6,000,000 steps. Including the curing
hold, the final `.npt_eq` state is reached after at most 11,000,000 steps, or
55 ns. The actual duration is shorter when the conversion target is detected
before the 5,000,000-step curing limit.

The bulk and film branches differ only in boundary and box control:

```text
Bulk system                         Film system
-----------                         -----------
boundary p p p                      boundary p p f
no confining walls                  repulsive walls at both z edges
x/y/z compression                   x/y compression; Lz fixed
isotropic NPT cooling               xy-coupled NPT cooling; Lz fixed
isotropic final NPT                 xy-coupled final NPT; Lz fixed
```

For a film, the requested thickness is the nominal 300 K wall-force-free
material thickness, not the complete box length. The wall-confined
`data.<case>.npt_eq` state is the end of network preparation. Moving the walls
away from the material for effectively free-surface dynamics is a separate,
nonreactive workflow described under [Extended layer dynamics](#extended-layer-dynamics).

The minimization uses the 800 K force field at fixed box dimensions. Its
conservative displacement limit relaxes severe initial overlaps, and film-wall
energies are included during this stage. If no conversion target is supplied,
the generator retains the original time-controlled route: initial relaxation,
compression with active crosslinking, a fixed-duration curing hold, cooling,
and final equilibration.

### Conversion control

Supplying `--target-conversion X` enables conversion-controlled curing. The
generator calculates

```text
maximum_new_bonds = min(strand functional groups,
                        crosslinker functional groups)
target_new_bonds  = floor(X/100 * maximum_new_bonds)
```

Moderator groups remain extra and are excluded from both stoichiometry and
this target. After the low-density 800 K equilibration, the system is compressed
without reactions to 0.5 g/cm3 over 1,000,000 steps. At that fixed curing
density, `fix bond/create` uses a default probability of 0.5 and a `fix halt`
condition checks the cumulative new-bond count every timestep during a curing
hold of at most 5,000,000 steps. When the target is detected, the hold ends and
both the halt and bond-creation fixes are removed. The cured network is then
compressed to the requested target density over 1,000,000 steps and
equilibrated for another 1,000,000 steps at 800 K without further reactions.
Finally, the timestep is reset before the parameter-change and cooling stages.
The same target setting can be written in a config file as
`target_conversion = 95`.

Because `fix bond/create` may form several independent bonds on its final
timestep, the realized count can exceed the target slightly. The 5,000,000-step
curing hold at 0.5 g/cm3 is an upper bound rather than a guarantee: if
geometric constraints produce a plateau below the requested conversion, curing
ends at that limit. The `.info` file records the target basis, compression and
curing schedule, integer bond target, upper bounds, and possible final-step
overshoot.

Every workflow finishes by writing an independent 1,000,000-step NVT
trajectory for MSD analysis.

The built-in MSD trajectory is produced only after `data.<case>.npt_eq` has
been written and is not part of the network-preparation time quoted above.
The longer, layer-resolved dynamics workflow is generated independently from
the same `.npt_eq` and `.info` pair.

The generated Slurm file is a template for the Iowa State Nova environment.
Review its modules, partition, memory, wall time, and email before use on a
different cluster.

## Tensile tests

The tensile-test generator starts from the frozen network in
`data.<case>.npt_eq` and the matching version-3 `<case>.info`. It does not
minimize, create bonds, or change network connectivity. Direct use is:

```bash
./bin/tensile_test_generator data.CASE.npt_eq CASE.info
```

The default `auto` mode creates independent `x`, `y`, and `z` tests for a
bulk system and the two in-plane `x` and `y` tests for a film. Each production
geometry also has `run_tensile.sh` beside its original `run.sh`. Once the
original simulation has completed, the usual command has no arguments:

```bash
bash simulations/02_linear_40_high_xlink/film_4Ree/run_tensile.sh
```

The runner derives the case name from `model.conf`, locates the normal
`<case>/data.<case>.npt_eq` and `<case>/<case>.info` outputs, and creates
`<case>/tensile/<direction>/`. Each direction contains a LAMMPS input,
versioned tensile `.info`, and Nova submission script. The run performs 5M
steps of 300 K equilibration followed by engineering-rate deformation at
`2e-9 fs^-1` to 20% strain. It writes block-averaged stress, stretch, box,
and volume data without a trajectory dump. Film stresses are normalized by
the nominal material volume `Lx*Ly*H`, not the full wall-containing box.

Nonstandard source locations can be supplied without editing the runner:

```bash
NPT_EQ_FILE=/path/data.CASE.npt_eq \
INFO_FILE=/path/CASE.info \
TENSILE_OUTPUT_DIR=/path/tensile \
  bash simulations/ARCHITECTURE/GEOMETRY/run_tensile.sh
```

## Extended layer dynamics

The original model workflow retains its established 1M-step MSD trajectory.
For the longer layer-resolved dynamics study, use the independent
`run_layer_dynamics.sh` beside `run.sh` and `run_tensile.sh`:

```bash
bash simulations/02_linear_40_high_xlink/film_4Ree/run_layer_dynamics.sh
```

The runner automatically reads the same final `.npt_eq` and version-3
`.info` pair. It creates `<case>/layer_dynamics/` containing a 300 K NVT
input, a versioned layer-dynamics `.info`, and a Nova submission script. For
a film, the material-adjacent source walls are not recreated and the
nonperiodic z box is expanded by 50 A at each boundary without remapping
atoms. Matching 300 K repulsive guard walls are placed only at the expanded
box edges. The film therefore retains two free material surfaces and wide
vacuum buffers, while the remote walls prevent detached chains from leaving
the box. Bulk boxes remain unchanged and periodic.

The source `.npt_eq` velocities are retained and network connectivity is
unchanged. After wall removal and expansion, the system equilibrates for 1M
steps (5 ns) and writes `data.<case>.free_surface_eq`. A 20,000-step (100 ps)
Debye-Waller segment is then written every 20 steps (0.1 ps) to
`dump.debye_waller.<case>.lammpstrj`. After that segment, the dynamics clock
and MSD origin are reset for an independent 5M-step (25 ns) production in
`dump.layer_dynamics.<case>.lammpstrj`. MSD frames are written every 1,000
steps (5 ps) through the first production 1M steps and every 5,000 steps
(25 ps) through the remaining 4M steps. The two expected trajectories contain
1,001 and 1,801 frames, respectively, with `x y z ix iy iz`. For example,
to request 20M total steps while retaining the same sampling schedule:

```bash
bash simulations/02_linear_40_high_xlink/film_4Ree/run_layer_dynamics.sh \
  --production-steps 20000000
```

Nonstandard source locations use `NPT_EQ_FILE` and `INFO_FILE`, as in the
tensile runner. Set `LAYER_DYNAMICS_OUTPUT_DIR` to override the output folder.

## Analysis

The topology-first analyzers read `data.<case>.npt_eq` and the matching
version-3 `<case>.info` file. They reduce linear, ring, star, and grafted
architectures to a common effective-strand graph, classify network defects,
calculate strand conformation and directional periodic-image shortest paths,
and export native Z1+ input with a companion graph mapping.

The Z1+ export keeps network dangling and self-loop paths by default. It
removes shared reaction-site beads from ring arcs and shared center beads from
star arms. Grafted paths are partitioned in reacted-graft order: each path
contains one functional side chain and the following backbone interval, while
the last path contains only the final side-chain beads.

The second-stage basic network analyzer reports per-strand `Ree`, `Rg`,
contour and shape metrics, grouped distributions, density, conversion,
connectivity, and structural elasticity estimates. `Lpp` remains explicitly
unavailable until a validated primitive-path result is supplied. The third
analyzer reports component density, local conversion, reaction, defect,
contour, conformation, and orientation profiles along z. It also writes
folded wall-distance and wall/core summary tables, optionally reads
`Z1+SP.dat` for kink and primitive-path profiles, and can use separate
high-frequency Debye-Waller and long-time MSD dumps for
origin-layer-resolved dynamics. A fourth fixed-lag analyzer uses every valid
time origin to extract the standard 10 ns MSD and Einstein ratio together with
the 10 ps Debye-Waller displacement, preserving directional and layer-resolved
values for uncertainty analysis. See
[`Analysis/README.md`](Analysis/README.md) for commands, definitions, output
columns, and the film/Z1+ boundary caveat.

## Repository layout

- `Generator/pdms_elastomer_generator.cpp`: generic model generator;
- `Generator/tensile_test_generator.cpp`: post-equilibration tensile-test
  generator for bulk and in-plane film loading;
- `Generator/layer_dynamics_generator.cpp`: extended free-surface trajectory
  generator for origin-layer-resolved dynamics;
- `examples/01_default/`: reproducible current-default sample;
- `examples/02_ring_bifunctional/`: ring strands with two regular reactive sites;
- `examples/03_ring_tetrafunctional/`: ring strands with four random reactive sites;
- `examples/04_star_3arm/`: three-arm, one-center star strands;
- `examples/05_star_4arm/`: four-arm, one-center star strands;
- `examples/06_star_6arm/`: six-arm, two-center star strands;
- `examples/07_star_8arm/`: eight-arm, three-center star strands;
- `examples/08_grafted_comb/`: sparse, comb-like grafted strands;
- `examples/09_grafted_bottlebrush/`: dense bottlebrush strands with a graft on every backbone bead;
- `examples/10_default_film/`: film counterpart of the default linear sample, with a thickness placeholder;
- `simulations/`: six-case production matrix containing bulk and
  `2Ree`, `4Ree`, and `8Ree` confined-film cases for the manuscript study;
- `Generator/pdms_filler_component.hpp`: neutral PDMS filler builder;
- `Analysis/`: topology reduction, static network properties, Z1+ export,
  z profiles, and layer dynamics;
- `tests/smoke_test.sh`: model-generator and analyzer regression checks;
- `tests/tensile_generator_test.sh`: bulk, film, and auto-detection tensile
  generator checks;
- `tests/layer_dynamics_generator_test.sh`: extended bulk/film dynamics and
  auto-detection checks.

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
