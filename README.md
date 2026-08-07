# PDMS Elastomer

This repository generates coarse-grained PDMS elastomer models with an optional
neutral PDMS-chain filler. It combines the V22 and V35 formulations in one
generator and includes the standard post-processing tools used with the
generated `.info` metadata.

## Build

```bash
make
```

The executables are written to `bin/`. A direct generator build is also
possible:

```bash
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic \
    Generator/pdms_elastomer_generator.cpp -o pdms_elastomer_generator
```

## Generate a model

V22 is the default formulation and no filler is added unless both filler
arguments are supplied:

```bash
./bin/pdms_elastomer_generator
./bin/pdms_elastomer_generator --formulation V35
./bin/pdms_elastomer_generator --formulation V22 \
    --filler-length 32 --filler-wt 10
./bin/pdms_elastomer_generator --formulation V35 \
    --filler-length 32 --filler-wt 5 --thickness 238.3
```

Each invocation creates a case-named subdirectory containing four files:

- `data.<case>`: initial LAMMPS data file
- `in.<case>`: equilibration, crosslinking, and MSD-production input
- `submit.<case>.sh`: one-node, 96-task Slurm script
- `<case>.info`: JSON metadata used by the analyzers

Omitting `--thickness` produces a cubic, fully periodic bulk system. Supplying
it fixes `Lz` at the requested film thickness, keeps `Lx = Ly`, uses a
nonperiodic z direction, and writes separate lower and upper wall fixes.

## Formulation defaults

| Setting | V22 | V35 |
|---|---:|---:|
| Network-strand length `N1` | 128 | 384 |
| Network-strand count `M1` | 900 | 306 |
| Crosslinker length `N2` | 32 | 32 |
| Crosslinker functionality | 8 | 4 |
| Reactive-site distribution | random | random |
| Network initial shape | straight | folded | 
| Bond-creation probability | 0.1 | 0.5 |
| Moderator `N4`, `M4` | 5, 6 | 5, 6 |

The crosslinker count is always stoichiometric:

```text
M2 = 2*M1/functionality
```

Functionality must be from 3 through 16 and must divide `2*M1` exactly. All
initial covalent bonds are type 1. Bond type 2 is reserved for bonds created by
LAMMPS during crosslinking.

## Generator options

```text
--formulation V22|V35
--n1 N --m1 M
--n2 N
--n4 N --m4 M
--filler-length N
--filler-wt X
--filler-seed N
--filler-min-separation X
--functionality F
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

For regular crosslinkers, reactive sites are placed at 1, 3, 5, and so on.
Random placement samples distinct sites reproducibly from `--crosslink-seed`.
The neutral PDMS filler is component 3 and is initially packed in the central
40% of the box while avoiding the network, crosslinkers, moderators, other
filler chains, and the box boundaries.

## Simulation template

The generated input first relaxes the low-density model at 800 K with a
repulsive Lennard-Jones cutoff at the potential minimum. Compression and
crosslinking then proceed together, followed by cooling and final equilibration
at 300 K. The final stage is an independent 1,000,000-step 300 K NVT run with
1,001 expected trajectory frames (`x y z ix iy iz`) for MSD analysis.

The generated files are templates for the Iowa State Nova environment. Review
the module names, partition, memory, wall time, and email before submitting on a
different cluster.

## Analysis

See [Analysis/README.md](Analysis/README.md) for the z-density profile, final
snapshot/network report, and filler/all-system MSD workflows.

## Repository layout

- `Generator/pdms_elastomer_generator.cpp`: unified V22/V35 command-line tool
- `Generator/pdms_filler_component.hpp`: boundary-safe neutral PDMS filler builder
- `Analysis/z_profile.cpp`: overlapping-window density and composition profiles
- `Analysis/final_snapshot_analyzer.cpp`: conversion, topology, and end-to-end report
- `Analysis/msd_analyzer.cpp`: bead or molecular-COM MSD and diffusion fitting
- `tests/smoke_test.sh`: compile and output-namespace regression checks

