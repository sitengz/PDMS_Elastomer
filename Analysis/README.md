# Network analysis workflow

The analyzers use the final LAMMPS data snapshot, normally
`data.<case>.npt_eq`, together with the matching version-3 `<case>.info` file.
They recognize linear, ring, star, and grafted component-1 architectures.

Build both programs from the repository root:

```bash
make
```

## 1. Topology analyzer

```bash
./bin/topology_analyzer data.CASE.npt_eq CASE.info
```

The reduction is based on actual bond-type-2 reactions in the final snapshot:

- a linear molecule remains one effective strand;
- a ring is cut at its reacted sites, so two reactions make two parallel
  arcs and one reaction is recorded as a dangling loop;
- every star arm is a center-to-arm-end strand;
- reacted grafted side-chain ends are ordered along the graft sequence and
  consecutive ends are connected by their unique intramolecular path;
- short crosslinkers and moderator/crosslinker clusters are collapsed to
  chemical junction nodes.

The analyzer reports conversion, defects, degree distribution, connected
components, parallel edges, self-loops, cycle rank, contour length,
end-to-end vectors, and directional shortest paths from each graph node to
its periodic image. Bulk samples analyze x, y, and z. Films analyze x and y
only. Increase the lifted-cell search if a sparse network needs it:

```bash
./bin/topology_analyzer data.CASE.npt_eq CASE.info \
    --image-search-bound 3
```

Default outputs in `analysis_<case>/` are:

- `topology_report.<case>.txt`;
- `effective_strands.<case>.tsv`;
- `network_nodes.<case>.tsv`;
- `parent_molecules.<case>.tsv`;
- `directional_self_paths.<case>.tsv`;
- `config.<case>.Z1` and `config.<case>.Z1.map.tsv`.

### Z1+ export

The native Z1 file contains the selected effective strands in unwrapped
contour order. Reaction bonds and crosslinker geometry are excluded. The
mapping table retains the original strand, molecule, endpoint atoms, and
reduced junction IDs because native Z1 format does not store chemical graph
connectivity.

By default only active strands between distinct junctions are exported. Other
choices are available for controlled comparisons:

```bash
--z1-selection active
--z1-selection active-and-dangling
--z1-selection all-linearizable
```

By default the Z1 file preserves the physical coordinates and the complete
`Lx Ly Lz` box lengths from the LAMMPS snapshot. Uniform coordinate and box
scaling is opt-in for PPA workflows that impose a maximum bond length, for
example `--z1-max-bond 1.4`. Any applied scale is recorded in the mapping and
report. Native three-line Z1 files do not encode `p p f`; film confinement or
fixed surface objects must be validated in the installed Z1+ version before
interpreting film primitive paths.

## 2. Distribution and layer dynamics

Static z profiles need only the final snapshot:

```bash
./bin/network_profile_analyzer data.CASE.npt_eq CASE.info \
    --bin-width 5
```

`network_z_profile.<case>.tsv` reports reaction-bond, junction,
active-strand, dangling-end, dangling-loop, self-loop, and isolated-parent
distributions. Markers are bond or contour midpoints except dangling ends and
isolated parent centers.

Add the production trajectory for layer-resolved dynamics:

```bash
./bin/network_profile_analyzer data.CASE.npt_eq CASE.info \
    --trajectory dump.msd.lammpstrj --bin-width 5 --frame-stride 1
```

`layer_dynamics.<case>.tsv` groups component-1 beads by their first-frame z
layer and reports x, y, z, in-plane, and total MSD relative to that frame.
Whole-system center-of-mass drift is removed. The generated trajectory's
`x y z ix iy iz` columns are sufficient; `xu yu zu` is also accepted.

The layer MSD is an origin-layer observable, not a time-origin-averaged MSD.
For films, x/y or the in-plane sum is normally the relevant mobility measure,
while z reflects confinement.
