# Network analysis workflow

The analyzers use the final LAMMPS data snapshot, normally
`data.<case>.npt_eq`, together with the matching version-3 `<case>.info` file.
They recognize linear, ring, star, and grafted component-1 architectures.

Build all three programs from the repository root:

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
contour order. Reaction bonds and crosslinker geometry are excluded. Linear
contours are unchanged. Ring reaction-site endpoints are removed from both
neighboring arc contours, and star center beads are removed from every arm,
so one source bead cannot occur in multiple ring or star Z1 chains.

For each grafted parent, reacted side-chain ends are numbered by increasing
backbone graft position. Path `i < n` starts at reacted end `i`, follows its
side chain and the backbone, and stops at the bead immediately before graft
`i+1`. Path `n` contains only the final reacted side-chain beads and excludes
its backbone grafting bead.
Backbone tails outside the first-to-last reacted range and nonfunctional side
chains are omitted. This partitions the exported grafted contour without
duplicating any source bead, including when only one graft has reacted.

The mapping table retains the original strand, molecule, endpoint atoms,
reduced junction IDs, ordered reacted-site index, next graft atom, excluded
atom IDs, and the exact atom IDs written for each Z1 chain. This makes the
partitioning and overlap checks auditable because native Z1 format does not
store chemical graph connectivity.

The default `network` selection keeps active, dangling, dangling-loop, and
self-loop paths from parent molecules with at least one reacted site. Wholly
isolated parents are omitted. This retains strands that may contribute through
transient entanglements at lower deformation. Other choices are available for
controlled comparisons:

```bash
--z1-selection network
--z1-selection active
--z1-selection active-and-dangling
--z1-selection all-linearizable
```

`active-and-dangling` includes both ordinary dangling strands and one-anchor
ring loops, but not self-loops. `all-linearizable` also includes isolated
linear/star paths. Zero-bead paths produced by adjacent ring reaction sites
are skipped and counted in the topology report.

By default the Z1 file preserves the physical coordinates and the complete
`Lx Ly Lz` box lengths from the LAMMPS snapshot. Uniform coordinate and box
scaling is opt-in for PPA workflows that impose a maximum bond length, for
example `--z1-max-bond 1.4`. Any applied scale is recorded in the mapping and
report. Native three-line Z1 files do not encode `p p f`; film confinement or
fixed surface objects must be validated in the installed Z1+ version before
interpreting film primitive paths.

## 2. Basic network information

Run the static network analyzer after topology reduction:

```bash
./bin/basic_network_analyzer data.CASE.npt_eq CASE.info
```

It uses exactly the same architecture-to-strand reduction as the topology
analyzer. For each effective strand it reports contour length (`Lc`),
end-to-end vector and distance (`Ree`), radius of gyration (`Rg`), gyration
tensor and eigenvalues, shape anisotropy, straightness/tortuosity,
orientation, winding, and contour mass. Statistics are grouped by network
status, parent topology, and membership in the largest connected component.
`Rg` and its tensor are bead-number-weighted over the effective contour.

System-level outputs include snapshot/component validation, mass density,
velocity-derived temperature when a `Velocities` section is present,
conversion, junction and defect counts, degree distribution, cycle rank,
active-strand density, and affine/phantom structural modulus estimates. The
modulus values are topology-based estimates rather than measured mechanical
properties. Their temperature defaults to the final temperature in the info
file and can be overridden with `--modulus-temperature`.

Default outputs in `analysis_<case>/` are:

- `basic_network_report.<case>.txt`;
- `network_statistics.<case>.tsv`;
- `strand_properties.<case>.tsv`;
- `strand_statistics.<case>.tsv`;
- `strand_histograms.<case>.tsv`;
- `junction_properties.<case>.tsv`.

`Lpp` is deliberately written as `NaN` with an unavailable status. A final
snapshot contains no primitive-path contour, and neither `Lc`, `Ree`, nor a
graph shortest path is a valid substitute. A future validated primitive-path
reader can populate this field without changing the table schema.

## 3. Distribution and layer dynamics

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
