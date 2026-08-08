# Analysis workflow

Build all analyzers from the repository root with `make`, or compile an
individual source using C++17. Each analyzer reads the generated JSON `.info`
file so component molecule ranges and simulation timing remain tied to the
specific sample.

## 1. Distribution along z

```bash
./bin/z_profile data.CASE.npt_eq CASE.info \
    --sample-spacing 1 --window-width 5
```

`--sample-spacing` controls the distance between reported points, while
`--window-width` controls the counting window. A window wider than the spacing
therefore produces overlapping samples and a smoother profile without reducing
the number of points. The output columns are numeric and suitable for MATLAB.
For samples without filler, filler-enrichment quantities are written as `NaN`.

Default outputs are created in a folder named after the sample:

- `z_profile.<case>.dat`
- `z_regions.<case>.dat`
- `z_report.<case>.txt`

## 2. Final snapshot, crosslinking, and end-to-end distance

```bash
./bin/final_snapshot_analyzer data.CASE.npt_eq CASE.info
```

The report checks composition and topology, calculates final density and any
velocity-derived temperature, reports crosslink conversion and network
connectivity, and summarizes network-strand end-to-end distances. The
headerless per-strand table can be read directly by MATLAB.

Default outputs:

- `<case>/final_snapshot_report.<case>.txt`
- `<case>/strand_end_to_end.<case>.dat`

## 3. MSD and diffusion

```bash
./bin/msd_analyzer dump.msd.lammpstrj CASE.info
```

The default observable is the time-averaged, mass-weighted center-of-mass MSD
of each filler chain with whole-system center-of-mass drift removed. Other
choices include the complete system, individual beads, and raw MSD:

```bash
./bin/msd_analyzer dump.msd.lammpstrj CASE.info \
    --selection all --particle beads --averaging raw
```

For three-dimensional bulk diffusion, the fitted relation is `D = slope/6`.
The default time-averaged fit uses 50–90% of the lag-time range, dropping the
sparsely sampled last 10%. The zero-lag row is omitted from the numeric output.
For films, inspect directional MSD components because confinement makes a
single isotropic 3D coefficient potentially misleading.

Default outputs:

- `<case>/msd.<case>.dat`
- `<case>/msd_report.<case>.txt`

Run any executable with `--help` for every available option and column detail.

## 4. Binary elastomer–filler morphology

`binary_morphology_analyzer` treats each sample as a two-population system:

- filler: every bead in component 3;
- elastomer matrix: network strands, crosslinkers, and moderators.

Membership is determined from molecule-ID ranges in the matching `.info`
file, not from atom type. This is required because the PDMS filler and neutral
matrix beads share force-field types.

The analyzer retains three complementary descriptors:

1. a voxelized local filler-fraction field and finite-count-corrected
   segregation index;
2. the binary concentration structure factor `Scc(q)`; and
3. a filler-chain contact graph with connected-component and graph-clustering
   statistics.

These are morphology descriptors, not an automatic thermodynamic phase
classifier. Conclusions should be based on their agreement, sensitivity to
grid spacing/contact cutoff, visual inspection, and—when available—stability
over time.

### Run

The final equilibrated data file is sufficient for static morphology:

```bash
./bin/binary_morphology_analyzer data.CASE.npt_eq CASE.info
```

The generated MSD trajectory can instead show morphology evolution:

```bash
./bin/binary_morphology_analyzer dump.msd.lammpstrj CASE.info \
    --frame-stride 10
```

Important controls are:

```text
--grid-spacing X      target voxel spacing, default 8 A
--contact-cutoff X    interchain filler-bead cutoff, default 8 A
--q-max X             maximum analyzed q, default 80% of grid Nyquist
```

The FFT requires power-of-two grid dimensions, so realized voxel dimensions
can differ from the requested spacing and are recorded in the report. Use the
same requested spacing and contact cutoff when comparing samples.

Default outputs are written to the case-named folder:

- `morphology_metrics.<case>.dat`
- `morphology_structure_factor.<case>.dat`
- `morphology_field_final.<case>.dat`
- `morphology_report.<case>.txt`

### Local composition

For occupied voxel `i`, the local filler bead fraction is

```text
phi_i = N_filler,i / (N_filler,i + N_matrix,i)
```

The analyzer calculates the bead-count-weighted variance of `phi_i`, subtracts
the fixed-composition random-label variance, and normalizes by the variance of
complete segregation. An index near zero is consistent with random labels at
the selected voxel scale; positive values indicate excess local segregation;
values approaching one indicate nearly pure local cells. Slightly negative
values can arise from finite-sample fluctuations.

Bulk samples report a 3D index. Film samples also report an XY index after
integrating through z; use the XY result for lateral domains and `z_profile`
for wall-normal layering.

The 10-column final field contains:

```text
ix iy iz x_center y_center z_center filler_beads matrix_beads total_beads filler_fraction
```

Visualize it in MATLAB with:

```matlab
plot_binary_field_3d('morphology_field_final.CASE.dat', ...
    'MinFillerFraction', 0.1, 'MaxMarkerArea', 180)
```

### Concentration structure factor

The normalized binary concentration mode is

```text
Scc(q) = |x_matrix*rho_filler(q) - x_filler*rho_matrix(q)|^2
         / (N*x_filler*x_matrix)
```

A random-label reference is approximately one. Enhanced low-q intensity is
consistent with large-scale composition fluctuations. The preliminary length
scale is `2*pi/q_peak`; if the peak is in the first q shell, this value is
box-limited and should be treated as a lower bound. Bulk systems use radial 3D
shells, while films use in-plane `q_xy` after integration through z.

### Filler-chain contact graph

Two filler chains are connected when any pair of their beads lies within the
contact cutoff under the applicable periodic boundaries. Intrachain contacts
are excluded. The analyzer reports:

- number of connected components;
- largest component in chains and as a fraction of all filler chains;
- unique chain–chain contacts and mean contacts per chain;
- graph transitivity; and
- mean local clustering coefficient over nodes with at least two neighbors.

For chain `i`, the local clustering coefficient is

```text
C_i = closed neighbor pairs / all neighbor pairs
```

A large connected component measures long-range connectivity. A high
clustering coefficient instead measures local triangle closure; a long bundle
can have a large component but low clustering. Both depend on the operational
contact cutoff and should be reported with it.

The 22-column metrics table contains, in order:

```text
frame_index timestep time_ns filler_markers matrix_markers filler_fraction
segregation_3d observed_variance random_variance segregation_xy
low_q_mean q_peak domain_size peak_scc peak_at_lowest_q
cluster_count largest_cluster_chains largest_cluster_fraction
unique_chain_contacts mean_contacts_per_chain graph_transitivity
mean_local_clustering
```

The structure-factor table contains `frame_index timestep time_ns q modes
Scc`. All tables are headerless numeric output for direct MATLAB import; the
text report records definitions and parameter values.

### References

1. M. D. Lefebvre, M. Olvera de la Cruz, and K. R. Shull, “Phase
   Segregation in Gradient Copolymer Melts,” *Macromolecules* **37**,
   1118–1123 (2004). [doi:10.1021/ma035141a](https://doi.org/10.1021/ma035141a)
2. A. B. Bhatia and D. E. Thornton, “Structural Aspects of the Electrical
   Resistivity of Binary Alloys,” *Physical Review B* **2**, 3004–3012
   (1970). [doi:10.1103/PhysRevB.2.3004](https://doi.org/10.1103/PhysRevB.2.3004)
3. B. Chu and B. S. Hsiao, “Small-Angle X-ray Scattering of Polymers,”
   *Chemical Reviews* **101**, 1727–1762 (2001).
   [doi:10.1021/cr9900376](https://doi.org/10.1021/cr9900376)
4. K. S. Schweizer and J. G. Curro, “Microscopic theory of the structure,
   thermodynamics, and apparent chi parameter of polymer blends,” *Physical
   Review Letters* **60**, 809–812 (1988).
   [doi:10.1103/PhysRevLett.60.809](https://doi.org/10.1103/PhysRevLett.60.809)
5. D. J. Watts and S. H. Strogatz, “Collective dynamics of small-world
   networks,” *Nature* **393**, 440–442 (1998).
   [doi:10.1038/30918](https://doi.org/10.1038/30918)
6. R. Yang et al., “Graph Theoretical Description of Phase Transitions in
   Complex Multiscale Phases with Supramolecular Assemblies,” *Advanced
   Science* **11**, e2402464 (2024).
   [doi:10.1002/advs.202402464](https://doi.org/10.1002/advs.202402464)

The papers motivate local-composition, scattering, and graph descriptors. The
finite-particle correction, voxel implementation, and chain-contact graph are
explicit reproducible choices of this software rather than equations claimed
from those sources.
