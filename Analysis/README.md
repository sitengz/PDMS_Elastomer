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

