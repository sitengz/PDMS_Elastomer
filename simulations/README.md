# Architecture and film-thickness simulation matrix

This directory contains the controlled simulation set for comparing precursor
architecture at nominally constant crosslink density. It is separate from
`examples/`: examples demonstrate generator features, whereas these folders
define the production cases intended for the paper.

## Controlled composition

Cases `01` and `03`--`06` contain 750 four-functional, 32-bead crosslinkers
and 3,000 functional groups on each stoichiometric side. Case `02` doubles
the linear-chain population while halving its contour length, giving 1,500
crosslinkers and 6,000 functional groups on each side. All systems start from
the common placement density `0.05 g/cm^3`, target `0.8 g/cm^3` and 95%
conversion, and contain no moderator or filler. The lower placement density is
required for the 160-bead rings and dense bottlebrushes to generate without
overlap; reaching 95% before compression must be confirmed in the bulk pilots.

| Case | Precursor population | Precursor beads | Total beads |
|:--|:--|--:|--:|
| `01_linear_reference` | 1,500 linear chains, functionality 2 | 120,000 | 144,000 |
| `02_linear_40_high_xlink` | 3,000 linear 40-bead chains, functionality 2 | 120,000 | 168,000 |
| `03_ring_4functional` | 750 rings, functionality 4 | 120,000 | 144,000 |
| `04_star_4arm` | 750 equal four-arm stars | 120,750 | 144,750 |
| `05_comb_4sidechain` | 750 four-side-chain combs | 120,000 | 144,000 |
| `06_bottlebrush_4functional` | 750 four-functional bottlebrushes | 120,000 | 144,000 |

The star intentionally contains 161 beads per precursor so all four arms can
contain exactly 40 beads. Its 0.52% larger total bead count is recorded rather
than hidden by using unequal arms.

## Geometry matrix

Each architecture contains four geometries:

- `bulk`;
- `film_2Ree`;
- `film_4Ree`;
- `film_8Ree`.

Run the bulk system first. After equilibration, use
`basic_network_analyzer` to determine the representative bulk effective-
strand end-to-end distance, `Ree`, in angstrom. Use one stated definition
consistently across all architectures; the recommended primary definition is
the mean `Ree` of effective strands in the largest connected network.

Then generate the three films, for example:

```bash
BULK_REE_ANGSTROM=42.5 \
  bash simulations/01_linear_reference/film_2Ree/run.sh
```

The runner calculates the nominal material thickness as the requested
multiple of `BULK_REE_ANGSTROM`. It refuses to run when the bulk value is
missing, so no unresolved placeholder can silently become a production case.

## Count and statistical scope

The matrix contains six bulk and eighteen film simulations. It uses one
realization per state point. The large number of precursors provides strong
within-system structural sampling, but one realization does not estimate
between-realization uncertainty. If two architectures differ only weakly in
the primary conclusions, add independent seeds for those specific cases
before making a statistical claim.

## Film boundary condition

The current generator creates films between two repulsive `wall/lj126`
boundaries. These are confined films, not vacuum-terminated free-standing
films. The requested thickness is the nominal 300 K wall-force-free material
thickness; the complete box length includes both wall cutoff layers.

## Common workflow

1. Generate and run each bulk system.
2. Verify 95% conversion and equilibration.
3. Run topology and basic-network analysis.
4. Record the chosen bulk `Ree` for each architecture.
5. Generate and run that architecture's `2Ree`, `4Ree`, and `8Ree` films.
6. Apply the same topology, profile, Z1+, and mechanical-analysis definitions
   to every system.
