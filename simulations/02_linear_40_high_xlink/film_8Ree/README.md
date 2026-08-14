# Short-linear high-crosslinker-loading reference: 8Ree film

This geometry uses the nominal material thickness

```text
H = 8 * Ree_bulk
```

where `Ree_bulk` comes from the matching equilibrated bulk network. The value
is intentionally not hard-coded.

Run from the repository root after replacing `42.5` with the measured bulk
value in angstrom:

```bash
BULK_REE_ANGSTROM=42.5 \
  bash simulations/02_linear_40_high_xlink/film_8Ree/run.sh
```

The runner validates `Ree_bulk`, calculates `H`, and passes it through
`--thickness`. The generator writes `linear_40_high_xlink_film_8Ree/`.

The current geometry uses two repulsive walls normal to `z` and periodic
boundaries in `x` and `y`. Report it as a confined film unless a later
generator revision introduces vacuum-terminated free surfaces.
