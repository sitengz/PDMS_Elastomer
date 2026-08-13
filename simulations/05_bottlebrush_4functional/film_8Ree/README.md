# Four-functional bottlebrushes: 8Ree film

This geometry uses the nominal material thickness

```text
H = 8 * Ree_bulk
```

where `Ree_bulk` comes from the matching equilibrated bulk network. The
value is intentionally not hard-coded.

Run from the repository root after replacing `42.5` with the measured bulk
value in angstrom:

```bash
BULK_REE_ANGSTROM=42.5 \
  bash simulations/05_bottlebrush_4functional/film_8Ree/run.sh
```

The runner validates `Ree_bulk`, calculates `H`, and passes it through
`--thickness`. The generator writes `bottlebrush_4functional_film_8Ree/`.

The current geometry uses two repulsive walls normal to `z` and periodic
boundaries in `x` and `y`. Report it as a confined film unless a later
generator revision introduces vacuum-terminated free surfaces.
