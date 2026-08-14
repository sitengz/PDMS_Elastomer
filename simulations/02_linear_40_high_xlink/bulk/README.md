# Short-linear high-crosslinker-loading reference: bulk

This bulk reference determines the effective-strand end-to-end distance used
to construct the three matching films.

Run from the repository root:

```bash
bash simulations/02_linear_40_high_xlink/bulk/run.sh
```

The generator writes `linear_40_high_xlink_bulk/`. Submit its LAMMPS input,
verify the target conversion and final equilibration, and run
`topology_analyzer` and `basic_network_analyzer` on the final `.npt_eq`
snapshot and matching `.info` file.

Record the mean `Ree` for effective strands in the largest connected network.
Use the same primary definition as the other cases and pass the value to the
three film runners as `BULK_REE_ANGSTROM`.
