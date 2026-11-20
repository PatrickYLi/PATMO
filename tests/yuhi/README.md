# Yuhi case inputs and generation

From this case directory:

```bash
./compile_yuhi.sh
make -C ../../build
```

The script generates model sources and copies input data into the project's
`build` directory. It does not run the model integration.

## External reaction rates

In `reaction_network.xlsx`, use `read:<filename>` in the `k` column for a rate
profile stored in this case directory. Ordinary Fortran rate expressions remain
unchanged. The filename is explicit, so changing the order of reactions does not
change which file supplies a reaction's rate.

| Reaction | k |
|---|---|
| 33COS + OH -> CO2 + 33SH | `read:33COS_OH_rate.txt` |
| 34COS + OH -> CO2 + 34SH | `read:34COS_OH_rate.txt` |
| 36COS + OH -> CO2 + 36SH | `read:36COS_OH_rate.txt` |
| 33COS + O -> CO + 33SO | `read:33COS_O_rate.txt` |
| 34COS + O -> CO + 34SO | `read:34COS_O_rate.txt` |
| 36COS + O -> CO + 36SO | `read:36COS_O_rate.txt` |
| 32SO2 + OH -> 32HSO3 | `read:32SO2_OH_rate.txt` |
| 33SO2 + OH -> 33HSO3 | `read:33SO2_OH_rate.txt` |
| 34SO2 + OH -> 34HSO3 | `read:34SO2_OH_rate.txt` |
| 36SO2 + OH -> 36HSO3 | `read:36SO2_OH_rate.txt` |
| 32H2SO4 -> 32SO2 + OH + OH | `read:32H2SO4_photorate.txt` |
| 33H2SO4 -> 33SO2 + OH + OH | `read:33H2SO4_photorate.txt` |
| 34H2SO4 -> 34SO2 + OH + OH | `read:34H2SO4_photorate.txt` |
| 36H2SO4 -> 36SO2 + OH + OH | `read:36H2SO4_photorate.txt` |

Each file contains one finite, nonnegative numeric value per atmospheric layer,
without a header, in the same layer order as the model. The number of values
must equal `cellsNumber` in the generated `options.opt` (currently 60). The
script preserves the values and their units; it does not interpolate or reverse
these rate profiles.

The script adds referenced files to `copylist.pcp` and generates file opens
before the cell loop in `computeRates`, one read per layer, and closes after the
loop. Each call starts reading from the beginning, matching the v4.2 behavior.
Bare `read` entries are rejected because they do not identify a data file.

## Isotope aerosol templates

`templates/patmo_ode.f90` uses the four isotope aerosol formation blocks from
`yuhi_v4.2/patmo_ode.f90`, retaining layers 13–34 and the original threshold
conditions and equations. `templates/patmo_parameters.f90` supplies the matching
`va32/33/34/36` and `pa32/33/34/36` declarations from v4.2.

The compile script selects these case templates during generation. Changes
belong in the workbook or these templates; files under `build` are regenerated.
