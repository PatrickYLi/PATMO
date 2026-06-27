# README

This folder is a volcanic-eruption variant of `modern_sulfur_cycle`.

It keeps the same reaction network, atmospheric profile, solar spectrum, rainout,
and aerosol inputs, then enables `patmo_volc.f90` through `volcano_events.dat`.

## Event File

`volcano_events.dat` uses one event per row:

```text
start_day duration_day center_km sigma_km so2_column_flux ash_tau_550 ash_lifetime_day ash_settling_cm_s ash_wavelength_exp
```

- `so2_column_flux`: SO2 column source in molecules cm-2 s-1, distributed over a Gaussian plume centered at `center_km`.
- `ash_tau_550`: effective 550 nm volcanic ash optical depth. Increase this value to test thicker ash shading.
- `ash_lifetime_day`: e-folding time after the eruption duration ends. Use a negative value for no decay.
- `ash_settling_cm_s`: downward settling speed of the ash plume center.
- `ash_wavelength_exp`: wavelength scaling, with optical depth proportional to `(lambda / 550 nm)^(-ash_wavelength_exp)`.

## Usage

From the repository root:

```bash
./tests/modern_sulfur_volcano/compile_modern_sulfur_volcano.sh
```

Enter `modern_sulfur_volcano` when prompted. Then compile in `./build`:

```bash
make
```

The run writes the usual sulfur-cycle outputs plus:

- `volcano_initial_state.dat`
- `volcano_final_state.dat`

These files contain altitude, SO2 source rate, local ash optical depth, and cumulative ash optical depth.
