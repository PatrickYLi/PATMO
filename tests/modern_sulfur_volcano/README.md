# README

This folder is a volcanic-eruption variant of `modern_sulfur_cycle`.

It keeps the same reaction network, atmospheric profile, solar spectrum, rainout,
and aerosol inputs, then enables `patmo_volc.f90` through `volcano_events.dat`.

## Event File

`volcano_events.dat` uses one active event per non-comment row. The preferred
format is `key=value`, for example:

```text
event_id=thin_ash start_day=0d0 duration_day=20d0 plume_center_km=25d0 plume_sigma_km=2.5d0 so2_flux_cm2_s=5d11 ash_tau_550=0.5d0 ash_lifetime_day=120d0 ash_settling_cm_s=0.2d0 ash_lambda_exponent=1d0
```

- `event_id`: free label for your own bookkeeping.
- `start_day`: eruption start time from model t=0 in days.
- `duration_day`: active SO2 injection and full ash-load interval in days.
- `plume_center_km`: center altitude of the SO2/ash plume in km.
- `plume_sigma_km`: Gaussian vertical width of the plume in km.
- `so2_flux_cm2_s`: SO2 column source in molecules cm-2 s-1, distributed over a Gaussian plume centered at `plume_center_km`.
- `ash_tau_550`: effective 550 nm volcanic ash optical depth. Increase this value to test thicker ash shading.
- `ash_lifetime_day`: e-folding time after the eruption duration ends. Use a negative value for no decay.
- `ash_settling_cm_s`: downward settling speed of the ash plume center.
- `ash_lambda_exponent`: wavelength scaling, with optical depth proportional to `(lambda / 550 nm)^(-ash_lambda_exponent)`.

The Fortran reader still accepts the older numeric-column format:

```text
start_day duration_day center_km sigma_km so2_column_flux ash_tau_550 ash_lifetime_day ash_settling_cm_s ash_wavelength_exp
```

For readability, prefer the named `key=value` form in new experiments.

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
