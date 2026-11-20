#!/usr/bin/env bash
set -euo pipefail

# --- Locate this case and the PATMO root independently of the caller's cwd ---
BASE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${BASE_DIR}/../.." && pwd)"
USER_SUBDIR="$(basename -- "$BASE_DIR")"

if [[ ! -f "${PROJECT_ROOT}/patmo" ]]; then
  echo "PATMO root not found (missing patmo): $PROJECT_ROOT" >&2
  exit 1
fi

# PATMO resolves its data, source, and build paths relative to the project root.
cd -- "$PROJECT_ROOT"

# PATMO uses the database in the project root's data/ directory.
if [[ ! -d "$PROJECT_ROOT/data" ]]; then
  echo "PATMO database directory not found: $PROJECT_ROOT/data" >&2
  exit 1
fi

for REQUIRED_DATA in \
  parseReference/reference.dat \
  thermochemistry/thermochemistry.dat \
  kida/kida.uva.2014.dat \
  kida/kida_reac_tb_2015-03-19_1.csv; do
  if [[ ! -r "$PROJECT_ROOT/data/$REQUIRED_DATA" ]]; then
    echo "Required PATMO database file missing or unreadable: $PROJECT_ROOT/data/$REQUIRED_DATA" >&2
    exit 1
  fi
done
if [[ ! -d "$PROJECT_ROOT/data/xsecs_swri" ]]; then
  echo "PATMO cross-section directory missing: $PROJECT_ROOT/data/xsecs_swri" >&2
  exit 1
fi

INFILE="${BASE_DIR}/reaction_network.xlsx"
OUTFILE="${BASE_DIR}/reaction_network.ntw"

echo "[*] Using base folder : $BASE_DIR"
echo "[*] Input Excel       : $INFILE"
echo "[*] Output NTW        : $OUTFILE"

# -------------------------------
# ------ Dependency check -------
# -------------------------------
need_cmd() { command -v "$1" >/dev/null 2>&1; }

echo "[*] Checking python3 ..."
if ! need_cmd python3; then
  echo "python3 not found. Please install it: sudo apt install -y python3"
  exit 1
fi

echo "[*] Checking pip ..."
if ! python3 -m pip --version >/dev/null 2>&1; then
  echo "[*] Installing python3-pip ..."
  if need_cmd sudo; then
    sudo apt update && sudo apt install -y python3-pip
  else
    echo "pip not available and no sudo permission. Please install python3-pip manually." >&2
    exit 1
  fi
fi

echo "[*] Checking Python dependencies (pandas, openpyxl) ..."
if ! python3 - <<'PY' >/dev/null 2>&1
try:
    import pandas, openpyxl  # noqa
except Exception:
    raise SystemExit(1)
PY
then
  echo "[*] Installing pandas and openpyxl ..."
  if ! python3 -m pip install --user --upgrade pandas openpyxl >/dev/null 2>&1; then
    echo "[*] Retrying install with --break-system-packages ..."
    python3 -m pip install --user --upgrade --break-system-packages pandas openpyxl
  fi
fi

# -------------------------------
# Reaction Network Python Converter
# -------------------------------
python3 - "$INFILE" "$OUTFILE" <<'PY'
# -*- coding: utf-8 -*-
import sys
from pathlib import Path
import pandas as pd

EXPECTED_HEADER_LEADING = "@format:idx"
INTERNAL_COLUMNS = ["IDX","R1","R2","R3","P1","P2","P3","RATE"]
EXPORT_HEADER = ["@format:idx","R","R","R","P","P","P","rate"]
HEADER_COMMENT = "#@var:T=Tgas"

def find_columns(df):
    df.columns = [c.strip() if isinstance(c, str) else c for c in df.columns]
    if EXPECTED_HEADER_LEADING not in df.columns:
        raise ValueError(f"Required column '{EXPECTED_HEADER_LEADING}' not found.")
    idx_col = EXPECTED_HEADER_LEADING
    r_cols = [c for c in df.columns if isinstance(c,str) and (c=="R" or c.startswith("R."))]
    p_cols = [c for c in df.columns if isinstance(c,str) and (c=="P" or c.startswith("P."))]
    if "k" in df.columns:
        rate_col = "k"
    elif "rate" in df.columns:
        rate_col = "rate"
    else:
        raise ValueError("Rate column 'k' or 'rate' not found.")
    r_cols = (r_cols + [None, None, None])[:3]
    p_cols = (p_cols + [None, None, None])[:3]
    return idx_col, r_cols, p_cols, rate_col

def s(series):
    series = series.astype("string").fillna("")
    return series.apply(lambda x: x.strip() if isinstance(x,str) else x)

def build_output_df(df, idx_col, r_cols, p_cols, rate_col):
    out = pd.DataFrame()
    out["IDX"]  = s(df[idx_col])
    out["R1"]   = s(df[r_cols[0]]) if r_cols[0] else ""
    out["R2"]   = s(df[r_cols[1]]) if r_cols[1] else ""
    out["R3"]   = s(df[r_cols[2]]) if r_cols[2] else ""
    out["P1"]   = s(df[p_cols[0]]) if p_cols[0] else ""
    out["P2"]   = s(df[p_cols[1]]) if p_cols[1] else ""
    out["P3"]   = s(df[p_cols[2]]) if p_cols[2] else ""
    out["RATE"] = s(df[rate_col])
    return out[INTERNAL_COLUMNS]

def main():
    if len(sys.argv) < 3:
        raise SystemExit("Usage: python3 - <in.xlsx> <out.ntw>")
    in_path = Path(sys.argv[1]).expanduser().resolve()
    out_path = Path(sys.argv[2]).expanduser().resolve()
    if not in_path.exists():
        raise SystemExit(f"Input file not found: {in_path}")
    df = pd.read_excel(in_path, dtype="string", engine="openpyxl")
    if df.empty:
        raise SystemExit("Excel file is empty.")
    idx_col, r_cols, p_cols, rate_col = find_columns(df)
    out_df = build_output_df(df, idx_col, r_cols, p_cols, rate_col)
    tmp_path = out_path.with_suffix(".tmp")
    try:
        out_df.to_csv(tmp_path, index=False, header=EXPORT_HEADER, encoding="utf-8", line_terminator="\n")
    except TypeError:
        out_df.to_csv(tmp_path, index=False, header=EXPORT_HEADER, encoding="utf-8", lineterminator="\n")
    with open(out_path, "w", encoding="utf-8") as f_out, open(tmp_path, "r", encoding="utf-8") as f_in:
        f_out.write(f"{HEADER_COMMENT}\n")
        f_out.writelines(f_in.readlines())
    tmp_path.unlink(missing_ok=True)
    print(f"File generated successfully: {out_path}")

if __name__ == "__main__":
    main()
PY

echo "Reaction network conversion completed successfully."

# -------------------------------
# Write copylist.pcp file
# -------------------------------
COPYLIST_FILE="${BASE_DIR}/copylist.pcp"

echo "[*] Creating copylist.pcp ..."
cat > "$COPYLIST_FILE" <<EOF
test.f90
profile.dat
solar_flux.txt
32H2SO4_photorate.txt
32SO2_OH_rate.txt
33COS_O_rate.txt
33COS_OH_rate.txt
33H2SO4_photorate.txt
33SO2_OH_rate.txt
34COS_O_rate.txt
34COS_OH_rate.txt
34H2SO4_photorate.txt
34SO2_OH_rate.txt
36COS_O_rate.txt
36COS_OH_rate.txt
36H2SO4_photorate.txt
36SO2_OH_rate.txt
H2SO4_wetdep.txt
partial_32H2SO4.txt
partial_33H2SO4.txt
partial_34H2SO4.txt
partial_36H2SO4.txt
Rainout-CS2.txt
Rainout-H2S.txt
Rainout-H2SO4.txt
Rainout-HSO3.txt
Rainout-OCS.txt
Rainout-SO2.txt
Rainout-SO4.txt
SO4_deposition_rate.txt
vapor_32H2SO4.txt
vapor_33H2SO4.txt
vapor_34H2SO4.txt
vapor_36H2SO4.txt
EOF
echo "copylist.pcp created at: $COPYLIST_FILE"

# -------------------------------
# Write test.f90 file
# -------------------------------
#TEST_F90_FILE="${BASE_DIR}/test.f90"
#
#echo "[*] Creating test.f90 ..."
#cat > "$TEST_F90_FILE" <<'EOF'
#program test
#  use patmo
#  use patmo_commons
#  use patmo_constants
#  use patmo_parameters
#  implicit none
#  real*8::dt,x(speciesNumber),t,tend,imass,one_year
#  real*8::heff(chemSpeciesNumber)
#  real*8::dep(chemSpeciesNumber)
#  integer::icell,i,j
#  real*8::convergence = 100.0
#  character (len = *), parameter :: FILE_NAME = "budget.nc"
#  integer :: ncid
#  integer, parameter :: NDIMS = 2
#  integer :: bud_varid
#  integer :: alt_dimid, bud_dimid, spec_dimid
#
#  call patmo_init()
#  call patmo_loadInitialProfile("profile.dat",unitH="km",unitX="1/cm3")
#  call patmo_setFluxBB()
#  call patmo_setGravity(9.8d2)
#  wetdep(:,:) = 0d0 
#
#  imass = patmo_getTotalMass()
#  print*,"mass:",imass
# 
#  dt = secondsPerDay 
#  tend = secondsPerDay*365*20d0 
#  one_year = secondsPerDay*365
#  t = 0d0
#
#  do
#     call patmo_run(dt,convergence)
#     t = t + dt        
#     if (t==tend .or. abs(convergence) < 1e-10) then
#        call patmo_run(dt,convergence)
#        call patmo_dumpDensityToFile(29,t,patmo_idx_O2)
#     endif  
#     print '(F11.2,a2)',t/tend*1d2," %"
#     if(t>=tend) exit
#  end do
#  
#  imass = patmo_getTotalMass()
#  print *,"mass:",imass
#  call patmo_dumpHydrostaticProfile("hydrostatEnd.out")
#  call patmo_dumpJValue("jvalue.dat")
#end program test
#EOF
#echo "test.f90 created at: $TEST_F90_FILE"

# -------------------------------
# Convert settings.xlsx → options.opt
# -------------------------------
SETTINGS_XLSX="${BASE_DIR}/settings.xlsx"
OPTIONS_OPT="${BASE_DIR}/options.opt"

echo "[*] Converting settings.xlsx to options.opt ..."
python3 - "$SETTINGS_XLSX" "$OPTIONS_OPT" "$BASE_DIR" <<'PY'
import sys
from pathlib import Path
import pandas as pd
import os

if len(sys.argv) < 4:
    raise SystemExit("Usage: python3 - <settings.xlsx> <options.opt> <BASE_DIR>")

in_path = Path(sys.argv[1]).expanduser().resolve()
out_path = Path(sys.argv[2]).expanduser().resolve()
base_dir = Path(sys.argv[3]).expanduser().resolve()

if not in_path.exists():
    raise SystemExit(f"settings.xlsx not found: {in_path}")

df = pd.read_excel(in_path, dtype="string", engine="openpyxl")
expected = ["Parameter", "unit", "input"]
if not all(col in df.columns for col in expected):
    raise SystemExit("settings.xlsx must have columns: Parameter, unit, input")
df = df.fillna("")
relative_path = os.path.relpath(base_dir / "reaction_network.ntw", Path.cwd()).replace("./", "")
with open(out_path, "w", encoding="utf-8") as f:
    f.write(f"network = {relative_path}\n")
    for _, row in df.iterrows():
        param = str(row["Parameter"]).strip()
        value = str(row["input"]).strip() if isinstance(row["input"], str) else ""
        f.write(f"{param} = {value}\n")
print(f"options.opt created at: {out_path}")
PY
echo "options.opt created successfully."

# -------------------------------
# Convert profile.xlsx → profile.dat
# -------------------------------
PROFILE_XLSX="${BASE_DIR}/profile.xlsx"
PROFILE_DAT="${BASE_DIR}/profile.dat"

echo "[*] Converting profile.xlsx to profile.dat ..."
python3 - "$PROFILE_XLSX" "$PROFILE_DAT" <<'PY'
import sys
from pathlib import Path
import pandas as pd
import numpy as np

if len(sys.argv) < 3:
    raise SystemExit("Usage: python3 - <profile.xlsx> <profile.dat>")

in_path = Path(sys.argv[1]).expanduser().resolve()
out_path = Path(sys.argv[2]).expanduser().resolve()

if not in_path.exists():
    raise SystemExit(f"profile.xlsx not found: {in_path}")

df = pd.read_excel(in_path, engine="openpyxl")
df = df.where(pd.notnull(df), "")

cols_keep = [0, 1]
for i, col in enumerate(df.columns):
    if i not in cols_keep:
        df[col] = df[col].apply(lambda x: f"{float(x):.4E}" if isinstance(x,(int,float,np.integer,np.floating)) else x)

total_columns = len(df.columns)
env_params = 5
species_count = total_columns - env_params

tmp_path = out_path.with_suffix(".tmp")
try:
    df.to_csv(tmp_path, sep="\t", index=False, header=True, encoding="utf-8", line_terminator="\n")
except TypeError:
    df.to_csv(tmp_path, sep="\t", index=False, header=True, encoding="utf-8", lineterminator="\n")

with open(out_path, "w", encoding="utf-8") as fout:
    fout.write(f"{env_params}\t{species_count}\n")
    with open(tmp_path, "r", encoding="utf-8") as fin:
        fout.writelines(fin.readlines())
tmp_path.unlink(missing_ok=True)
print(f"profile.dat created with header line: {env_params} {species_count}")
PY
echo "profile.dat created successfully."

# -------------------------------
# Build solar_flux.txt from solar_flux.xlsx + options.opt
# -------------------------------
SOLAR_XLSX="${BASE_DIR}/solar_flux.xlsx"
SOLAR_TXT="${BASE_DIR}/solar_flux.txt"
OPTIONS_OPT="${BASE_DIR}/options.opt"

echo "[*] Building solar_flux.txt in increasing energy order (decreasing wavelength) ..."

python3 - "$SOLAR_XLSX" "$SOLAR_TXT" "$OPTIONS_OPT" <<'PY'
import sys
from pathlib import Path
import pandas as pd
import numpy as np

if len(sys.argv) < 4:
    raise SystemExit("Usage: python3 - <solar_flux.xlsx> <solar_flux.txt> <options.opt>")

in_xlsx = Path(sys.argv[1]).expanduser().resolve()
out_txt = Path(sys.argv[2]).expanduser().resolve()
opt_path = Path(sys.argv[3]).expanduser().resolve()

if not in_xlsx.exists():
    raise SystemExit(f"solar_flux.xlsx not found: {in_xlsx}")
if not opt_path.exists():
    raise SystemExit(f"options.opt not found: {opt_path}")

params = {}
with open(opt_path,"r",encoding="utf-8") as f:
    for line in f:
        line=line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" in line:
            k,v=line.split("=",1)
            params[k.strip()] = v.strip()

def get_num(key,alt_keys=(),cast=float,required_name=None):
    keys_to_try = (key,)+tuple(alt_keys)
    for kk in keys_to_try:
        if kk in params and params[kk] != "":
            try:
                return cast(params[kk])
            except Exception:
                raise SystemExit(f"Invalid value for '{kk}': {params[kk]}")
    readable = required_name or key
    alts = ", ".join(alt_keys) if alt_keys else ""
    hint = f" (or {alts})" if alts else ""
    raise SystemExit(f"Missing '{readable}'{hint} in options.opt")

lam_min = get_num("wavelengMin",cast=float)
lam_max = get_num("wavelengMax",cast=float,required_name="wavelengMax")
NBIN = get_num("photoBinsNumber",cast=int)

if NBIN <= 1:
    raise SystemExit("photoBinsNumber must be > 1")
if not np.all(np.isfinite([lam_min, lam_max])) or min(lam_min, lam_max) <= 0:
    raise SystemExit("Wavelength limits must be finite and positive (nm)")
if lam_min > lam_max:
    lam_min, lam_max = lam_max, lam_min
if lam_min == lam_max:
    raise SystemExit("Wavelength limits must span a nonzero interval")

df = pd.read_excel(in_xlsx,engine="openpyxl")
cols = {c.strip():c for c in df.columns if isinstance(c,str)}
wcol = next((cols[c] for c in cols if c.lower().startswith("wavelength")),None)
icol = next((cols[c] for c in cols if c.lower().startswith("irradiance")),None)
if wcol is None or icol is None:
    raise SystemExit("Expected two columns: 'Wavelength (nm)' and 'Irradiance photon/(cm2 s nm)'")

df = df[[wcol,icol]].rename(columns={wcol:"lambda_nm",icol:"Iph"})
df = df.dropna(how="all")
df["lambda_nm"] = pd.to_numeric(df["lambda_nm"],errors="coerce")
df["Iph"] = pd.to_numeric(df["Iph"],errors="coerce")
invalid = ~np.isfinite(df[["lambda_nm", "Iph"]].to_numpy(dtype=float)).all(axis=1)
if invalid.any():
    rows = ", ".join(str(i + 2) for i in df.index[invalid][:5])
    raise SystemExit(f"Non-numeric, missing, or non-finite solar data in Excel row(s): {rows}")
if len(df) < 2:
    raise SystemExit("Solar spectrum must contain at least two valid data points")
df = df.sort_values("lambda_nm")

lam_src = df["lambda_nm"].to_numpy(dtype=float)
I_src = df["Iph"].to_numpy(dtype=float)
if np.any(lam_src <= 0) or np.any(I_src < 0):
    raise SystemExit("Solar wavelengths must be positive and photon fluxes nonnegative")
if np.any(np.diff(lam_src) <= 0):
    raise SystemExit("Solar wavelengths must be unique; resolve duplicate rows in Excel")
if lam_src[0] > lam_min or lam_src[-1] < lam_max:
    raise SystemExit(
        f"Solar spectrum [{lam_src[0]:g}, {lam_src[-1]:g}] nm does not cover "
        f"the requested [{lam_min:g}, {lam_max:g}] nm range"
    )

# Match src_py/patmo_photochemistry.py:createLogMetric, not a linear wavelength grid.
# Construct the grid here because build/xsecs/photoMetric.dat is generated later.
hc_eV_nm = 4.135667662e-15 * 2.99792458e10 * 1e7
energy_edges = 10.0 ** np.linspace(
    np.log10(hc_eV_nm / lam_max),
    np.log10(hc_eV_nm / lam_min),
    NBIN + 1,
)
lambda_long = hc_eV_nm / energy_edges[:-1]
lambda_short = hc_eV_nm / energy_edges[1:]

# Preserve the BASE endpoint-average convention and photon/(cm2 s nm) units.
# Coverage is checked above; np.interp clamps only endpoint roundoff here.
I_out = 0.5 * (
    np.interp(lambda_long, lam_src, I_src)
    + np.interp(lambda_short, lam_src, I_src)
)
# Rows already follow increasing energy (decreasing wavelength): do not reverse.
# The fixed integration width is applied by Fortran, not in this spectrum file.

np.savetxt(out_txt,I_out,fmt="%.8e")
print(f"solar_flux.txt created at: {out_txt} ({NBIN} bins, {lam_max:g} -> {lam_min:g} nm)")
PY

echo "solar_flux.txt created successfully."

# -------------------------------
# Generate PATMO files for the case containing this script
# -------------------------------
echo "[*] Running PATMO test with: -test=\"${USER_SUBDIR}\" ..."
python3 - "$USER_SUBDIR" <<'PY'
from collections import Counter
import csv
import math
from pathlib import Path
import re
import runpy
import shutil
import sys

case_name = sys.argv[1]
case_dir = Path("tests") / case_name
template_dir = case_dir / "templates"
for name in ("patmo_ode.f90", "patmo_parameters.f90"):
    if not (template_dir / name).is_file():
        raise SystemExit(f"Missing Yuhi template: {template_dir / name}")

# A rate written as read:<filename> denotes one scalar per atmospheric layer.
# Validate these inputs before PATMO starts replacing files in build/.
options_text = (case_dir / "options.opt").read_text()
cell_match = re.search(r"^cellsNumber\s*=\s*(\d+)\s*$", options_text, re.M)
if not cell_match or int(cell_match[1]) < 1:
    raise SystemExit("options.opt must specify a positive cellsNumber")
cells_number = int(cell_match[1])
external_files = []
rate_column = None
with (case_dir / "reaction_network.ntw").open(newline="") as handle:
    for row in csv.reader(handle):
        if not row or not row[0].strip() or row[0].startswith("#"):
            continue
        if row[0].startswith("@format:"):
            rate_column = [part.strip().lower() for part in row].index("rate")
            continue
        if row[0].startswith("@"):
            continue
        if rate_column is None or len(row) <= rate_column:
            raise SystemExit("Invalid reaction_network.ntw format")
        rate = row[rate_column].strip()
        if rate.lower() == "read":
            raise SystemExit(f"Reaction {row[0]}: replace bare read with read:<filename> in Excel")
        if rate.lower().startswith("read:"):
            name = rate[5:]
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name):
                raise SystemExit(f"Invalid external-rate filename: {name!r}")
            external_files.append(name)

for name in dict.fromkeys(external_files):
    path = case_dir / name
    if not path.is_file():
        raise SystemExit(f"External-rate file not found: {path}")
    rows = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if len(rows) != cells_number:
        raise SystemExit(f"{path}: expected {cells_number} rate values, found {len(rows)}")
    for layer, value in enumerate(rows, 1):
        try:
            number = float(value.replace("D", "e").replace("d", "e"))
        except ValueError:
            raise SystemExit(f"{path}: layer {layer} must contain one numeric rate")
        if not math.isfinite(number) or number < 0:
            raise SystemExit(f"{path}: layer {layer} must contain a finite, nonnegative rate")

# Keep newly referenced rate files in the generated copy list as well.
copylist_path = case_dir / "copylist.pcp"
copylist = copylist_path.read_text().splitlines()
for name in dict.fromkeys(external_files):
    if name not in copylist:
        copylist.append(name)
copylist_path.write_text("\n".join(copylist) + "\n")

# Use case-local templates without replacing the shared generator sources.
sys.path.insert(0, str(Path("src_py").resolve()))
import patmo_string

default_file_replace = patmo_string.fileReplaceBuild
def build_case_template(source_file, destination_file, *args, **kwargs):
    source = Path(source_file)
    override = template_dir / source.name
    if source.parent == Path("src_f90") and override.is_file():
        source_file = str(override)
    return default_file_replace(source_file, destination_file, *args, **kwargs)

patmo_string.fileReplaceBuild = build_case_template
Path("build").mkdir(exist_ok=True)
sys.argv = ["patmo", f"-test={case_name}"]
try:
    runpy.run_path("patmo", run_name="__main__")
except SystemExit as exc:
    # Some generator errors call sys.exit() without a nonzero status.
    if exc.code is None or exc.code == 0:
        raise SystemExit("PATMO stopped before generation completed; see the error above.")
    raise

shutil.copyfile(template_dir / "patmo_parameters.f90", "build/patmo_parameters.f90")

# PATMO supplies the current reaction indices. Convert its file-reference
# assignments to readers, so reordering the Excel rows cannot misassign rates.
rates_path = Path("build/patmo_rates.f90")
rates_source = rates_path.read_text()
assignment = re.compile(
    r"^[ \t]*krate\(icell,\s*(\d+)\)\s*=\s*read:([^\s]+)[ \t]*$",
    re.M | re.I,
)
matches = list(assignment.finditer(rates_source))
if Counter(match[2] for match in matches) != Counter(external_files):
    raise SystemExit("Generated external-rate references do not match reaction_network.ntw")

if matches:
    opens, closes = [], []
    replacements = {}
    for slot, match in enumerate(matches, 1):
        reaction_index, name = match[1], match[2]
        unit = f"rateFileUnits({slot})"
        opens.extend([
            f"    open(newunit={unit},file='{name}',status='old',action='read',iostat=rateFileStatus)",
            "    if (rateFileStatus /= 0) then",
            f"      write(*,*) 'Cannot open external-rate file: {name}'",
            "      error stop 1",
            "    end if",
        ])
        replacements[match.start()] = "\n".join([
            f"      read({unit},*,iostat=rateFileStatus) krate(icell,{reaction_index})",
            "      if (rateFileStatus /= 0) then",
            f"        write(*,*) 'Cannot read {name} at layer ', icell",
            "        error stop 1",
            "      end if",
        ])
        closes.append(f"    close({unit})")
    rates_source = assignment.sub(lambda match: replacements[match.start()], rates_source)

    def insert_once(pattern, addition, before=False):
        global rates_source
        anchor = re.compile(pattern, re.M | re.I)
        if len(list(anchor.finditer(rates_source))) != 1:
            raise SystemExit(f"Cannot locate external-rate insertion point: {pattern}")
        rates_source = anchor.sub(
            lambda match: addition + "\n" + match[0] if before else match[0] + "\n" + addition,
            rates_source,
        )

    insert_once(r"^[ \t]*implicit\s+none[ \t]*$",
                f"    integer :: rateFileUnits({len(matches)}), rateFileStatus")
    insert_once(r"^[ \t]*do\s+icell\s*=\s*1\s*,\s*cellsNumber[ \t]*$",
                "\n".join(opens), before=True)
    insert_once(r"^[ \t]*end\s+subroutine\s+computeRates[ \t]*$",
                "\n".join(closes), before=True)
    rates_path.write_text(rates_source)
print(f"[*] Generated external-rate readers for {len(matches)} reactions.")
PY

echo "[*] Generated files are in: ${PROJECT_ROOT}/build"
