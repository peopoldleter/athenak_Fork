# AthenaK User Guide: Binary Neutron Star Simulations with Elliptica

This guide explains how to set up and run binary neutron star (BNS) merger
simulations in AthenaK using initial data provided by the
[Elliptica](https://github.com/rashti-alireza/Elliptica) initial data solver
and its companion `Elliptica_ID_Reader` library.

---

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Building AthenaK with Elliptica Support](#building-athenak-with-elliptica-support)
3. [Elliptica Initial Data File](#elliptica-initial-data-file)
4. [Configuring the Input File](#configuring-the-input-file)
5. [Running the Simulation](#running-the-simulation)
6. [Understanding the Outputs](#understanding-the-outputs)
7. [Restarting a Simulation](#restarting-a-simulation)
8. [Common Pitfalls](#common-pitfalls)
9. [Further Reading](#further-reading)

---

## Prerequisites

Before proceeding, make sure you have:

- **AthenaK source code** (this repository).
- **Kokkos** already installed or available as the bundled submodule
  (`kokkos/` inside the repository).
- **Elliptica_ID_Reader** library compiled and placed inside the repository
  root as `Elliptica_ID_Reader/`, with the following layout expected by CMake:
  ```
  Elliptica_ID_Reader/
  ├── include/
  │   └── elliptica_id_reader_lib.h
  └── lib/
      └── libelliptica_id_reader.a
  ```
- An Elliptica **initial data file** (the output of an Elliptica BNS solve).
- A C++17-capable compiler (GCC ≥ 8, Clang ≥ 7, or NVCC ≥ 11).
- *(Optional)* MPI for multi-node runs.
- *(Optional)* An EOS table in the format expected by Elliptica if you want a
  realistic equation of state instead of a piecewise-polytrope.

---

## Building AthenaK with Elliptica Support

The problem generator for Elliptica BNS data is `elliptica_bns`. You must
specify it at CMake configure time. The build system will automatically link
the `Elliptica_ID_Reader` library and enable OpenMP (required by the reader).

### Step 1 – Create a build directory

```bash
cd /path/to/athenak
mkdir build_elliptica && cd build_elliptica
```

### Step 2 – Configure with CMake

**CPU-only build** (serial or OpenMP):

```bash
cmake \
  -DPROBLEM=elliptica_bns \
  -DKokkos_ENABLE_OPENMP=ON \
  -DAthena_ENABLE_OPENMP=ON \
  -DAthena_ENABLE_MPI=ON \
  ..
```

**GPU build** (NVIDIA CUDA):

```bash
cmake \
  -DPROBLEM=elliptica_bns \
  -DKokkos_ENABLE_CUDA=ON \
  -DKokkos_ARCH_AMPERE86=ON \   # replace with your GPU architecture
  -DAthena_ENABLE_MPI=ON \
  ..
```

> **Note:** The Elliptica_ID_Reader library performs its interpolation on the
> CPU even when AthenaK runs on a GPU. The interpolated data is then copied to
> the GPU by AthenaK's problem generator. This is handled automatically.

### Step 3 – Build

```bash
make -j$(nproc)
```

This produces the executable `src/athena` inside the build directory.

---

## Elliptica Initial Data File

Elliptica writes its solution to a directory (referred to here as
`<elliptica_output_dir>`). The path to that directory is passed to AthenaK
via the `initial_data_file` parameter in the input file (see below).

The directory typically contains several field files and a parameter file. You
do not need to modify these files; AthenaK reads them through the
`Elliptica_ID_Reader` API.

If you are using a realistic (tabulated) equation of state, Elliptica requires
the path to the EOS table. AthenaK will forward this path to the reader when
you set `initial_data_table` in the `<problem>` block.

---

## Configuring the Input File

An example input file is provided at
[`inputs/dyngr/elliptica_bns.athinput`](../inputs/dyngr/elliptica_bns.athinput).

Below is a block-by-block explanation of all relevant parameters.

### `<comment>` – Human-readable description

```ini
<comment>
problem = Elliptica BNS initial data
```

### `<job>` – Output file prefix

```ini
<job>
basename = bns
```

All output files will be named `bns.<type>.<number>`.

### `<mesh>` – Global mesh and domain

```ini
<mesh>
nghost  = 4      # Number of ghost cells; must be ≥ 4 for WENOZ reconstruction
nx1     = 192    # Total zones in x-direction
x1min   = -1536  # Domain left edge (geometric units: G=c=M_sun=1)
x1max   =  1536  # Domain right edge
ix1_bc  = outflow
ox1_bc  = outflow
# Repeat for nx2/x2min/x2max/ix2_bc/ox2_bc and nx3/x3min/x3max/ix3_bc/ox3_bc
```

> **Tip:** The domain size should be much larger than the gravitational-wave
> extraction radii and the neutron star separation. A typical choice covers at
> least 1000 M in each direction for a standard BNS.

### `<meshblock>` – Per-block resolution

```ini
<meshblock>
nx1 = 32
nx2 = 32
nx3 = 32
```

MeshBlock size controls how work is distributed across MPI ranks and GPU
thread blocks. Sizes that are multiples of 32 work best on GPUs.

### `<mesh_refinement>` – Adaptive mesh refinement

```ini
<mesh_refinement>
refinement         = adaptive
max_nmb_per_rank   = 100   # Increase for large runs
num_levels         = 7     # Number of refinement levels
refinement_interval = 1    # Check refinement every N time steps
```

### `<amr_criterion0>` – Refinement criterion

```ini
<amr_criterion0>
method = user
```

The `user` method invokes the custom refinement condition defined in
`elliptica_bns.cpp`, which delegates to the Z4c AMR tracker.

### `<z4c_amr>` – Z4c AMR tracker settings

```ini
<z4c_amr>
method = tracker
```

This tracks the neutron stars and keeps the finest grid centred on them.

### `<refined_region1>` / `<refined_region2>` – Static initial refinement boxes

These blocks pre-refine the regions containing each neutron star before the
simulation begins, improving the quality of the initial data interpolation.

```ini
<refined_region1>
level = 6
x1min =   6   # Centred on star 1 (adjust to your binary separation)
x1max =  26
x2min = -10
x2max =  10
x3min = -10
x3max =  10

<refined_region2>
level = 6
x1min = -26   # Centred on star 2
x1max =  -6
x2min = -10
x2max =  10
x3min = -10
x3max =  10
```

### `<time>` – Time integration settings

```ini
<time>
evolution  = dynamic   # Must be "dynamic" for GRMHD
integrator = rk3       # Third-order Runge-Kutta
cfl_number = 0.25      # CFL < 0.3 is recommended for BNS
nlim       = -1        # No cycle limit; use tlim
tlim       = 5000      # End time in geometric units
ndiag      = 1
```

### `<coord>` – Coordinate and spacetime settings

```ini
<coord>
general_rel = true    # Enable GR
excise      = false   # No BH excision (appropriate for BNS)
m           = 0.0     # Irrelevant for BNS; set to 0
a           = 0.0
```

### `<mhd>` – MHD and equation of state

```ini
<mhd>
eos         = ideal          # Required syntax; actual EOS set by dyn_eos
gamma       = 2.0            # Used only as a placeholder

dyn_eos     = piecewise_poly # Piecewise polytrope (change to "ideal" for Gamma-law)
dyn_error   = reset_floor    # Error handling policy in primitive recovery
reconstruct = wenoz           # WENOZ reconstruction (recommended for BNS)
rsolver     = llf            # LLF Riemann solver (robust choice)

dfloor      = 1.28e-21       # Density floor (geometric units)
tfloor      = 1.58255e-19    # Temperature floor
gamma_max   = 20.025         # Maximum Lorentz factor

# Piecewise polytrope parameters (SLy EOS example):
pwp_poly_rmd          = 5.340089573220712e+23
pwp_density_pieces_0  = 0.000000000000000e+00
pwp_density_pieces_1  = 1.136499842827129e+17
pwp_density_pieces_2  = 5.011872336272714e+17
pwp_density_pieces_3  = 9.999999999999999e+17
pwp_gamma_pieces_0    = 1.356920000000000e+00
pwp_gamma_pieces_1    = 3.456000000000000e+00
pwp_gamma_pieces_2    = 3.011000000000000e+00
pwp_gamma_pieces_3    = 1.425000000000000e+00
pwp_gamma_thermal     = 1.7
```

> **Tabulated EOS:** If you generated Elliptica data with a tabulated EOS,
> set `dyn_eos = table` here and provide the table path via
> `initial_data_table` in the `<problem>` block.

### `<adm>` – ADM variables block (leave empty)

```ini
<adm>
```

### `<z4c>` – Z4c numerical relativity parameters

```ini
<z4c>
# Gauge
lapse_oplog     = 2.0    # 1+log slicing parameter
lapse_harmonicf = 1.0
lapse_harmonic  = 0.0
lapse_advect    = 1.0
shift_eta       = 0.3    # Gamma-driver damping; ~1/M_ADM
shift_advect    = 1.0

# Kreiss-Oliger dissipation
diss            = 0.5

# Constraint damping
chi_div_floor   = 1e-05
damp_kappa1     = 0.02
damp_kappa2     = 0.0

# Gravitational-wave extraction
nrad_wave_extraction = 3
extraction_radius_1  = 200
extraction_radius_2  = 400
extraction_radius_3  = 600
waveform_dt          = 2.5

# Neutron star tracker (positions in geometric units)
nco          = 2
co_0_type    = NS
co_0_x       = 16.0     # Initial x-position of star 1 (from Elliptica output)
co_0_radius  = 8.2
co_1_type    = NS
co_1_x       = -16.0    # Initial x-position of star 2
co_1_radius  = 8.2
```

> **Important:** The neutron star tracker positions (`co_0_x`, `co_1_x`) must
> match the positions reported by Elliptica in its output. Check the Elliptica
> log or the initial data parameter file for the coordinate separation.

### `<problem>` – Elliptica-specific parameters

```ini
<problem>
initial_data_file  = /path/to/elliptica_output_dir   # REQUIRED
initial_data_table = /path/to/eos_table.h5            # Optional; only for tabulated EOS
user_hist          = true
```

- `initial_data_file` is the path to the Elliptica output directory containing
  the field data files. This is the only mandatory parameter.
- `initial_data_table` is the path to the EOS table used during the Elliptica
  solve. Omit this line if you used a polytropic EOS.

### `<output*>` – Output configuration

```ini
<output1>
file_type = hst         # History: max density, min lapse (every 0.375 M)
dt        = 0.375

<output2>
file_type = bin         # 2-D equatorial slice of primitive MHD variables
variable  = mhd_w
dt        = 1.5
slice_x3  = 0.0

<output3>
file_type = bin
variable  = z4c         # Z4c spacetime variables
dt        = 1.5
slice_x3  = 0.0

<output4>
file_type = bin
variable  = adm         # ADM metric variables
dt        = 1.5
slice_x3  = 0.0

<output5>
file_type = rst         # Restart checkpoint (every 375 M)
dt        = 375.0
```

---

## Running the Simulation

On a university HPC cluster the recommended way to run AthenaK is to submit
the simulation as a **SLURM batch job** via `sbatch`. This lets the scheduler
allocate the requested resources, keeps the job alive after you log out, and
allows you to chain restart jobs automatically.

### Step 1 – Copy and adapt the job script

Two template job scripts are provided below. Copy the one that matches your
hardware, save it (e.g. as `run_bns.sh`), and adjust every line marked with
`# EDIT`.

---

#### Template A – CPU / MPI (multi-core or multi-node)

```bash
#!/bin/bash
#SBATCH --job-name=athenak_bns          # Name shown in the queue
#SBATCH --partition=standard            # EDIT: your cluster's partition/queue
#SBATCH --nodes=4                       # EDIT: number of nodes
#SBATCH --ntasks-per-node=32            # EDIT: MPI ranks per node (≤ cores/node)
#SBATCH --cpus-per-task=1               # OpenMP threads per rank (increase if needed)
#SBATCH --mem=0                         # Use all available memory on each node
#SBATCH --time=48:00:00                 # EDIT: wall-clock limit (HH:MM:SS)
#SBATCH --output=bns_%j.out             # stdout  (%j = job ID)
#SBATCH --error=bns_%j.err              # stderr

# --- Environment ---------------------------------------------------------
module purge
module load gcc/12 openmpi/4.1          # EDIT: your site's MPI module

# If AthenaK was built with OpenMP support, set the thread count here:
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PROC_BIND=close

# --- Paths ---------------------------------------------------------------
ATHENAK_BUILD=/path/to/build_elliptica  # EDIT
ATHINPUT=/path/to/inputs/dyngr/elliptica_bns.athinput  # EDIT
ELLIPTICA_DATA=/path/to/elliptica_output_dir           # EDIT

# --- Run -----------------------------------------------------------------
cd ${ATHENAK_BUILD}

srun ./src/athena \
  -i ${ATHINPUT} \
  problem/initial_data_file=${ELLIPTICA_DATA}
```

> **Tip – MeshBlock count ≥ MPI ranks:**
> Make sure `(nx1/mb_nx1) × (nx2/mb_nx2) × (nx3/mb_nx3)` is at least equal to
> `--nodes × --ntasks-per-node`, otherwise AthenaK will abort at startup.

---

#### Template B – GPU / CUDA (one MPI rank per GPU)

```bash
#!/bin/bash
#SBATCH --job-name=athenak_bns_gpu      # Name shown in the queue
#SBATCH --partition=gpu                 # EDIT: your GPU partition/queue
#SBATCH --nodes=2                       # EDIT: number of GPU nodes
#SBATCH --ntasks-per-node=4             # EDIT: GPUs per node (= MPI ranks per node)
#SBATCH --gpus-per-task=1               # One GPU per MPI rank
#SBATCH --cpus-per-task=8               # EDIT: CPU cores available to each rank
#SBATCH --mem=0
#SBATCH --time=48:00:00                 # EDIT: wall-clock limit
#SBATCH --output=bns_gpu_%j.out
#SBATCH --error=bns_gpu_%j.err

# --- Environment ---------------------------------------------------------
module purge
module load gcc/12 cuda/12.2 openmpi/4.1   # EDIT: your site's modules

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# --- Paths ---------------------------------------------------------------
ATHENAK_BUILD=/path/to/build_elliptica  # EDIT
ATHINPUT=/path/to/inputs/dyngr/elliptica_bns.athinput  # EDIT
ELLIPTICA_DATA=/path/to/elliptica_output_dir           # EDIT

# --- Run -----------------------------------------------------------------
cd ${ATHENAK_BUILD}

# srun automatically sets CUDA_VISIBLE_DEVICES for each rank when
# --gpus-per-task=1 is used; no manual binding needed in most cases.
srun ./src/athena \
  -i ${ATHINPUT} \
  problem/initial_data_file=${ELLIPTICA_DATA}
```

---

### Step 2 – Submit the job

```bash
sbatch run_bns.sh
```

Monitor the job with:

```bash
squeue -u $USER          # show your running/pending jobs
scontrol show job <JOBID> # detailed job information
```

---

### Step 3 – Restarting within SLURM (checkpoint chains)

When the wall-clock limit is reached SLURM will kill the job. Use AthenaK's
restart capability to continue automatically by submitting a restart job from
inside the job script:

```bash
#!/bin/bash
#SBATCH --job-name=athenak_bns_restart
#SBATCH --partition=standard            # EDIT
#SBATCH --nodes=4                       # EDIT (must match original run)
#SBATCH --ntasks-per-node=32            # EDIT
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --time=48:00:00
#SBATCH --output=bns_restart_%j.out
#SBATCH --error=bns_restart_%j.err

module purge
module load gcc/12 openmpi/4.1         # EDIT

ATHENAK_BUILD=/path/to/build_elliptica # EDIT
OUTPUT_DIR=/path/to/output             # EDIT: directory containing .rst files

cd ${ATHENAK_BUILD}

# Pick the latest restart file automatically
RST_FILE=$(ls -t ${OUTPUT_DIR}/bns.rst.?????.rst 2>/dev/null | head -1)

if [ -z "${RST_FILE}" ]; then
  echo "ERROR: no restart file found in ${OUTPUT_DIR}" >&2
  exit 1
fi

echo "Restarting from ${RST_FILE}"
srun ./src/athena -r ${RST_FILE}
```

You can also self-chain a restart by having SLURM submit the next job before
the current one ends:

```bash
# At the end of your original run_bns.sh, add:
sbatch run_bns_restart.sh
```

---

### Command-line parameter overrides

Any input-file parameter can be overridden on the command line (works with
both `sbatch` scripts and interactive runs):

```bash
srun ./src/athena -i elliptica_bns.athinput \
  time/tlim=10000 \
  problem/initial_data_file=/new/path/to/elliptica_data
```

---

### Quick interactive test (single node, no scheduler)

For short test runs on a login node or a single interactive allocation:

```bash
cd /path/to/build_elliptica
# Serial / OpenMP
./src/athena -i /path/to/inputs/dyngr/elliptica_bns.athinput

# MPI (use srun inside an salloc session, not mpirun, on most clusters)
srun -n <N_ranks> ./src/athena -i /path/to/inputs/dyngr/elliptica_bns.athinput
```

---

## Understanding the Outputs

### History file (`.hst`)

Written every `dt` time units. Contains two columns beyond the standard header:

| Column | Quantity |
|--------|----------|
| `rho-max` | Maximum rest-mass density across the grid |
| `alpha-min` | Minimum lapse function (drops toward 0 near a collapsing core) |

Use the maximum density to track stellar oscillations or collapse onset.

### Binary data files (`.bin`)

Written by the `bin` output type. Each file contains one or more variables for
a single output time. These can be read with the Python utilities in
`vis/python/`:

```python
import sys
sys.path.insert(0, '/path/to/athenak/vis/python')
import athena_read
data = athena_read.athdf('bns.mhd_w.00100.bin')
```

**Key variables** for BNS analysis:

| Variable name | Description |
|---|---|
| `mhd_w` | All primitive variables (density, velocity, pressure, B-field) |
| `mhd_w_d` | Rest-mass density ρ |
| `mhd_w_vx/vy/vz` | Fluid 3-velocity components |
| `mhd_w_p` | Gas pressure |
| `mhd_bcc1/2/3` | Cell-centred magnetic field components |
| `z4c` | Full set of Z4c evolution variables |
| `adm` | ADM lapse α, shift β^i, spatial metric g_ij, extrinsic curvature K_ij |

### Gravitational wave output

The Z4c module computes the Newman-Penrose scalar Ψ₄ at the extraction radii
defined by `extraction_radius_1/2/3`. The waveforms are written as additional
history columns. Use a standard GW post-processing toolkit (e.g., `kuibit` or
`GWpy`) to extract the strain from Ψ₄.

### Restart files (`.rst`)

Written by the `rst` output type. Restart from a checkpoint with:

```bash
./src/athena -r bns.rst.00010.rst
```

---

## Restarting a Simulation

To restart from the last checkpoint:

```bash
./src/athena -r bns.rst.<number>.rst
```

To change parameters at restart, append `-i <input_file>` and override
specific blocks on the command line. The restart file always takes precedence
for mesh and physics state; only output and time limit overrides are safe to
change.

---

## Common Pitfalls

| Problem | Likely cause | Fix |
|---------|-------------|-----|
| Fatal: `BNS data requires <mhd> and <z4c> blocks` | Missing `<mhd>` or `<z4c>` block in input file | Add both blocks |
| Code exits immediately with no output | `initial_data_file` path is wrong | Check the path and make sure the Elliptica output directory exists |
| `The velocity is superluminal!` warnings | Physical; the reader will rescale automatically | Monitor the count; if frequent, consider a stricter density floor |
| Simulation blows up early | CFL too large, or floors too small | Reduce `cfl_number` (try 0.2), increase `dfloor` |
| AMR tracker loses a star | Initial tracker positions wrong | Set `co_0_x` / `co_1_x` to match the Elliptica coordinate separation exactly |
| MPI hangs at startup | More ranks than MeshBlocks | Ensure `(nx1/mb_nx1) × (nx2/mb_nx2) × (nx3/mb_nx3) ≥ N_ranks` |
| OpenMP link error when building | Elliptica reader requires OpenMP | Add `-DAthena_ENABLE_OPENMP=ON` or CMake will add it automatically |

---

## Further Reading

- [AthenaK wiki](https://github.com/IAS-Astrophysics/athenak/wiki) – general
  build and run instructions.
- [Stone et al. (2024)](https://ui.adsabs.harvard.edu/abs/2024arXiv240916053S/abstract) –
  AthenaK framework paper.
- [Zhu et al. (2024)](https://ui.adsabs.harvard.edu/abs/2024arXiv240910383Z/abstract) –
  Z4c numerical relativity solver.
- [Fields et al. (2024)](https://ui.adsabs.harvard.edu/abs/2024arXiv240910384F/abstract) –
  GR hydro and MHD in dynamical spacetimes.
- [Elliptica](https://github.com/rashti-alireza/Elliptica) – initial data
  solver used to generate the BNS initial data.
- `src/pgen/elliptica_bns.cpp` – the AthenaK problem generator that reads
  Elliptica data; consult this file for implementation details.
- `inputs/dyngr/elliptica_bns.athinput` – the annotated example input file
  provided with this repository.
