# QDWithClustersMinimalV1 RK Tableaux

This repository contains **exported explicit Runge-Kutta tableaux** from the research project:

`Low Stage High Order Explicit Runge--Kutta Methods via Q- and D-Conditions: General Theory and Efficient Recursive Construction`

It is a data-and-usage repository: the paper is still in preparation, but the exported tableaux, file formats, and example scripts are already usable.

Use this repository if you want **precomputed low-stage, high-order, and very-high-order explicit RK methods** for **smooth, nonstiff ODE initial value problems**, and you mainly need the final Butcher tableau data `A`, `b`, `c` plus a minimal Julia workflow to load and test them.

This repository is **not a full solver package**. The included scripts focus on **fixed-step usage and convergence checks**. If your goal is stiffness handling, dense output, production solver features, or reconstructing the methods from scratch, this repository by itself is not enough.

Adaptive control is part of the broader research direction and will be discussed in the paper, but this repository does not currently provide standalone adaptive-step example code.

The **standard exported methods currently included here** have the following order-stage pairs:

- order 4 with 4 stages
- order 6 with 8 stages
- order 8 with 14 stages
- order 10 with 22 stages
- order 12 with 32 stages
- order 14 with 44 stages
- order 16 with 58 stages
- order 18 with 74 stages
- order 20 with 92 stages

**Typical use cases:**

- extremely high-accuracy time integration of smooth, nonstiff ODEs
- time discretization matched to high-order spatial schemes
- asymptotic convergence and error-regime studies
- accuracy-cost studies for very-high-order explicit methods
- direct reuse in fixed-step research solvers
- comparison between standard and optimized exported methods

## What Problem This Repository Solves

The original RK construction workflow is separate from this repository. Here the focus is much narrower: most users only need the final Butcher tableau and a minimal way to run it:

- stage matrix `A`
- weight vector `b`
- node vector `c`

This repository packages those results in reusable files and includes standalone Julia scripts so you can:

- load `csv` or `jld2`
- apply one tableau to a simple test problem
- run a convergence study with dyadic step sizes `h = 2^{-k}`

## Quick Start

Run all commands from the repository root.

### 1. Test a standard exported method from CSV

```bash
julia scripts/use_exported_rk_tableau.jl --input="standard methods/csv/QDWithClustersMinimalV1_p8_s14_bigfloat512.csv" --precision=BigFloat --prec=512 --k=4 --tfinal=1.0
```

### 2. Test a standard exported method from JLD2

```bash
julia scripts/use_exported_rk_tableau.jl --input="standard methods/jld2/QDWithClustersMinimalV1_p8_s14_bigfloat512.jld2" --precision=BigFloat --prec=512 --k=4 --tfinal=1.0
```

### 3. Run a convergence test

```bash
julia scripts/convergence_test_exported_rk_tableau.jl --input="standard methods/csv/QDWithClustersMinimalV1_p8_s14_bigfloat512.csv" --precision=BigFloat --prec=512 --k-range=2:8 --tfinal=1.0
```

### 4. Use the optimized method files

Optimized JLD2 file stored with `DoubleFloats.DoubleFloat{Float64}` coefficients:

```bash
julia scripts/use_exported_rk_tableau.jl --input="optimized methods/jld2/optimized_rk_p8_s14_double64_trial2.jld2" --k=4 --tfinal=1.0
```

The corresponding CSV file is also available:

```bash
julia scripts/use_exported_rk_tableau.jl --input="optimized methods/csv/optimized_rk_p8_s14_double64_trial2.csv" --precision=Double64 --k=4 --tfinal=1.0
```

## Repository Layout

```text
standard methods/
  csv/
  jld2/

optimized methods/
  csv/
  jld2/

scripts/
  use_exported_rk_tableau.jl
  convergence_test_exported_rk_tableau.jl
```

## Which File Should You Use?

Use `standard methods/...` if you want the standard exported tableaux from this repository.

Use `optimized methods/...` if you specifically want the optimized method files included here.

Use `csv` if:

- you want a human-readable format
- you want easier interoperability with non-Julia code
- you want explicit control over how coefficients are parsed, such as `Float64`, `BigFloat`, or `Double64`

Use `jld2` if:

- you want to preserve the Julia-side stored numeric types
- you want metadata alongside `A`, `b`, and `c`
- you are working entirely in Julia

## Fast Usage Guide

### Single Run Script

The main script is [`scripts/use_exported_rk_tableau.jl`](scripts/use_exported_rk_tableau.jl).

It:

- loads one exported tableau
- solves the test problem `u' = -u`, `u(0) = 1`
- compares the numerical result against `u(t) = e^{-t}`
- prints the absolute error

Example:

```bash
julia scripts/use_exported_rk_tableau.jl --input="standard methods/csv/QDWithClustersMinimalV1_p4_s4_bigfloat256.csv" --precision=BigFloat --prec=256 --k=4 --tfinal=1.0
```

### Convergence Test Script

The convergence script is [`scripts/convergence_test_exported_rk_tableau.jl`](scripts/convergence_test_exported_rk_tableau.jl).

It:

- runs the same test problem for a range of `k`
- uses `h = 2^{-k}`
- prints absolute errors and estimated convergence rates

Example:

```bash
julia scripts/convergence_test_exported_rk_tableau.jl --input="standard methods/jld2/QDWithClustersMinimalV1_p8_s14_bigfloat512.jld2" --precision=BigFloat --prec=512 --k-range=2:8 --tfinal=1.0
```

## Command-Line Arguments

Both scripts use `--key=value` syntax.

### Shared Arguments

- `--input=...`
  Path to a tableau file. Must point to a `.csv` or `.jld2` file.
- `--tfinal=...`
  Final time. Default is `1.0`.
- `--u0=...`
  Initial value. Default is `1.0`.

### Step Size Arguments

For `scripts/use_exported_rk_tableau.jl`:

- `--k=4`
  Uses `h = 2^{-4} = 1/16`.
- `--h=...`
  Directly specify the step size. You may also write `--h=2^-4`.

If neither `--k` nor `--h` is given, the script defaults to `k = 4`.

For `scripts/convergence_test_exported_rk_tableau.jl`:

- `--k-range=2:8`
  Runs `k = 2, 3, ..., 8`.

If omitted, the default is `2:8`.

### Precision Arguments

- `--precision=Float64`
- `--precision=BigFloat`
- `--precision=Double64`
- `--prec=256`
- `--prec=512`

Notes:

- `Float64` and `BigFloat` are Julia built-in numeric types. They do not require separate package installation.
- `Double64` is provided by `DoubleFloats`. If needed, the scripts will automatically install and load `DoubleFloats`.
- `--prec=...` is only relevant when using `BigFloat`.
- for `csv` files, `--precision` controls how coefficients are parsed
- for `jld2` files, omitting `--precision` means "use the coefficient type stored in the file"
- for `jld2` files, if you explicitly pass `--precision=Float64`, `--precision=BigFloat`, or `--precision=Double64`, the script converts the loaded tableau to that type before solving

## Automatic Package Handling

The scripts automatically install and load missing packages when needed.

### `JLD2`

If you use a `.jld2` input and `JLD2` is missing, the script will automatically install it.

### `DoubleFloats`

Some `jld2` files, such as:

```text
optimized methods/jld2/optimized_rk_p8_s14_double64_trial2.jld2
```

store coefficients as `DoubleFloats.DoubleFloat{Float64}`.

If such a file is used and `DoubleFloats` is missing, the script will automatically install and load `DoubleFloats`.

The same automatic handling also applies when you use `csv` files with `--precision=Double64`.

This applies to both:

- [`scripts/use_exported_rk_tableau.jl`](scripts/use_exported_rk_tableau.jl)
- [`scripts/convergence_test_exported_rk_tableau.jl`](scripts/convergence_test_exported_rk_tableau.jl)

## File Naming

Typical filenames look like:

```text
QDWithClustersMinimalV1_p8_s14_bigfloat512.csv
QDWithClustersMinimalV1_p8_s14_bigfloat512.jld2
```

Meaning:

- `p8`: method order is 8
- `s14`: number of stages is 14
- `bigfloat512`: coefficients were exported with `BigFloat` at 512-bit precision

## What Is Stored

Each method is stored as an explicit Runge-Kutta Butcher tableau:

```text
  c | A
  -------
    | b^T
```

where:

- `A` is the stage matrix of size `s x s`
- `b` is the weight vector of length `s`
- `c` is the node vector of length `s`

## CSV Format

CSV files are stored row by row as:

```text
c, A[1], A[2], ..., A[s], b
```

So a method with `s` stages has:

- `s` data rows
- `s + 2` columns

Example:

```csv
c, A[1], A[2], A[3], A[4], b
0.0,0.0,0.0,0.0,0.0,0.16666666666666666
0.5,0.5,0.0,0.0,0.0,0.3333333333333333
0.5,0.0,0.5,0.0,0.0,0.3333333333333333
1.0,0.0,0.0,1.0,0.0,0.16666666666666666
```

Reconstruction rule:

- first column is `c`
- middle `s` columns are `A`
- last column is `b`

## JLD2 Format

Each `*.jld2` file stores:

- `A`
- `b`
- `c`
- `metadata`
- `description`

The `metadata` dictionary typically contains:

- `order`
- `stages`
- `method`
- `precision`
- `date`

## Minimal Julia Usage

The examples below are intentionally minimal and show the basic file formats only. The provided scripts in [`scripts/use_exported_rk_tableau.jl`](scripts/use_exported_rk_tableau.jl) and [`scripts/convergence_test_exported_rk_tableau.jl`](scripts/convergence_test_exported_rk_tableau.jl) contain more robust type handling, conversion logic, and automatic package loading.

### Load from CSV

```julia
using DelimitedFiles

function load_rk_tableau_csv(path::String, T::Type=BigFloat)
    lines = readlines(path)
    rows = [split(line, ',') for line in lines[2:end] if !isempty(strip(line))]
    s = length(rows)

    c = Vector{T}(undef, s)
    A = Matrix{T}(undef, s, s)
    b = Vector{T}(undef, s)

    for i in 1:s
        c[i] = parse(T, strip(rows[i][1]))
        for j in 1:s
            A[i, j] = parse(T, strip(rows[i][j + 1]))
        end
        b[i] = parse(T, strip(rows[i][end]))
    end

    return A, b, c
end
```

### Load from JLD2

```julia
using JLD2

data = load("standard methods/jld2/QDWithClustersMinimalV1_p8_s14_bigfloat512.jld2")
A = data["A"]
b = data["b"]
c = data["c"]
metadata = data["metadata"]
```

## Using the Tableaux in Your Own Code

Once you have `A`, `b`, and `c`, you can plug them into any fixed-step explicit Runge-Kutta solver.

Example scalar stepper:

```julia
function rk_step(f, u, t, h, A, b, c)
    s = length(b)
    k = Vector{typeof(u)}(undef, s)

    for i in 1:s
        stage_state = u
        for j in 1:i-1
            stage_state += h * A[i, j] * k[j]
        end
        k[i] = f(stage_state, t + c[i] * h)
    end

    u_next = u
    for i in 1:s
        u_next += h * b[i] * k[i]
    end
    return u_next
end
```

For vector-valued ODEs, the same formula applies as long as `u`, `k[i]`, and `f(u, t)` support vector arithmetic.

## Precision Notes

- `jld2` preserves Julia-side stored coefficient types
- `csv` stores decimal strings and is easier to inspect manually
- if you need maximum fidelity within Julia, `jld2` is usually preferable
- if you need portability or explicit parse control, `csv` is usually preferable

For `BigFloat` workflows, it is recommended to set precision before parsing or solving:

```julia
setprecision(BigFloat, 512)
```

For convergence tests, dyadic step sizes are recommended:

```text
h = 2^{-k}
```

This reduces extra rounding effects compared with decimal values such as `0.1`.

## Scope

This repository provides exported explicit RK tableaux plus minimal fixed-step example scripts. It does not include:

- embedded error estimators
- dense output
- stiffness handling
- the original optimization / construction pipeline

It also does not yet include standalone adaptive-step example code, even though adaptive control will be part of the broader research framework described in the paper.

