# Exported RK Tableaux

This repository contains precomputed explicit Runge-Kutta tableaux for `QDWithClustersMinimalV1`.

The files are provided in two formats:

- `*.csv`: human-readable text format
- `*.jld2`: native Julia binary format with metadata

These files are intended for users who want to apply the methods directly, without access to the original construction code.

## File Naming

Typical filenames look like:

```text
QDWithClustersMinimalV1_p8_s14_bigfloat512.csv
QDWithClustersMinimalV1_p8_s14_bigfloat512.jld2
```

Meaning:

- `p8`: method order is 8
- `s14`: number of stages is 14
- `bigfloat512`: coefficients were stored as `BigFloat` with 512-bit precision

## What Is Stored

Each method is represented by the standard explicit Runge-Kutta Butcher tableau:

```text
  c | A
  -------
    | b^T
```

where:

- `A` is the stage matrix, of size `s x s`
- `b` is the weight vector, of length `s`
- `c` is the node vector, of length `s`

## CSV Format

The CSV files store the tableau row by row in the following layout:

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

To reconstruct the tableau:

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

### 1. Load from CSV

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

### 2. Load from JLD2

```julia
using JLD2

data = load("QDWithClustersMinimalV1_p8_s14_bigfloat512.jld2")
A = data["A"]
b = data["b"]
c = data["c"]
metadata = data["metadata"]
```

## Time Stepping Formula

Given an ODE

```math
u'(t) = f(u,t),
\qquad u(t_n) = u_n,
```

an explicit Runge-Kutta step with step size `h` is:

```math
k_i = f\left(u_n + h \sum_{j=1}^{i-1} a_{ij} k_j,\; t_n + c_i h\right), \qquad i=1,\dots,s
```

```math
u_{n+1} = u_n + h \sum_{i=1}^s b_i k_i
```

## Standalone Example Script

This repository includes a self-contained Julia example:

- [`scripts/use_exported_rk_tableau.jl`](scripts/use_exported_rk_tableau.jl)
- [`scripts/convergence_test_exported_rk_tableau.jl`](scripts/convergence_test_exported_rk_tableau.jl)

It does not depend on the original construction code. It can:

- load `csv`
- load `jld2`
- apply one of these RK methods to a simple test ODE
- run a small convergence study for `h = 2^{-k}`

Example:

```bash
julia scripts/use_exported_rk_tableau.jl --input="standard methods/csv/QDWithClustersMinimalV1_p8_s14_bigfloat512.csv" --precision=BigFloat --prec=512 --k=4 --tfinal=1.0
```

or

```bash
julia scripts/use_exported_rk_tableau.jl --input="standard methods/jld2/QDWithClustersMinimalV1_p8_s14_bigfloat512.jld2" --precision=BigFloat --prec=512 --k=4 --tfinal=1.0
```

Here `--k=4` means `h = 2^{-4} = 1/16`. The script also accepts `--h=2^-4`.

Notes:

- run the commands from the repository root
- `--input` must be written as `--input=...`
- because the folder name contains a space (`standard methods`), quote the path as shown above
- if you want to use `Float64`, choose a file whose coefficients were exported for that precision, for example:

```bash
julia scripts/use_exported_rk_tableau.jl --input="optimized methods/optimized_rk_p8_s14_double64_trial2.jld2" --precision=Float64 --k=4 --tfinal=1.0
```

The example solves:

```math
u' = -u, \qquad u(0)=1
```

and compares the numerical result with the exact solution `u(t)=e^{-t}`.

For a small convergence test over `k = 2:8`, run:

```bash
julia scripts/convergence_test_exported_rk_tableau.jl --input="standard methods/csv/QDWithClustersMinimalV1_p8_s14_bigfloat512.csv" --precision=BigFloat --prec=512 --k-range=2:8 --tfinal=1.0
```

It prints a table with:

- `k`
- `h = 2^{-k}`
- absolute error at `t = tfinal`
- estimated rate `log2(err_k / err_{k+1})`

## Using the Tableaux in Your Own Code

The only ingredients you need are `A`, `b`, and `c`. Once loaded, you can plug them into any fixed-step explicit Runge-Kutta solver.

For scalar ODEs:

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

For vector-valued ODEs, the same formula applies as long as `u`, `k[i]`, and `f(u,t)` support vector arithmetic.

## Precision Notes

- `jld2` preserves Julia numeric types directly, including `BigFloat`
- `csv` stores decimal strings and is easier to inspect manually
- if you need maximum fidelity in Julia, prefer `jld2`
- if you need portability to other languages, prefer `csv`

For `BigFloat` workflows, it is recommended to set precision before parsing:

```julia
setprecision(BigFloat, 512)
```

For convergence tests, it is recommended to use dyadic step sizes:

```text
h = 2^{-k}
```

This reduces extra rounding effects compared with decimal inputs such as `0.1`.

## Dependencies

For CSV usage:

- Julia standard library `DelimitedFiles`

For JLD2 usage:

- Julia package `JLD2`

The standalone script [`scripts/use_exported_rk_tableau.jl`](scripts/use_exported_rk_tableau.jl) will automatically install `JLD2` on first use if it is missing.

If you prefer to install it manually:

```julia
import Pkg
Pkg.add("JLD2")
```

## Scope

These files provide fixed explicit RK tableaux only. They do not include:

- adaptive step-size control
- embedded error estimators
- dense output
- stiffness handling

If you need those features, use these tableaux as the core stepping coefficients inside your own solver framework.
