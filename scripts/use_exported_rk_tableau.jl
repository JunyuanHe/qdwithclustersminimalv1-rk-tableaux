using DelimitedFiles
import Pkg

const JLD2_UUID = Base.UUID("033835bb-8acc-5ee8-8aae-3f567f8a3819")

function parse_kv_args(args)
    options = Dict{String, String}()
    for arg in args
        if startswith(arg, "--") && occursin("=", arg)
            key, value = split(arg[3:end], "=", limit = 2)
            options[key] = value
        end
    end
    return options
end

function parse_number(::Type{Float64}, s::AbstractString)
    return parse(Float64, s)
end

function parse_number(::Type{BigFloat}, s::AbstractString)
    return parse(BigFloat, s)
end

function pow2_step(::Type{T}, k::Int) where {T}
    k < 0 && error("k must be nonnegative in h = 2^{-k}")
    return ldexp(one(T), -k)
end

function parse_step_size(T::Type, options::Dict{String, String})
    if haskey(options, "k")
        return pow2_step(T, parse(Int, options["k"]))
    end

    raw_h = get(options, "h", "")
    if isempty(raw_h)
        return pow2_step(T, 4)
    end

    h_str = replace(strip(raw_h), ' ' => "")
    m = match(r"^2\^-(\d+)$", h_str)
    if m !== nothing
        return pow2_step(T, parse(Int, m.captures[1]))
    end

    return parse_number(T, raw_h)
end

function maybe_load_jld2()
    try
        Base.eval(Main, :(import JLD2))
        pkgid = Base.PkgId(JLD2_UUID, "JLD2")
        return Base.loaded_modules[pkgid]
    catch err
        if err isa ArgumentError && occursin("Package JLD2 not found", sprint(showerror, err))
            println("JLD2 is not installed in the current Julia environment. Installing it now...")
            try
                Pkg.add("JLD2")
                Base.eval(Main, :(import JLD2))
                pkgid = Base.PkgId(JLD2_UUID, "JLD2")
                return Base.loaded_modules[pkgid]
            catch install_err
                error(
                    "Automatic installation of JLD2 failed. " *
                    "Please make sure Julia can access the package registry and GitHub, " *
                    "then rerun the script. Original installation error: $(install_err)"
                )
            end
        end
        error("Failed to load JLD2 automatically. Original error: $(err)")
    end
end

function load_rk_tableau_csv(path::String, T::Type)
    lines = readlines(path)
    isempty(lines) && error("CSV file is empty: $path")
    length(lines) < 2 && error("CSV file has no data rows: $path")

    rows = [split(line, ',') for line in lines[2:end] if !isempty(strip(line))]
    s = length(rows)
    width = length(rows[1])
    width == s + 2 || error("Invalid RK CSV shape: expected $(s + 2) columns, got $width")
    any(length(row) != width for row in rows) && error("Inconsistent CSV row width in $path")

    c = Vector{T}(undef, s)
    A = Matrix{T}(undef, s, s)
    b = Vector{T}(undef, s)

    for i in 1:s
        row = rows[i]
        c[i] = parse_number(T, strip(row[1]))
        for j in 1:s
            A[i, j] = parse_number(T, strip(row[j + 1]))
        end
        b[i] = parse_number(T, strip(row[end]))
    end

    return (A = A, b = b, c = c, metadata = Dict{String, Any}())
end

function load_rk_tableau_jld2(path::String)
    JLD2_mod = maybe_load_jld2()
    data = Base.invokelatest(JLD2_mod.load, path)
    return (
        A = data["A"],
        b = data["b"],
        c = data["c"],
        metadata = get(data, "metadata", Dict{String, Any}()),
    )
end

function load_rk_tableau(path::String; T::Type = BigFloat)
    if endswith(lowercase(path), ".csv")
        return load_rk_tableau_csv(path, T)
    elseif endswith(lowercase(path), ".jld2")
        return load_rk_tableau_jld2(path)
    else
        error("Unsupported file extension for $path. Use .csv or .jld2")
    end
end

function rk_step(f, u, t, h, A, b, c)
    s = length(b)
    k = Vector{typeof(u)}(undef, s)

    for i in 1:s
        stage_state = u
        for j in 1:(i - 1)
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

function solve_rk(f, u0, tspan, h, A, b, c)
    t0, tf = tspan
    t = t0
    u = u0

    ts = typeof(t0)[t0]
    us = typeof(u0)[u0]

    while t < tf
        hstep = min(h, tf - t)
        u = rk_step(f, u, t, hstep, A, b, c)
        t += hstep
        push!(ts, t)
        push!(us, u)
    end

    return ts, us
end

function main(args)
    options = parse_kv_args(args)
    input = get(options, "input", "")
    isempty(input) && error("Missing --input=path/to/tableau.csv or .jld2")

    precision_name = get(options, "precision", "BigFloat")
    T = precision_name == "Float64" ? Float64 :
        precision_name == "BigFloat" ? BigFloat :
        error("Unsupported precision: $precision_name. Use Float64 or BigFloat.")

    if T == BigFloat
        precision_bits = parse(Int, get(options, "prec", "256"))
        setprecision(BigFloat, precision_bits)
    end

    h = parse_step_size(T, options)
    tfinal = parse_number(T, get(options, "tfinal", "1.0"))
    u0 = parse_number(T, get(options, "u0", "1.0"))

    tableau = load_rk_tableau(input; T = T)
    A, b, c = tableau.A, tableau.b, tableau.c

    f(u, t) = -u
    ts, us = solve_rk(f, u0, (zero(T), tfinal), h, A, b, c)

    exact = exp(-tfinal)
    err = abs(us[end] - exact)

    println("input     = ", abspath(input))
    println("stages    = ", length(b))
    if haskey(tableau.metadata, "order")
        println("order     = ", tableau.metadata["order"])
    end
    println("tfinal    = ", tfinal)
    println("h         = ", h)
    println("u_num     = ", us[end])
    println("u_exact   = ", exact)
    println("abs error = ", err)
    println()
    println("First few time points:")
    for i in 1:min(length(ts), 5)
        println("  t = ", ts[i], ", u = ", us[i])
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
