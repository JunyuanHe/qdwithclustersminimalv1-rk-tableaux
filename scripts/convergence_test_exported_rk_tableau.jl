using Printf

include(joinpath(@__DIR__, "use_exported_rk_tableau.jl"))

function parse_k_range(spec::String)
    text = replace(strip(spec), ' ' => "")
    m = match(r"^(\d+):(\d+)$", text)
    m === nothing && error("Invalid --k-range format: $spec. Expected start:stop, e.g. 2:8")
    k_start = parse(Int, m.captures[1])
    k_stop = parse(Int, m.captures[2])
    k_start <= k_stop || error("Invalid --k-range: start must be <= stop")
    return collect(k_start:k_stop)
end

function solve_test_problem(T::Type, h, A, b, c; tfinal, u0)
    f(u, t) = -u
    _, us = solve_rk(f, u0, (zero(T), tfinal), h, A, b, c)
    exact = exp(-tfinal)
    err = abs(us[end] - exact)
    return us[end], exact, err
end

function print_error_table(T::Type, ks, errors)
    println()
    println("k    h                    abs error            rate")
    println("--------------------------------------------------------")
    for (idx, k) in enumerate(ks)
        h = pow2_step(T, k)
        rate_str = if idx == 1
            "-"
        else
            @sprintf("%.6f", log2(errors[idx - 1] / errors[idx]))
        end
        error_str = @sprintf("%.8e", errors[idx])
        println(rpad(string(k), 5), rpad(string(h), 21), rpad(error_str, 21), rate_str)
    end
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

    ks = parse_k_range(get(options, "k-range", "2:8"))
    tfinal = parse_number(T, get(options, "tfinal", "1.0"))
    u0 = parse_number(T, get(options, "u0", "1.0"))

    tableau = load_rk_tableau(input; T = T)
    A, b, c = tableau.A, tableau.b, tableau.c

    values = Vector{T}(undef, length(ks))
    exacts = Vector{T}(undef, length(ks))
    errors = Vector{T}(undef, length(ks))

    for (idx, k) in enumerate(ks)
        h = pow2_step(T, k)
        u_num, u_exact, err = solve_test_problem(T, h, A, b, c; tfinal = tfinal, u0 = u0)
        values[idx] = u_num
        exacts[idx] = u_exact
        errors[idx] = err
    end

    println("input     = ", abspath(input))
    println("stages    = ", length(b))
    if haskey(tableau.metadata, "order")
        println("order     = ", tableau.metadata["order"])
    end
    println("tfinal    = ", tfinal)
    println("u0        = ", u0)
    println("k range   = ", first(ks), ":", last(ks))
    println("problem   = u' = -u, u(0) = ", u0)
    println("exact     = exp(-tfinal)")

    print_error_table(T, ks, errors)

    println()
    println("finest h  = ", pow2_step(T, last(ks)))
    println("u_num     = ", values[end])
    println("u_exact   = ", exacts[end])
    println("abs error = ", errors[end])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
