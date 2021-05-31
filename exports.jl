include("clamped.jl")

# ====== Converting to MAT files =======

using MAT

function write_mat(; N)
    sys = clamped_canonical(N=100)
    dat = Dict("A" => sys.A, "b" => sys.c)
    matwrite("CB21_$(N).mat", dat; compress = true)
end

write_mat(N=100)
write_mat(N=500)
write_mat(N=1000)

# ====== Converting to SX files =======

function write_sxmodel_clamped_constant(; N, U = Interval(9900, 10100))
    aux = clamped_canonical(N=N)

    Aext = [aux.A aux.c; zeros(1, 2N+1)]
    Sac = @system(x' = Aext * x);

    return writesxmodel("CB21C_$N.xml", Sac)
end

function write_sxmodel_clamped_timevarying(; N, U = Interval(9900, 10100))
    aux = clamped_canonical(N=N)

    X = Universe(N)
    A = aux.A
    B = hcat(aux.c)
    Sac = @system(x' = A*x + B*u, x ∈ X, u ∈ U);

    return writesxmodel("CB21F_$N.xml", Sac)
end

write_sxmodel_clamped_constant(N=100)
write_sxmodel_clamped_constant(N=500)
write_sxmodel_clamped_constant(N=1000)

write_sxmodel_clamped_timevarying(N=100)
write_sxmodel_clamped_timevarying(N=500)
write_sxmodel_clamped_timevarying(N=1000)