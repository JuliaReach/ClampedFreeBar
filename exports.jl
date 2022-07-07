include("clamped.jl")

# ====== Converting to MAT files =======

using MAT

function write_mat(; N, damped)
    name = damped ? "CB21d_$N.mat" : "CB21_$N.mat"
    sys = clamped(N=N, a=1e-6, b=1e-6, damped=damped, homogeneize=false)
    dat = Dict("A" => sys.A, "b" => sys.c)
    return matwrite(name, dat; compress = true)
end

#=
for N in [100, 500, 1000]
    for damped in [true, false]
        write_mat(N=N, damped=damped)
    end
end
=#

# ====== Converting to SX files =======

using SpaceExParser

function write_sxmodel_clamped_constant(; N, damped)
    name = damped ? "CB21Cd_$N.xml" : "CB21C_$N.xml"
    sys = clamped(N=N, a=1e-6, b=1e-6, damped=damped, homogeneize=true)
    return writesxmodel(name, sys)
end

function write_sxmodel_clamped_timevarying(U = Interval(0.99, 1.01); N, damped)
    X = Universe(N)
    name = damped ? "CB21Fd_$N.xml" : "CB21F_$N.xml"
    sys =  clamped(N=N, a=1e-6, b=1e-6, damped=damped, homogeneize=false)    
    A = sys.A
    B = hcat(sys.c)
    S = @system(x' = A*x + B*u, x ∈ X, u ∈ U)
    return writesxmodel(name, S)
end

#=
for func in [write_sxmodel_clamped_constant, write_sxmodel_clamped_timevarying] 
for N in [100, 500, 1000]
    for damped in [true, false]
        func(N=N, damped=damped)
    end
end
end
=#
