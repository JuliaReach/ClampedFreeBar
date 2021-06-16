include("clamped.jl")

# ====== Converting to MAT files =======

using MAT

function write_mat(; N, damped=true)

    if damped
        sys = clamped(N=N, a=1e-6, b=1e-6, damped=true, homogeneize=false)    
        name = "CB21d_$N.mat"
       
    else
        sys = clamped(N=N, damped=false, homogeneize=false)    
        name = "CB21_$N.mat"
    end
    
    dat = Dict("A" => sys.A, "b" => sys.c)
    matwrite(name, dat; compress = true)
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

# WARNING requires manual modification of the xml
function write_sxmodel_clamped_constant(; N, damped=true)
    if damped
        sys = clamped(N=N, a=1e-6, b=1e-6, damped=true, homogeneize=true)
        name = "CB21Cd_$N.xml"
    else
        sys = clamped(N=N, damped=false, homogeneize=true)
        name = "CB21C_$N.xml"
    end
    writesxmodel(name, sys)
end

function write_sxmodel_clamped_timevarying(; N, U = Interval(9900, 10100), damped=true)
    
    X = Universe(N)

    if damped
        sys = clamped(N=N, a=1e-6, b=1e-6, damped=true, homogeneize=false)
        name = "CB21Fd_$N.xml"

    else
        sys = clamped(N=N, damped=false, homogeneize=false)
        name = "CB21F_$N.xml"
    end
    
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
