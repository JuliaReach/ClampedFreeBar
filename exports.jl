include("clamped.jl")

using MAT

function write_mat(; N)
    sys = clamped_canonical(N=100)
    dat = Dict("A" => sys.A, "b" => sys.c)
    matwrite("CB21_$(N).mat", dat; compress = true)
end

write_mat(N=100)
write_mat(N=500)
write_mat(N=1000)
