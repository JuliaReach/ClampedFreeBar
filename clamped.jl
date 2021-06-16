#=
This model is taken from [1]. The theoretical derivations of the analytic solution can be found in [2].

[1] Malakiyeh, Mohammad Mahdi, Saeed Shojaee, and Klaus-Jürgen Bathe. "The Bathe time integration method revisited for prescribing desired numerical dissipation." Computers & Structures 212 (2019): 289-298.

[2] Mechanical Vibrations, Gerardin et al, page 250-251.
=#

using ReachabilityAnalysis, LinearAlgebra, LazySets
using ReachabilityAnalysis: solve, discretize
using LazySets.Arrays
using SparseArrays

LazySets.set_ztol(Float64, 1e-15)
LazySets.set_atol(Float64, 1e-15)
LazySets.set_rtol(Float64, 1e-15)

function clamped_matrices(; N=1000, E=30e6, ρ=7.3e-4, A=1, L=200)
    ℓ = L / N
    K = (E * A / ℓ) * SymTridiagonal(fill(2.0, N), fill(-1.0, N))
    K[end, end] = E * A / ℓ

    M = (ρ * A * ℓ / 2) * Diagonal(vcat(fill(2.0, N-1), 1.0))
    M[end, end] = ρ * A * ℓ / 2

    return M, K
end

function clamped_free(; N)
    M, K = clamped_matrices(N=N)
    C = zeros(N, N) # no damping
    sys = SecondOrderLinearContinuousSystem(M, C, K)
end

function clamped_forced(; N, F=10e3, E=30e6, A=1)
    M, K = clamped_matrices(N=N)
    C = zeros(N, N) # no damping
    F = vcat(zeros(N-1), F) # the right-most node is excited
    sys = SecondOrderAffineContinuousSystem(M, C, K, F)
end

function clamped(; a=0.0, b=0.0, # ignored if damped = false
                   constant=true,
                   N,            # number of elements
                   homogeneize,  # flag to homogeneize the system
                   damped)       # flag to consider damping C = a*K + b*M

    sys = clamped_forced(N=N)
    M = sys.M; K = sys.K; F = sys.b

    invM = inv(M)    
    ZN = zeros(N, N)
    IN = Matrix(1.0I, N, N)

    if damped
        C = a*K + b*M
        A = [ZN           IN       ;
            -invM*K       -invM*C]
    
    else
        A = [ZN           IN ;
            -invM*K       ZN]
    end

    f0 = vcat(zeros(N), invM * F)
    
    if homogeneize
        n = 2N
        Aext = zeros(n+1, n+1)
        Aext[1:n, 1:n] .= A
        Aext[1:n, n+1] .= f0
        Aext = sparse(Aext)
        S = @system(x' = Aext*x)

    else
        S = @system(x' = A*x + f0)
    end
    
    if !constant
        @assert homogeneize == false

        # model time-varying forcing
        X = Universe(statedim(S))
        ΔF = Interval(0.99, 1.01)
        S = @system(x' = S.A * x + S.c * u, u ∈ ΔF, x ∈ X)
    end

    return S
end

# "nominal" step size
const dtn = 9.88e-7

#=
# problem instances
s5 = clamped_forced(N=5);
s10 = clamped_forced(N=10);
s20 = clamped_forced(N=20);
s30 = clamped_forced(N=30);
s50 = clamped_forced(N=50);
s100 = clamped_forced(N=100);
s200 = clamped_forced(N=200);
s500 = clamped_forced(N=500);
s1000 = clamped_forced(N=1000);
=#

#=

# singleton initial condition
X0(N) = Singleton(zeros(2N))
prob(s) = InitialValueProblem(s, X0(statedim(s)))

# box initial condition
prob_box(s, e) = InitialValueProblem(s, X0(statedim(s)) + BallInf(2*statedim(s), e))

# solve with GLGM06
function solve_glg(s; δ=dtn, NSTEPS=10, model=Forward(inv=false))
    ivp = prob(s)
    alg = GLGM06(δ=δ, approx_model=model)
    sol = solve(ivp, NSTEPS=NSTEPS, alg=alg, homogeneize=true)
end

# solve with LGG09 for nodeidx
function solve_lgg(s; δ=dtn, NSTEPS=100, nodeidx, model=Forward(inv=false))
    ivp = prob(s)
    n = statedim(s)
    eplus = SingleEntryVector(nodeidx, 2n+1, 1.0) 
    eminus = SingleEntryVector(nodeidx, 2n+1, -1.0)

    dirs = [eminus, eplus]
 
    alg = LGG09(δ=δ, template=dirs, approx_model=model) 
    sol = solve(ivp, NSTEPS=NSTEPS, alg=alg, homogeneize=true)
end

# solve with BOX
# usage: sol = solve_box(s10, NSTEPS=100)
function solve_box(s; δ=dtn, NSTEPS=100, model=Forward(inv=false))

    # p = prob(s)
    # ph = homogeneize(normalize(p))
    # @time pd = discretize(ph, dt, model) # ~35 sec for 1000 nodes

    alg = BOX(δ=δ, approx_model=model)
    p = prob(s)
    sbox = solve(p, alg=alg, NSTEPS=NSTEPS, homogeneize=true)
    return sbox
end
=#

# Ec. (4.181) ref Geradin & Rixen (2015)
function _sol_analytic_clamped(x, t, Nmax)

    F = 10E3
    E = 30E6
    dens = 7.3E-4
    L = 200
    A = 1
    
    # reescalar
    nMasses = 1000    
    x = x*L/nMasses

    m = A*dens
    c = sqrt(E*A/m)

    α = 8*F*L/π^2/E/A

    β = π/2/L
    u = sin.( β * x ) * ( 1.0 .- cos.( β * c * t ) )
    v = sin.( β * x ) * β * c * sin.( β * c * t )

    for s in 2:Nmax
        u .+= (-1)^(s-1)/(2*s-1)^2 * sin.( (2*s-1)* β * x ) * ( 1.0 .- cos.( (2*s-1)* β * c * t ) )
        v .+= (-1)^(s-1)/(2*s-1)^2 * sin.( (2*s-1)* β * x ) * ( (2*s-1)* β * c * sin.( (2*s-1)* β * c * t ) )
    end

    return u*α, v*α

end
