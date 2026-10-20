# Compute SIC-POVM fiducials in dimension d by the algebraic method of Appleby, Flammia and Kopp
# (arXiv:2501.03970), with Flammia's Zauner.jl.
#
# A binary quadratic form, of conductor 1 unless that ghost has F_a symmetry (which Zauner.jl does not support),
# gives a "ghost" ψ, whose overlaps (d+1) ϕ'D_pψ / ϕ'ψ with its parity reversal ϕ are real Stark units,
# products of Shintani's double sine function; Newton's method takes ψ to a few thousand bits. On one orbit of
# the Galois group G, LLL then identifies those overlaps exactly over the ring class field H, through a normal
# basis for each cyclic factor of G. An automorphism of H that sends √D to -√D turns them into the overlaps of
# a SIC, up to a shift along G: the right shift is the one that gives a rank-one projector. Alternating
# projections fill in the overlaps off the orbit and bring the fiducial to double precision.
#
# The result is canonical. Of the h automorphisms that send √D to -√D, each giving a Galois conjugate SIC, it
# uses the one that maximizes σ(Σ_orbit u²), u the ghost overlaps. Of the fiducials U_A ψ, A commuting with U_Z's
# symplectic matrix, all eigenvectors of U_Z, it keeps the largest |ψ_k| at the first k where they differ.
#
# Writes fiducial_<d>_stark.txt in Scott and Grassl's layout, as search.py does: an eigenvector of their
# order-3 Zauner unitary U_Z. refine.py takes it to 150 digits. Each step's time goes to stderr.
#
#   julia --project=. stark.jl 197 198 203

using Zauner, LinearAlgebra, Printf
using Nemo: randseed!
using Random: seed!
using GenericFFT: ifft
using QuadGK: quadgk

τ(j, d) = cispi(mod(j * (d + 1), 2d) / d)   # τ^j, τ = -exp(iπ/d): Appleby's D_p is τ^(p1 p2) X^p1 Z^p2, p mod 2d for even d
ω(j, d) = cispi(2mod(j, d) / d)
pq(d) = (0:d-1) * (0:d-1)'
const lap = let t = time()   # lap(s) prints s and the seconds since the last lap
    s -> (println(stderr, "  ", s, ": ", round(Int, time() - t), " s"); flush(stderr); t = time())
end

"ϖ_n(z, t) = Π_{0≤j<n} (1 - e(z + jt)) for real z, t: a phase times sines that obey a three-term recurrence."
function qpoch(z, t, n)
    n < 0 && return (1 - e(z)) / qpoch(z, -t, 1 - n)
    p, a, b, c = one(z), 2sinpi(z), 2sinpi(z + t), 2cospi(t)
    for _ in 1:n
        p, a, b = p * a, b, c * b - a
    end
    return e(n * (2z + (n - 1) * t - 1) / 4) * p
end

"Replaces Zauner._sigma_s (long products, BigFloat quadrature): z first shifted into [0, β), the double sine in Float64."
function sigma(z, β)
    m = floor(Int, z / β)
    z -= m * β
    n = floor(Int, -z) + floor(Int, β / 2)
    w, b = Float64(z + n + 1), Float64(β)
    a = (b + 1) / 2 - w   # g0 is Zauner._g0, but a series below t = 1e-3, where _g0 cancels catastrophically
    g0(t) = t < 1e-3 ? 2a / b * expm1((a^2 - (b^2 + 1) / 4) * t^2 / 6 - (a^4 - (b^4 + 1) / 16) * t^4 / 180) / t^2 : Zauner._g0(w, b, 1.0, t)
    s = quadgk(g0, 0, 1)[1] - 2a / b + quadgk(t -> exp(-t) * Zauner._g(w, b, 1.0, t), 0, Inf)[1]
    return qpoch(z, β, m) * qpoch(z / β, -1 / β, -n) * e((6(z + n)^2 + 6(1 - β) * (z + n) + β^2 - 3β + 1) / (24β)) * exp(BigFloat(s))
end

"Replaces Zauner.ghost (too slow at d ≈ 200): its last σ_s and its ϖ_m(w, β) run to millions of terms; their ratio is σ_s(w + k, β)."
function stark_ghost(F)
    d, A, β = Int(F.d), F.A, F.x
    QA = QuadBin(A[2, 1], A[2, 2] - A[1, 1], -A[1, 2])
    ζ, c = -cispi(BigFloat(1) / d), -e(-BigFloat(rademacher(A)) / 24) / sqrt(BigFloat(d + 1))
    v = Zauner._get_periods(A, β)
    χ = zeros(Complex{BigFloat}, d, 2)
    χ[1, 1] = 1
    for (p1, p2) in Iterators.drop(Iterators.product(0:1, 0:d-1), 1)
        k = Int(divexact(A[1, 2] * p2 - (A[2, 2] - 1) * p1, d))
        s = iseven(d) * (p1 + p2 + p1 * p2)   # Zauner.ghost's sign for even d, (-1)^((1+p1)(1+p2)), over the -1 in c
        ν = (-1)^s * ζ^(-BigFloat(QA(p1, p2)) / (d * (d - 2))) * c * sigma((p2 * β - p1) / d + k, β)
        for i in 2:length(v)-2
            ν *= sigma((v[1] * p2 - v[2] * p1) / (d * v[i+1]), v[i] / v[i+1])
        end
        χ[p2+1, p1+1] = ζ^(p1 * p2) * real(ν)
    end
    χ = ifft(χ, 1)
    return sqrt(abs(χ[1, 1])) * circshift(cumprod(χ[:, 2] ./ χ[:, 1]), 1)
end

"Jacobian of Zauner._ghost_olp_func in closed form: each of its quartic sums has four terms with a given component."
function ghost_jacobian(z::Vector{T}) where {T}
    d, m = 1 + length(z) ÷ 2, 1:length(z)÷2
    ψ = Complex{T}[1; complex.(z[m], z[m.+(d-1)])]
    v(k) = ψ[mod.(k, d).+1]
    N = sum(real(ψ[mod(-k, d)+1] * conj(ψ[k+1])) for k in 0:d-1)
    J = zeros(T, length(z), length(z))
    for (row, n) in enumerate(Zauner._lex_first_inds(d))
        r, p, q = radix(n, [2, d, d])
        A = conj.(v(-m .- p) .* v(-m .- q)) .* v(m .+ (p + q)) .+ v(m .- (p + q)) .* conj.(v(q .- m) .* v(p .- m))
        B = v(-m .- p) .* conj.(v(m .+ (p - q))) .* v(q .- m) .+ v(-m .- q) .* conj.(v(m .+ (q - p))) .* v(p .- m)
        c = 4 * ((p == 0) + (q == 0)) * N / (d + 1)
        J[row, :] = r == 0 ? [real.(A .+ B) .- c .* real.(v(-m)); imag.(B .- A) .- c .* imag.(v(-m))] : [imag.(A .+ B); real.(A .- B)]
    end
    return J
end

"Replaces Zauner.precision_bump (ForwardDiff Jacobians at full precision, too slow here): Newton's method with the
closed-form Jacobian at 256 bits, renewed whenever a step gains under 100 bits."
function bump(ψ, bits)
    z, J, prev = re_im_proj(ψ), nothing, Inf
    while true
        r = setprecision(() -> Zauner._ghost_olp_func(z .+ 0), BigFloat, bits)
        res = maximum(abs, r)
        res > prev / 2^10 && return setprecision(() -> re_im_proj(z), BigFloat, bits)   # at the precision floor
        (J === nothing || res > prev / big(2)^100) && (J = setprecision(() -> lu(ghost_jacobian(z .+ 0)), BigFloat, 256))
        z = setprecision(() -> z .- setprecision(() -> J \ (r .+ 0), BigFloat, 256), BigFloat, bits)
        prev = res
    end
end

"Replaces Zauner.ghost_orbit (a wh call per point, minutes at this precision): (d+1) ϕ'D_p ψ / ϕ'ψ at P, by FFT."
function ghost_units(ψ, P)
    d, ϕ = length(ψ), circshift(reverse(ψ), 1)
    rows = [ifft(conj.(ϕ) .* circshift(ψ, a)) for a in 0:d-1]
    ζ = [(-cispi(one(real(ψ[1])) / d))^-j for j in 0:2d-1]
    return [d * (d + 1) * real(ζ[mod(a * b, 2d)+1] * rows[mod(a, d)+1][mod(b, d)+1]) / real(ϕ'ψ) for (a, b) in P]
end

hankel(x) = [x[mod(a + c, length(x))+1] for a in eachindex(x) .- 1, c in eachindex(x) .- 1]
"T multiplied along each axis m by the Hankel matrix of x[m]: out[a] = Σ_c T[c] Π_m x[m][a_m + c_m]."
hankel_prod(T, x) = foldl((T, m) -> (H = hankel(x[m]); mapslices(t -> H * t, T; dims=m)), 1:ndims(T); init=T)

"Trace forms G[a,b] = Σ_c x[a+c] x[b+c] and C[r,b] = Σ_a x[a+r] x[a] x[a+b]: integral, unlike C G^-1, the product by x[0]."
forms(x) = (H = hankel(x); (H * H, [x[mod(a + r, end)+1] * x[a+1] for r in eachindex(x) .- 1, a in eachindex(x) .- 1] * H))

"Replaces Zauner._dualize (two LLL runs per number, and a third more bits): √D → -√D on integers of H, by Babai rounding to
coordinates in an LLL basis U h, then exactly to coordinates in h. Replaces Zauner.signswitch, which takes the first such
automorphism Hecke lists (that order depends on Hecke's random choices), by the one that makes σ(s) largest."
function conjugator(F, bits, s)
    H = ghostclassfield(F)
    hb, e, (_, r) = lll_basis(maximal_order(H)), real_embeddings(H)[1], is_square_with_sqrt(H(F.D))
    real_at(z) = BigFloat.(real.(evaluation_function(e, bits + 64).(z)))
    h, gs = real_at(hb), [real_at(g.(hb)) for g in automorphism_list(H) if g(r) == -r]
    lap("ring class field H, degree $(length(h))")
    B, U = lll_with_transform(hcat(identity_matrix(ZZ, length(h)), matrix(ZZ, length(h), 1, round.(BigInt, ldexp.(h, bits)))))
    lap("LLL at $bits bits")
    Bf, Ut = BigFloat.(BigInt.(Matrix(B))), transpose(U)
    w = ldexp.((Bf * Bf') \ Bf[:, end], bits)
    coords(x) = BigInt.(vec(Matrix(Ut * matrix(ZZ, length(h), 1, round.(BigInt, x .* w)))))   # in floating point, U g cancels
    σs = [dot(coords(s), g) for g in gs]
    i = sortperm(σs; rev=true)
    length(i) == 1 || σs[i[1]] - σs[i[2]] > 1e-6 * abs(σs[i[1]]) || error("d = $(F.d): two conjugates of s tie")
    return x -> dot(coords(x), gs[i[1]])
end

"Replaces Zauner.pseudonecromancy's core, whose multiplication matrices have denominators that need far more bits: uses trace forms."
function sic_units(σ, u)
    x = ghost_basis(u)[1]
    y = map(x) do x
        G, C = (X -> σ.(X)).(forms(x))
        λ, V = eigen(Float64.(C / G))
        ComplexF64.(G) \ (V[:, 1] * λ[1] / V[1, 1])   # the conjugate normal basis, times G^-1
    end
    return hankel_prod(ComplexF64.(σ.(hankel_prod(u, x))), y)
end

"Π = I/d + Σ_p W_p D_p / (d(d+1)), which inverts W = overlaps(ψ) when Π = |ψ><ψ| is a SIC."
function projector(W)
    d = size(W, 1)
    R = (W .* τ.(pq(d), d)) * ω.(pq(d), d) / (d * (d + 1))
    B = [R[mod(j - k, d)+1, k+1] for j in 0:d-1, k in 0:d-1] + I / d
    return Hermitian((B + B') / 2)
end

"W_p = (d+1) <ψ|D_p|ψ>^* / <ψ|ψ> for every p."
function overlaps(ψ)
    d = length(ψ)
    V = [conj(ψ[mod(k + a, d)+1]) * ψ[k+1] for a in 0:d-1, k in 0:d-1]
    return conj.((d + 1) / norm(ψ)^2 * τ.(pq(d), d) .* (V * ω.(pq(d), d)))
end
sic_error(ψ) = maximum(abs, abs2.(overlaps(ψ))[2:end] .- (length(ψ) + 1)) / (length(ψ) + 1)^2

"Replaces Zauner.shift_search (a d×d matrix per point per shift, too slow here). Shifts go in order of how well W's rows 0
and 1, which give Π's diagonal and subdiagonal, fit |ψ_k|^2 and ψ_(k+1) ψ_k^*; one passes if Π's top eigenvalue is near the
fraction of W known. Alternating projections then fill in the rest, the job of rank_one_least_squares."
function sic(uu, P, Ks, d)
    # relabel the points by M = N^-1, N = [v Kv] of det 1: M K M^-1 = [0 -1; 1 -1], U_Z's symplectic matrix, for K the first
    # power of the stabilizer generator with trace -1; all mod dd = 2d for even d, where Clifford unitaries need matrices mod 2d
    dd = 2^iseven(d) * d
    N = first(N for K in Ks if mod(K[1] + K[4] + 1, dd) == 0
              for N in ([[a, b] K * [a, b]] for a in 0:dd-1, b in 0:dd-1) if mod(N[1] * N[4] - N[2] * N[3], dd) == 1)
    M = [N[4] -N[3]; -N[2] N[1]]
    pts = map(CartesianIndices(P)) do i
        p = mod.(M * P[i], dd)
        q = mod.(p, d)
        (q .+ 1, CartesianIndex(Base.tail(Tuple(i))), (-1)^div(p[1] * p[2] - q[1] * q[2], d))   # D_p = ±D_q
    end
    at(g, s) = conj(uu[CartesianIndex(mod.(Tuple(g) .- Tuple(s), size(uu)) .+ 1)])
    shifts = vec(CartesianIndices(uu))
    A = [zeros(ComplexF64, length(shifts), d) for _ in 1:2]
    for (p, g, σ) in filter(q -> q[1][1] <= 2, pts), (i, s) in enumerate(shifts)
        A[p[1]][i, p[2]] = σ * at(g, s) * τ((p[1] - 1) * (p[2] - 1), d) / (d * (d + 1))
    end
    D, S = 1 / d .+ real.(A[1] * ω.(pq(d), d)), A[2] * ω.(pq(d), d)
    fit = round.(d * vec(sum(abs.(abs2.(S) .- D .* circshift(D, (0, -1))), dims=2)), digits=6)   # rounded: ties go by index
    known = falses(d, d)
    foreach(q -> known[q[1]...] = true, pts)
    known[1, 1] = true   # W_0 = 0, for the I/d in Π
    for s in shifts[sortperm(fit)]
        W = zeros(ComplexF64, d, d)
        foreach(((p, g, σ),) -> W[p...] = σ * at(g, s), pts)
        λ, V = eigen(projector(W))
        λ[end] > count(known) / d^2 / 2 || continue
        ψ = V[:, end]
        for _ in 1:200   # hold W where known, else only |W| = √(d+1); release all once nearly a SIC
            C = overlaps(ψ)
            ψ = normalize(projector(ifelse.(known .& (sic_error(ψ) > 1e-9), W, sqrt(d + 1) * C ./ abs.(C))) * ψ)
        end
        sic_error(ψ) < 1e-13 && return ψ
    end
    error("d = $d: no shift gives a SIC")
end

"Of the fiducials U_A ψ, A in the centralizer of U_Z's symplectic matrix (mod 2d for even d), each an eigenvector of U_Z and
each one the shift search may find, the one with the largest |ψ_k| at the first k where they differ by more than 1e-9."
function canonical(ψ)
    d, W, best = length(ψ), overlaps(ψ), ψ
    dd = 2^iseven(d) * d
    larger(a, b) = (k = findfirst(i -> abs(abs(a[i]) - abs(b[i])) > 1e-9, eachindex(a)); k !== nothing && abs(a[k]) > abs(b[k]))
    for x in 0:dd-1, y in 0:dd-1
        mod(x^2 - x * y + y^2, dd) == 1 || continue
        V = zeros(ComplexF64, d, d)
        for a in 0:d-1, b in 0:d-1   # U_A ψ has W_p at A p, and D_(A p) = ±D_q for q = A p mod d
            p = mod.([x -y; y x-y] * [a, b], dd)
            q = mod.(p, d)
            V[q[1]+1, q[2]+1] = (-1)^div(p[1] * p[2] - q[1] * q[2], d) * W[a+1, b+1]
        end
        V[1, 1] = 0
        Π = projector(V)   # |φ><φ|: φ from its largest column, with φ_0 real and positive
        k = argmax(real(diag(Π)))
        φ = Π[:, k] * conj(Π[1, k]) / abs(Π[1, k]) / sqrt(real(Π[k, k]))
        larger(φ, best) && (best = φ)
    end
    return best
end

"The form of the smallest conductor q without F_a symmetry (which Zauner.jl does not support) that Zauner.jl's table, which
stops at d = 256, lists first for d: in the principal class, the reduced form whose stabilizer has the shortest word."
function Q(d)
    Δ, f = coredisc((d + 1) * (d - 3))
    q = first(q for q in sort(divisors(f)) if !(mod(d, 9) == 3 && mod(f ÷ q, 3) == 0))   # as Zauner.has_fa_symmetry
    return minimal_hj_stabilizer(reduced_hj_orbit(quadbinid(q^2 * Δ)), d)
end

"A SIC fiducial from the admissible tuple (d, Q(d))."
function fiducial(d)
    Qd = Q(d)
    randseed!(0); seed!(0)   # Hecke's random choices (FLINT's RNG and Julia's), fixed per d; the result does not depend on them
    F, dd = AdmissibleTuple(d, Qd), 2^iseven(d) * d
    S, G = stabilizer_elements(F), galois_elements(d, galois_normal_form(F)...)
    orbit(p) = [mod.(Int.(s * g) * p, dd) for s in S, g in G]
    seeds = Iterators.flatten((([1, 0], [0, 1], [1, 1]), ([a, b] for a in 0:dd-1, b in 0:dd-1)))
    free = Iterators.filter(p -> allunique(orbit(p)), seeds)
    isempty(free) && error("d = $d: no free orbit")
    P = orbit(first(free))
    lap("$(repr("text/plain", F)), $(length(P)) points")
    ψ = bump(setprecision(() -> stark_ghost(F), BigFloat, 128), 256)   # in double precision, ghost_basis can mistake 0 for a basis
    u = setprecision(() -> ghost_units(ψ, selectdim(P, 1, 1)), BigFloat, 256)
    x = ghost_basis(u)[1]
    top = maximum(maximum(abs, X) for X in (hankel_prod(u, x), Iterators.flatten(forms.(x))...))
    bits = 2F.h * (ceil(Int, log2(top)) + 40)   # LLL needs about this many bits to identify these numbers in H
    lap("ghost")
    sic_like(uu) = maximum(abs, abs2.(uu) .- (d + 1)) < 1e-3   # conjugated units have modulus √(d+1)
    uu = nothing
    for b in (bits, 2bits, 4bits)   # more, when these numbers are larger in the other embeddings of H
        ψ = bump(ψ, b + 200)
        lap("Newton's method to $(b + 200) bits")
        uu = setprecision(BigFloat, b + 200) do
            u = ghost_units(ψ, selectdim(P, 1, 1))
            return sic_units(conjugator(F, b, sum(u .^ 2)), u)
        end
        lap("√D → -√D")
        sic_like(uu) && break
    end
    sic_like(uu) || error("d = $d: conjugated units do not have modulus √(d+1), even at $(4bits) bits")
    ψ = sic(uu, P, [Int.(mod.(F.L^i, dd)) for i in 1:2F.k], d)   # L^k = I mod d, so L^2k = I mod 2d
    lap("shift search")
    ψ = canonical(ψ)
    lap("canonical representative")
    return ψ, F
end

"x as '+.<17 digits>e<exponent>', as search.py's scott_number writes it."
function scott(x)
    iszero(x) && return "+." * "0"^17 * "e0"
    m, e = split(@sprintf("%.16e", abs(x)), 'e')
    return (x > 0 ? "+." : "-.") * replace(m, "." => "") * "e$(parse(Int, e) + 1)"
end

for d in parse.(Int, ARGS)
    t = @elapsed ψ, F = fiducial(d)
    ψ = ψ * conj(ψ[1]) / abs(ψ[1]) / norm(ψ)   # Scott and Grassl's gauge: unit norm, ψ_0 real and positive
    U = cispi((d - 1) / 12) / sqrt(d) * cispi.(mod.(2pq(d) .+ (d + 1) * (0:d-1) .^ 2, 2d) / d)   # U_Z, as in verify.py, for every d
    residual = norm(U * ψ - (ψ' * U * ψ) * ψ)
    residual < 1e-12 || error("d = $d: not an eigenvector of U_Z (residual $residual)")
    write("fiducial_$(d)_stark.txt", join(scott(real(z)) * scott(k == 1 ? 0.0 : imag(z)) * "i\n" for (k, z) in enumerate(ψ)))
    @printf("d = %d, form %s: fiducial_%d_stark.txt in %.0f s; SIC error %.1e, |U_Z ψ - λψ| = %.1e\n", d, F.Q, d, t, sic_error(ψ), residual)
    flush(stdout)
end
