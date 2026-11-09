#Computer assisted proof of solutions for the 2D Swift Hohenberg equation and (in)stability with continuous symmetries.
# (I + Δ)²u + μu - γu² + u³ = 0
# The following code computes the solution and rigorously proves the results given in section 5 of
# "Stability of patterns in variational systems with continuous symmetries: a computer-assisted approach"  J.B. van den Berg and D. Blanco

# We computed the approximate solution by placing random spots in a domain of size d and running Newton's method. We provide the data here.
# From this we can check if the proof of the solution is verified or not and its (in)stability using Gershgorin circles. We essentially prove Theorems 5.1, 5.2, 5.3, and 5.4.
# To choose which theorem to prove, set sol on line 128: sol = 1, 2, 3 or 4 proves Theorem 5.1, 5.2, 5.3 or 5.4.

#####################################################################################################################################################################

# Needed packages
using RadiiPolynomial, LinearAlgebra, DelimitedFiles

#####################################################################################################################################################################


#################################### Helper functions #################################################

# Reads n1 n2 real(a_n) imag(a_n) and returns the complex coefficient vector.
# The coefficient matrix a[n1,n2] is stored column by column.
function load_coefficients(name)
    C = readdlm(joinpath(@__DIR__, name * ".txt"), Float64;
                comments=true, comment_char='%')
    return complex.(C[:,3], C[:,4])
end

# Builds an operator 𝒱 such that ||u||_1,ν = ||𝒱u||_1.
function _build_𝒱(ν,N)
    𝒱 = interval.(zeros((2N+1)^2))
    r = 1
    for n₁ = -N:N 
        for n₂ = -N:N 
            𝒱[r] = ν^(interval(abs(n₁)+abs(n₂)))
            r += 1 
        end
    end
    return 𝒱
end

# Computes the laplacian as a vector for more efficient operations.
function _laplacian_vectorf(N,d)
    Δ = interval.(zeros((2N+1)^2))
    r = 0
    for n₁ = -N:N 
        for n₂ = -N:N 
            r += 1
            Δ[r] = -(exact(n₁)*(π/d))^2 - (exact(n₂)*(π/d))^2 
        end
    end
    return Δ 
end

# The map F̄ defined for intervals. 
function F̄i!(F̄i,L,x,μ,γ,∂₁u₀,∂₁ū₀,∂₂u₀,∂₂ū₀)
    β₁,β₂,u = eachcomponent(x)
    component(F̄i,1) .= sum(u.*∂₁ū₀)
    component(F̄i,2) .= sum(u.*∂₂ū₀)
    project!(component(F̄i,3),L.*u - γ*u^2 + u^3 + β₁[1]*∂₁u₀ + β₂[1]*∂₂u₀)
    return F̄i
end

# Builds the operator 𝒟 defined in Lemma 4.6 of Section 4.3.
function _build_𝒟(N,d)
    𝒟 = interval.(zeros((2N+1)^2))
    r = 0
    for n₁ = -N:N 
        for n₂ = -N:N 
            r += 1
            𝒟[r] = interval(abs(n₁)+abs(n₂))*(π/d)
        end
    end
    𝒟[div((2N+1)^2-1,2)+1] = interval(1)
    return 𝒟
end

# Computes the Gershgorin radii of a disc.
function _gershgorin_radiusi(D)
    _rows = size(D)[1]
    R = interval.(zeros(_rows))
    for i = 1:_rows 
        for k = 1:(i-1)
            R[i] += abs(D[i,k])
        end
        for k = (i+1):_rows 
            R[i] += abs(D[i,k])
        end
    end
    return R 
end

# Computes the Gershgorin discs.
function _gershgorin_discsi(R,E)
    l = length(R)
    _discs = interval.(complex.(zeros(l)))
    for i = 1:l 
        _discs[i] = real(E[i]) + interval(-sup(R[i]), sup(R[i]))
    end
    return _discs 
end

# Verifies whether we have a stable or unstable solution.
function _gershgorin_verify(_discs,Ri)
    l = length(_discs)
    _negative = 0 # the number of negative eigenvalues of Ā^N
    for k = 1:l
        if inf(real(_discs[k])) <= 0 <= sup(real(_discs[k]))
            println("A Gershgorin disc meets the imaginary axis, so no conclusion is possible")
            return
        end
        if sup(real(_discs[k])) < 0
            _negative += 1
        end
    end
    if _negative == 2
        println("Solution is stable")
    else
        println("Solution is unstable, number of unstable directions: ", _negative - 2)
    end
    k = _negative + 1
    println("Discs 1, 2 and ", k, " (centre, radius): ", (mid(real(_discs[1])), sup(Ri[1])), (mid(real(_discs[2])), sup(Ri[2])), (mid(real(_discs[k])), sup(Ri[k])))
end

################### PROOF OF SOLUTIONS : MAIN CODE #################################################################################################################################################
# To run the code, click the run button in the terminal.
# Below is the data for four solutions. To prove one of them, set sol below.
# Make sure to download the necessary file that includes the data you wish to prove!
sol = 1 # the solution to prove: 1, 2, 3 or 4
if sol == 1 # Theorem 5.1
N₀ = 16
N = 10
d = interval(7)
ν = interval(1.1)
r₁_star = interval(5e-8)
r₂_star = interval(5e-8) 
u₀_coeffs = load_coefficients("u_stab1")
elseif sol == 2 # Theorem 5.2
N₀ = 27
N = 20
d = interval(12)
ν = interval(1.09) 
r₁_star = interval(1e-10)
r₂_star = interval(1e-8)
u₀_coeffs = load_coefficients("u_stab2")
elseif sol == 3 # Theorem 5.3
N₀ = 18
N = 13
d = interval(9)
ν = interval(1.1)#Sol 1
r₁_star = interval(1e-9)
r₂_star = interval(6e-8) 
u₀_coeffs = load_coefficients("u_unstab1")
elseif sol == 4 # Theorem 5.4
N₀ = 35
N = 27
d = interval(15)
ν = interval(1.08) 
r₁_star = interval(4e-9)
r₂_star = interval(4e-9) 
u₀_coeffs = load_coefficients("u_unstab2")
end

# Creating needed objects
μ = interval(0.25)
γ = interval(1.6)
fourier = Fourier(N,π/d) ⊗ Fourier(N,π/d) # space for operators 
fourier0 = Fourier(N₀,π/d) ⊗ Fourier(N₀,π/d) # space for sequences
pfourier = ParameterSpace() × ParameterSpace() × fourier
w₀ = Sequence(ParameterSpace() × ParameterSpace() × fourier0, [interval(0) ; interval(0) ; interval.(u₀_coeffs)])
u₀ = component(w₀,3)
∂₁ū₀ = conj(Derivative(1,0)*u₀)
∂₂ū₀ = conj(Derivative(0,1)*u₀)
∂₁u₀ = Derivative(1,0)*u₀
∂₂u₀ = Derivative(0,1)*u₀

# Verify the unfolding parameters are 0. 
M₀ = [norm(∂₁u₀,2)^2 sum(∂₁ū₀.*∂₂u₀) ; sum(∂₂ū₀.*∂₁u₀) norm(∂₂u₀,2)^2 ]
κ = π/(exp(interval(1))*d*log(ν))
M_κ = κ*r₂_star*[norm(∂₁u₀,2) norm(∂₂u₀,2); norm(∂₁u₀,2) norm(∂₂u₀,2)]
M_check = M₀ + M_κ*interval(-1,1)
@show 0 < inf(real(M_check[1,1]*M_check[2,2] - M_check[1,2]*M_check[2,1]))

# Choice of Banach Space
ℓ¹_ν = Ell1(GeometricWeight(ν))
X_ν = NormedCartesianSpace((ℓ¹(),ℓ¹(),ℓ¹_ν),ℓ¹())

##################################### Construction of Ā ######################
Δ = _laplacian_vectorf(N,d)
L = (exact(1) .+ Δ).^2 .+ μ

# Building the approximate inverse.
v₀ = -exact(2)*γ*u₀ + exact(3)*u₀^2
𝕧₀ᴺ = project(Multiplication(v₀),fourier,fourier,Complex{Interval{Float64}})
𝕧̄₀ᴺ = LinearOperator(pfourier, pfourier, complex.(interval.(zeros(dimension(pfourier),dimension(pfourier)))))
project!(component(𝕧̄₀ᴺ,3,3),𝕧₀ᴺ)
component(𝕧̄₀ᴺ,1,3)[1,:] .= conj(project(∂₁u₀,fourier))
component(𝕧̄₀ᴺ,2,3)[1,:] .= conj(project(∂₂u₀,fourier))
project!(component(𝕧̄₀ᴺ,3,1), ∂₁u₀)
project!(component(𝕧̄₀ᴺ,3,2), ∂₂u₀)
L̄ = LinearOperator(pfourier,pfourier, complex.(interval.(zeros(dimension(pfourier),dimension(pfourier)))))
project!(component(L̄,3,3), LinearOperator(fourier,fourier,Diagonal(L)))

𝒱 = _build_𝒱(ν,N)
𝒱⁻¹ = interval(ones((2N+1)^2))./𝒱
Ā = inv(mid.(L̄+𝕧̄₀ᴺ)) # note that we do the inverse using floats as this is an approximate inverse.
Ā = interval.((Ā + Ā')/2) # We symmetrize Ā so that it is self-adjoint.
norm_Ā = opnorm([interval(1) ; interval(1) ; 𝒱].*Ā.*[interval(1) ; interval(1) ; 𝒱⁻¹]',1)
ζ_N = exact(N+1)*(π/d)
L_N = abs((ζ_N^2-exact(1))^2 + μ)
@show sup(1 + abs(μ)) <= inf(ζ_N^2)

##################################### Computation of Y₀ⁱ for i = 1,2 ######################
Δ_long = _laplacian_vectorf(N₀,d)
L_long = (exact(1) .+ Δ_long).^2 .+ μ
G = -γ*u₀^2 + u₀^3
ĀF̄ = Ā*project(F̄i!(similar(w₀),L_long,w₀,μ,γ,∂₁u₀,∂₁ū₀,∂₂u₀,∂₂ū₀),pfourier)
Y₀¹ = norm(component(ĀF̄,1:2),1)
Y₀² = norm(𝒱.*component(ĀF̄,3),1) + exact(1)/L_N*norm(_build_𝒱(ν,3N₀).*(L_long.*(u₀-project(u₀,fourier))+G-project(G,fourier)),1)
@show Y₀¹,Y₀²

############################## Computation of W₂₂ⁱ for i = 1,2 #############################
𝒱_N_N0 = _build_𝒱(ν,N+N₀) # Builds 𝒱 with N + N₀ coefficients in order to take into account the projection on q.
𝒱_N₀ = _build_𝒱(ν,N₀) #We also need 𝒱 with N₀ coefficients to compute ||q||_1,ν.
𝒱⁻¹_N_N0 = interval(ones((2*(N+N₀)+1)^2))./𝒱_N_N0
fourier_N0_N = Fourier(N+N₀,π/d) ⊗ Fourier(N+N₀,π/d)
q = -exact(2)*γ + exact(6)*u₀
𝕢 = project(Multiplication(q),fourier_N0_N,fourier,Complex{Interval{Float64}})
Ā₁₂ = component(Ā,1:2,3)
Ā₂₂ = component(Ā,3,3)
W₂₂¹(r₂) = opnorm(Ā₁₂*𝕢.*𝒱⁻¹_N_N0',1) + exact(3)*opnorm(Ā₁₂.*𝒱⁻¹',1)*r₂
W₂₂²(r₂) = (opnorm(𝒱.*Ā₂₂*𝕢.*𝒱⁻¹_N_N0',1) + interval(1)/L_N*norm(𝒱_N₀.*q,1)) + exact(3)*(opnorm(𝒱.*Ā₂₂.*𝒱⁻¹',1) + exact(1)/L_N)*r₂

############################# Computation of Z₁ⁱ for i = 1,2 ###########################
v₀ᴺ = project(v₀,fourier)
∂₁u₀ᴺ = project(∂₁u₀,fourier)
∂₂u₀ᴺ = project(∂₂u₀,fourier)
∂u₀ᴺ = LinearOperator(ParameterSpace() × ParameterSpace(),fourier, interval.(complex.(zeros((2N+1)^2,2))))
component(∂u₀ᴺ,1) .= coefficients(∂₁u₀ᴺ)
component(∂u₀ᴺ,2) .= coefficients(∂₂u₀ᴺ)
Z₁¹ᴺ = opnorm(UniformScaling(interval(1))-Ā₁₂*∂u₀ᴺ,1)

Z₁²ᴺ = opnorm(𝒱.*Ā₂₂*∂u₀ᴺ,1)

Z₁¹_∞ =  opnorm(Ā₁₂.*𝒱⁻¹',1)*(norm(𝒱_N₀.*(∂₁u₀ᴺ - ∂₁u₀),1) + norm(𝒱_N₀.*(∂₂u₀ᴺ - ∂₂u₀),1))
Z₁²_∞ = (opnorm(𝒱.*Ā₂₂.*𝒱⁻¹',1) + exact(1)/L_N)*(norm(𝒱_N₀.*(∂₁u₀ᴺ - ∂₁u₀),1) + norm(𝒱_N₀.*(∂₂u₀ᴺ - ∂₂u₀),1))

Z₁¹ = Z₁¹ᴺ + Z₁¹_∞
Z₁² = Z₁²ᴺ + Z₁²_∞
@show Z₁¹,Z₁²
############################# Computation of Z₂ⁱ for i = 1,2 ###########################
𝒱_2N = _build_𝒱(ν,2N) # We need 𝒱 of size 2N to handle the projection.
𝒱⁻¹_2N = interval(ones((4N+1)^2))./𝒱_2N
fourier2 = Fourier(2N,π/d) ⊗ Fourier(2N,π/d)
𝕧₀ᴺ_2N = project(Multiplication(v₀),fourier2,fourier,Complex{Interval{Float64}})
Ā₁₁ = component(Ā,1:2,1:2)
L_full = LinearOperator(fourier,fourier, Diagonal(L))
Z₂¹ᴺ = opnorm((Ā₁₁*(∂u₀ᴺ') + Ā₁₂*(L_full+𝕧₀ᴺ_2N)).*𝒱⁻¹_2N',1)

Ā₂₁ = component(Ā,3,1:2)
Z₂²ᴺ = opnorm(𝒱.*(UniformScaling(interval(1)) - Ā₂₁*(∂u₀ᴺ') - Ā₂₂*(L_full+𝕧₀ᴺ_2N)).*𝒱⁻¹_2N',1) + exact(1)/L_N * norm(𝒱.*v₀ᴺ,1)

𝒱_2N₀ = _build_𝒱(ν,2N₀) #We need 𝒱 of 2N₀ coefficients since v₀ has 2N₀ coefficients.
Z₂¹_∞ = opnorm(Ā₁₁,1)*(norm(∂₁u₀ᴺ - ∂₁u₀,2) + norm(∂₂u₀ᴺ - ∂₂u₀,2)) + opnorm(Ā₁₂.*𝒱⁻¹',1)*norm(𝒱_2N₀.*(v₀ᴺ - v₀),1)
Z₂²_∞ = opnorm(𝒱.*Ā₂₁,1)*(norm(∂₁u₀ᴺ - ∂₁u₀,2) + norm(∂₂u₀ᴺ - ∂₂u₀,2)) + (opnorm(𝒱.*Ā₂₂.*𝒱⁻¹',1) + exact(1)/L_N)*norm(𝒱_2N₀.*(v₀ᴺ - v₀),1)

Z₂¹ = Z₂¹ᴺ + Z₂¹_∞ 
Z₂² = Z₂²ᴺ + Z₂²_∞
@show Z₂¹,Z₂²

################################### Computation of Z_qⱼⁱ for i,j = 1,2 #############################
𝒟 = _build_𝒟(N,d) # The operator 𝒟 defined in Lemma 4.6.

Z_q₁¹(r₂) = opnorm(Ā₁₂.*(𝒟.*𝒱⁻¹)',1)*r₂ 
Z_q₂¹(r₂) = opnorm(Ā₁₁,1)*exact(2)*κ*r₂
Z_q₁²(r₂) = (opnorm(𝒱.*(Ā₂₂.*𝒟').*𝒱⁻¹',1) + sqrt(interval(2))*ζ_N/L_N)*r₂ 
Z_q₂²(r₂) = opnorm(𝒱.*Ā₂₁,1)*exact(2)*κ*r₂

# We name the conditions to check p₁ and p₂ for convenience.
p₁(r₁,r₂) = Y₀¹ + Z₁¹*r₁ + r₂*(Z₂¹ + interval(0.5)*W₂₂¹(r₂)*r₂)
p₂(r₁,r₂) = Y₀² + Z₁²*r₁ + r₂*(Z₂² + interval(0.5)*W₂₂²(r₂)*r₂) 

# Checking the conditions of Theorems 2.2 and 3.2.
H₁¹(r₁,r₂) = Z₁¹
H₁²(r₁,r₂) = Z₁² 
H₂¹(r₁,r₂) = Z₂¹ + W₂₂¹(r₂)*r₂
H₂²(r₁,r₂) = Z₂² + W₂₂²(r₂)*r₂
@show sup(p₁(r₁_star,r₂_star)) < inf(r₁_star)
@show sup(p₂(r₁_star,r₂_star)) < inf(r₂_star)
δ₁ = H₁¹(r₁_star,r₂_star) + Z_q₁¹(r₂_star) + H₁²(r₁_star,r₂_star) + Z_q₁²(r₂_star)
δ₂ = H₂¹(r₁_star,r₂_star) + Z_q₂¹(r₂_star) + H₂²(r₁_star,r₂_star) + Z_q₂²(r₂_star)
δ = max(δ₁,δ₂) # δ from equation 3.3
@show sup(δ) < 1
################################### Computing the Gershgorin Discs #############################
V = interval.(eigvecs(coefficients(mid.(Ā)))) # The matrix V of approximate eigenvectors we will use for our psuedo-diagonalization.
D = inv(V)*coefficients(Ā)*V # Note that we compute an interval inverse of V to rigorously enclose V⁻¹.
Di = interval.(D)
Ri = _gershgorin_radiusi(Di)
_discs = _gershgorin_discsi(Ri,diag(Di))
_gershgorin_verify(_discs,Ri) # We now verify if the solution is (un)stable. 

# We now compute the spectral gap from Theorem 3.8.
_rho = interval(maximum(sup.(abs.(real.(_discs))))) # an upper bound for the spectral radius of Ā^N, from the Gershgorin discs
_spectral_gap = inf((exact(1)-δ)/((exact(1)-δ)*max(_rho, exact(1)/L_N) + δ*max(norm_Ā, exact(1)/L_N)))
@show _spectral_gap
