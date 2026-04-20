
# Rotate using Rodrigue's formula
@inline function rotation_matrix(Θ::AbstractVecOrMat{T}) where T
    
    Θ_norm = norm(Θ)
    if Θ_norm > 10*eps(T)
        sinΘ = sin(Θ_norm)
        SΘ = skew(Θ)
        R = ID3 + sinΘ/Θ_norm*SΘ +  2*(sin(Θ_norm/2)/Θ_norm)^2 * SΘ*SΘ
        return R
    else
        return SMatrix{3,3,T,9}(1,0,0,0,1,0,0,0,1)
    end

end

@inline function local_R⁰(x₁, x₂)

    v₁ = x₂ - x₁
    l = norm(v₁)
    v₁ = v₁/l
    e₂ = @SVector [0,1,0]
    e₃ = @SVector [0,0,1]

    if ( v₁[1] > 1e-8*l ) || ( v₁[2] > 1e-8*l )
        aux = cross(e₃, v₁)
        v₂ = aux/norm(aux)
    else
        T = eltype(v₁)
        v₂ = T.(e₂)
    end
    
    v₃ = cross(v₁, v₂)

    Rₑ = [v₁ v₂ v₃]

    return Rₑ

end

@inline function local_Rₑ_and_aux(x₁, x₂, R₁, R₂, Rₑ⁰E₂, lₙ)

    v₁ = (x₂ - x₁)/lₙ
    p₁ = R₁ * Rₑ⁰E₂
    p₂ = R₂ * Rₑ⁰E₂
    p = (p₁+p₂)/2
    v₃ = cross(v₁, p)
    v₃ = v₃ / norm(v₃)
    v₂ = cross(v₃, v₁)

    Rₑ = [v₁ v₂ v₃]

    ru₁ = -v₁'
    ru₂ = v₁'

    Rₑ¹²ᵀ = Rₑ[:, SOneTo(2)]'

    q₁, q₂ = Rₑ¹²ᵀ*p
    q¹¹, q¹² = Rₑ¹²ᵀ*p₁
    q²¹, q²² = Rₑ¹²ᵀ*p₂
    
    η = q₁/q₂
    η¹¹ = q¹¹/q₂
    η¹² = q¹²/q₂
    η²¹ = q²¹/q₂
    η²² = q²²/q₂

    Gᵀu₁ = @SMatrix [0 0 η/lₙ; 0 0 1/lₙ; 0 -1/lₙ 0]
    Gᵀu₂ = -Gᵀu₁
    GᵀΘ₁ = @SMatrix [η¹²/2 -η¹¹/2 0; 0 0 0; 0 0 0]
    GᵀΘ₂ = @SMatrix [η²²/2 -η²¹/2 0; 0 0 0; 0 0 0]

    D₃ = (ID3 - v₁*v₁')/lₙ

    return Rₑ, ru₁, ru₂, η, Gᵀu₁, GᵀΘ₁, Gᵀu₂, GᵀΘ₂, D₃
end

@inline function toangle(R::AbstractMatrix{T}) where T

    if abs((tr(R)-1)/2) < 1
        norm_v = max(acos((tr(R)-1)/2), 10*eps(T))
    else
        norm_v = 10*eps(T)
    end

    n_v = (1/(2*sin(norm_v))) * @SVector [R[3,2]-R[2,3], R[1,3]-R[3,1], R[2,1]-R[1,2]]

    return norm_v*n_v

end

@inline function skew(vec)

    return @SMatrix [0       -vec[3] vec[2] ;
                     vec[3]  0       -vec[1];
                     -vec[2] vec[1]  0      ]

end

# Get the inverse of skew symmetric matrix from toangle
@inline function Tₛ⁻¹(Θ::AbstractVector{T}) where T

    Θ_norm = norm(Θ)

    if Θ_norm < 10*eps(T)
        Tₛ⁻¹ = SMatrix{3,3,T,9}(1,0,0,0,1,0,0,0,1)
    else
        SΘ = skew(Θ)
        sinΘ2, cosΘ2 = sincos(Θ_norm/2)
        Tₛ⁻¹ = ID3 - SΘ/2 + (sinΘ2-Θ_norm/2*cosΘ2)/(sinΘ2*Θ_norm^2)*SΘ*SΘ
    end

    return Tₛ⁻¹

end

@inline function Tₛ(Θ::AbstractVector{T}) where T

    Θ_norm = norm(Θ)

    if Θ_norm < 10*eps(T)
        Tₛ = SMatrix{3,3,T,9}(1, 0, 0, 0, 1, 0, 0, 0, 1)
    else
        SΘ = skew(Θ)
        Tₛ = ID3 + 2*(sin(Θ_norm/2)/Θ_norm)^2*SΘ + (1-sin(Θ_norm)/Θ_norm)/Θ_norm^2*(SΘ*SΘ)
    end

    return Tₛ

end

@inline function compute_Kᵥ(Θ::AbstractVector{T}, v) where T

    Θnorm = norm(Θ)

    if Θnorm < 10*eps(T)
        Kᵥ = skew(v)/2
    else
        sinΘ, cosΘ = sincos(Θnorm)
        sinΘ2 = sin(Θnorm/2)
        aux = (2*sinΘ2/Θnorm)^2

        u = Θ/Θnorm
        uvᶜ = cross(u, v)
        uvᵈ = dot(u,v)
        UU = u*u'
        VU = v*u'
        UV = u*v'
    
        Kᵥ =  (cosΘ-sinΘ/Θnorm)/Θnorm * (VU-uvᵈ*UU) + 
              (1-sinΘ/Θnorm)/Θnorm * (UV - 2*uvᵈ*UU + uvᵈ*ID3) -
              (sinΘ/Θnorm-aux) * (uvᶜ*u') +
              aux * skew(v)/2
    
    end

    return Kᵥ

end

#  Compute Kint matrix
@inline function K̄ⁱⁿᵗ_beam(E, G, Iₒ, A, I₂₂, I₃₃, l₀)
    
    K̄ⁱⁿᵗū = A*E/l₀
    K̄ⁱⁿᵗΘ̅ = Diagonal(@SVector [G*Iₒ/l₀, 4*E*I₃₃/l₀, 4*E*I₂₂/l₀])
    K̄ⁱⁿᵗΘ̅Θ̅ = Diagonal(@SVector [-G*Iₒ/l₀, 2*E*I₃₃/l₀, 2*E*I₂₂/l₀])
    
    return K̄ⁱⁿᵗū, K̄ⁱⁿᵗΘ̅, K̄ⁱⁿᵗΘ̅Θ̅
    
end

@inline function compute_η_μ(Θ̄::AbstractVector{T}) where T

    Θ = norm(Θ̄)

    if Θ<10*eps(T)
        η = T(1/12)
        μ = T(1/360)
    else
        sinΘ, cosΘ = sincos(Θ)
        sinΘ2 = sin(Θ/2)
        η = ((2*sinΘ)-Θ*(1+cosΘ))/(2*(Θ^2)*sinΘ)
        μ = (Θ*(Θ+sinΘ)-8*(sinΘ2)^2)/(4*(Θ^4)*(sinΘ2)^2)
    end
    
    return η, μ

end

@inline function compute_K̄ₕ(Θ̅, M̄, Tₛ⁻¹Θ̅, η, μ)

    Θ̅M̄ᵀ = Θ̅*M̄'
    M̄Θ̅ᵀ = Θ̅M̄ᵀ'

    SΘ̅ = skew(Θ̅)
    SM̄ = skew(M̄)

    K̄ₕ =  ( η*(Θ̅M̄ᵀ - 2*M̄Θ̅ᵀ + dot(Θ̅, M̄)*ID3) + μ*(SΘ̅*SΘ̅*M̄Θ̅ᵀ) - SM̄/2 ) * Tₛ⁻¹Θ̅

    return K̄ₕ

end

@inline function Pmatrices(N₁, N₂, N₃, N₄, N₅, N₆, lₙ, η, η₁₁, η₁₂, η₂₁, η₂₂)

    P₁P¹ = @SMatrix [0 0 0; 0 (N₃+N₄)/lₙ 0; 0 0 (N₃+N₄)/lₙ]
    P₁P² = @SMatrix [0 0 0; 0 0 N₃; 0 -N₃ 0]
    P₁P³ = -P₁P¹
    P₁P⁴ = @SMatrix [0 0 0; 0 0 N₄; 0 -N₄ 0]

    P₂P¹ = @SMatrix [0 0 -η*(N₁+N₂)/lₙ; 0 0 -(N₅+N₆)/lₙ; 0 (N₅+N₆)/lₙ 0]
    P₂P² = @SMatrix [-(N₂*η₁₂)/2-N₁*(η₁₂/2 - 1) (η₁₁*(N₁+N₂))/2 0; 0 N₅ 0; 0 0 N₆]
    P₂P³ = -P₂P¹
    P₂P⁴ = @SMatrix [-(N₁*η₂₂)/2-N₂*(η₂₂/2 - 1) (η₂₁*(N₁+N₂))/2 0; 0 N₅ 0; 0 0 N₆]


    return P₁P¹, P₁P², P₁P³, P₁P⁴, P₂P¹, P₂P², P₂P³, P₂P⁴

end

function compute(u₁::AbstractVector{T}, u₂, R₁, R₂, ΔR₁, ΔR₂, u̇₁, u̇₂, ẇ₁, ẇ₂, ü₁, ü₂, ẅ₁, ẅ₂, constants, exact=true, isdynamic=true) where T

    # Superscript ¹ means matrix or vector associated to u₁
    # Superscript ² means matrix or vector associated to Θ₁
    # Superscript ³ means matrix or vector associated to u₂
    # Superscript ⁴ means matrix or vector associated to Θ₂

    init, gausspoints_info, beamproperties, contactparams, sdf, Δt, gausspoints = constants
    iB, X₁, X₂, l₀, Rₑ⁰ = init
    nG, ωG, zG = gausspoints_info
    @unpack K̄ⁱⁿᵗ, Jᵨ, Aᵨ, damping = beamproperties
    
    contactsearch = !(isnothing(sdf) || isnothing(contactparams))

    x₁ =  X₁ + u₁
    x₂ =  X₂ + u₂
    
    lₙ = norm(x₂ - x₁)

    ū = lₙ - l₀

    Rₑ, r¹, r³, η, Gᵀ¹, Gᵀ², Gᵀ³, Gᵀ⁴, D₃ = local_Rₑ_and_aux(x₁, x₂, R₁, R₂, Rₑ⁰[:,2], lₙ)


    R̅₁ = Rₑ' * R₁ * Rₑ⁰
    R̅₂ = Rₑ' * R₂ * Rₑ⁰

    Θ̅₁ = toangle(R̅₁)
    Θ̅₂ = toangle(R̅₂)

    if exact
        Tₛ⁻¹Θ̅₁ = Tₛ⁻¹(Θ̅₁)
        Tₛ⁻¹Θ̅₂ = Tₛ⁻¹(Θ̅₂)
    end


    P¹¹ = -Gᵀ¹ 
    P²¹ = P¹¹
    P¹² = ID3-Gᵀ²
    P²² = -Gᵀ²
    P¹⁴ = -Gᵀ⁴
    P²⁴ = ID3-Gᵀ⁴


    # B̄⁺ = [r; PEᵀ]
    B̄⁺¹ = r¹
    B̄⁺¹¹ = P¹¹ * Rₑ'
    B̄⁺²¹ = B̄⁺¹¹
    B̄⁺¹² = P¹² * Rₑ'
    B̄⁺²² = P²² * Rₑ'
    B̄⁺¹⁴ = P¹⁴ * Rₑ'
    B̄⁺²⁴ = P²⁴ * Rₑ'

    
    # B = B̄B̄⁺
    B¹ = B̄⁺¹
    B¹¹ = exact ?  Tₛ⁻¹Θ̅₁ * B̄⁺¹¹ : B̄⁺¹¹
    B¹² = exact ?  Tₛ⁻¹Θ̅₁ * B̄⁺¹² : B̄⁺¹²
    B¹⁴ = exact ?  Tₛ⁻¹Θ̅₁ * B̄⁺¹⁴ : B̄⁺¹⁴
    B²¹ = exact ?  Tₛ⁻¹Θ̅₂ * B̄⁺²¹ : B̄⁺²¹
    B²² = exact ?  Tₛ⁻¹Θ̅₂ * B̄⁺²² : B̄⁺²²
    B²⁴ = exact ?  Tₛ⁻¹Θ̅₂ * B̄⁺²⁴ : B̄⁺²⁴

    

    K̄ⁱⁿᵗū, K̄ⁱⁿᵗΘ̅, K̄ⁱⁿᵗΘ̅Θ̅ = K̄ⁱⁿᵗ

    # T̄ⁱⁿᵗ = K̄ⁱⁿᵗ D̄
    T̄ⁱⁿᵗū  = K̄ⁱⁿᵗū  * ū
    T̄ⁱⁿᵗΘ̅₁ = K̄ⁱⁿᵗΘ̅  * Θ̅₁ + K̄ⁱⁿᵗΘ̅Θ̅ * Θ̅₂
    T̄ⁱⁿᵗΘ̅₂ = K̄ⁱⁿᵗΘ̅Θ̅ * Θ̅₁ + K̄ⁱⁿᵗΘ̅  * Θ̅₂

    strain_energy = (ū*T̄ⁱⁿᵗū + dot(Θ̅₁, T̄ⁱⁿᵗΘ̅₁) + dot(Θ̅₂, T̄ⁱⁿᵗΘ̅₂))/2

    # Tⁱⁿᵗ = Bᵀ T̄ⁱⁿᵗ
    Tⁱⁿᵗ¹ = B¹'*T̄ⁱⁿᵗū + B¹¹'*T̄ⁱⁿᵗΘ̅₁ + B²¹'*T̄ⁱⁿᵗΘ̅₂
    Tⁱⁿᵗ² =             B¹²'*T̄ⁱⁿᵗΘ̅₁ + B²²'*T̄ⁱⁿᵗΘ̅₂
    Tⁱⁿᵗ³ = -Tⁱⁿᵗ¹
    Tⁱⁿᵗ⁴ =             B¹⁴'*T̄ⁱⁿᵗΘ̅₁ + B²⁴'*T̄ⁱⁿᵗΘ̅₂



    # Force
    Tⁱⁿᵗ = [Tⁱⁿᵗ¹; Tⁱⁿᵗ²; Tⁱⁿᵗ³; Tⁱⁿᵗ⁴]


    # [N̄ M̄⁺₁ M̄⁺₂] = B̄ᵀ T̄ⁱⁿᵗ
    N̄   = T̄ⁱⁿᵗū
    M̄⁺₁ = exact ? Tₛ⁻¹Θ̅₁' * T̄ⁱⁿᵗΘ̅₁  : T̄ⁱⁿᵗΘ̅₁
    M̄⁺₂ = exact ? Tₛ⁻¹Θ̅₂' * T̄ⁱⁿᵗΘ̅₂  : T̄ⁱⁿᵗΘ̅₂


    # Qₛ = Pᵀ [M̄⁺₁ M̄⁺₂]
    Qₛ¹ = P¹¹' * M̄⁺₁ + P²¹' * M̄⁺₂
    Qₛ² = P¹²' * M̄⁺₁ + P²²' * M̄⁺₂
    Qₛ⁴ = P¹⁴' * M̄⁺₁ + P²⁴' * M̄⁺₂
    

    # Q = S(Qₛ)
    Q¹ = skew(Qₛ¹)
    Q² = skew(Qₛ²)
    Q⁴ = skew(Qₛ⁴)

    a = @SVector [0, η*(M̄⁺₁[1] + M̄⁺₂[1])/lₙ + (M̄⁺₁[2] + M̄⁺₂[2])/lₙ, (M̄⁺₁[3] + M̄⁺₂[3])/lₙ]


    # DN̄ (DN̄¹¹ = DN̄³³ = -DN̄¹³ = -DN̄³¹)
    DN̄¹¹ = D₃*N̄

    #QGᵀ
    QGᵀ¹¹ = Q¹*Gᵀ¹
    QGᵀ¹² = Q¹*Gᵀ²
    QGᵀ¹⁴ = Q¹*Gᵀ⁴

    QGᵀ²² = Q²*Gᵀ²
    QGᵀ²⁴ = Q²*Gᵀ⁴

    QGᵀ⁴⁴ = Q⁴*Gᵀ⁴


    # EGa (diagonal)
    # Note: Rₑ*Ga = 0 for Θ indices because Rₑ*GᵀΘ' has only non-zero values in the first column and a = [0 ...]
    EGa¹ = Rₑ*Gᵀ¹'*a

    # EGar (EGar¹¹ = EGar³³ = -EGar³¹ = -EGar¹³)
    EGar¹¹ = EGa¹*r¹

    # Kₘ = DN̄ - EQGᵀEᵀ + EGar
    Kₘ¹¹ = DN̄¹¹ - Rₑ*QGᵀ¹¹*Rₑ' + EGar¹¹
    Kₘ¹² =      - Rₑ*QGᵀ¹²*Rₑ'
    Kₘ¹⁴ =      - Rₑ*QGᵀ¹⁴*Rₑ'
    
    Kₘ²² =      - Rₑ*QGᵀ²²*Rₑ'
    Kₘ²⁴ =      - Rₑ*QGᵀ²⁴*Rₑ'

    Kₘ⁴⁴ =      - Rₑ*QGᵀ⁴⁴*Rₑ'


    # K̃

    if exact

        η₁, μ₁ = compute_η_μ(Θ̅₁)
        η₂, μ₂ = compute_η_μ(Θ̅₂)

        M̄₁ = T̄ⁱⁿᵗΘ̅₁
        M̄₂ = T̄ⁱⁿᵗΘ̅₂

        K̄ₕ₁ = compute_K̄ₕ(Θ̅₁, M̄₁, Tₛ⁻¹Θ̅₁, η₁, μ₁)
        K̄ₕ₂ = compute_K̄ₕ(Θ̅₂, M̄₂, Tₛ⁻¹Θ̅₂, η₂, μ₂)

    end


    K̃¹¹ = exact ?  B̄⁺¹¹' * K̄ₕ₁ * B̄⁺¹¹  +  B̄⁺²¹' * K̄ₕ₂ * B̄⁺²¹  +  Kₘ¹¹   :   Kₘ¹¹ 
    K̃¹² = exact ?  B̄⁺¹¹' * K̄ₕ₁ * B̄⁺¹²  +  B̄⁺²¹' * K̄ₕ₂ * B̄⁺²²  +  Kₘ¹²   :   Kₘ¹² 
    K̃¹⁴ = exact ?  B̄⁺¹¹' * K̄ₕ₁ * B̄⁺¹⁴  +  B̄⁺²¹' * K̄ₕ₂ * B̄⁺²⁴  +  Kₘ¹⁴   :   Kₘ¹⁴ 

    K̃²² = exact ?  B̄⁺¹²' * K̄ₕ₁ * B̄⁺¹²  +  B̄⁺²²' * K̄ₕ₂ * B̄⁺²²  +  Kₘ²²   :   Kₘ²² 
    K̃²⁴ = exact ?  B̄⁺¹²' * K̄ₕ₁ * B̄⁺¹⁴  +  B̄⁺²²' * K̄ₕ₂ * B̄⁺²⁴  +  Kₘ²⁴   :   Kₘ²⁴ 

    K̃⁴⁴ = exact ?  B̄⁺¹⁴' * K̄ₕ₁ * B̄⁺¹⁴  +  B̄⁺²⁴' * K̄ₕ₂ * B̄⁺²⁴  +  Kₘ⁴⁴   :   Kₘ⁴⁴ 


    # BᵀKⁱⁿᵗB
    
    BᵀKⁱⁿᵗB¹¹ = B¹' * K̄ⁱⁿᵗū * B¹      +      B¹¹' * (K̄ⁱⁿᵗΘ̅ * B¹¹  +  K̄ⁱⁿᵗΘ̅Θ̅ * B²¹)  +  B²¹' * (K̄ⁱⁿᵗΘ̅Θ̅ * B¹¹  +  K̄ⁱⁿᵗΘ̅ * B²¹)    
    BᵀKⁱⁿᵗB¹² =                              B¹¹' * (K̄ⁱⁿᵗΘ̅ * B¹²  +  K̄ⁱⁿᵗΘ̅Θ̅ * B²²)  +  B²¹' * (K̄ⁱⁿᵗΘ̅Θ̅ * B¹²  +  K̄ⁱⁿᵗΘ̅ * B²²)     
    BᵀKⁱⁿᵗB¹⁴ =                              B¹¹' * (K̄ⁱⁿᵗΘ̅ * B¹⁴  +  K̄ⁱⁿᵗΘ̅Θ̅ * B²⁴)  +  B²¹' * (K̄ⁱⁿᵗΘ̅Θ̅ * B¹⁴  +  K̄ⁱⁿᵗΘ̅ * B²⁴)    
        
    BᵀKⁱⁿᵗB²² =                              B¹²' * (K̄ⁱⁿᵗΘ̅ * B¹²  +  K̄ⁱⁿᵗΘ̅Θ̅ * B²²)  +  B²²' * (K̄ⁱⁿᵗΘ̅Θ̅ * B¹²  +  K̄ⁱⁿᵗΘ̅ * B²²)      
    BᵀKⁱⁿᵗB²⁴ =                              B¹²' * (K̄ⁱⁿᵗΘ̅ * B¹⁴  +  K̄ⁱⁿᵗΘ̅Θ̅ * B²⁴)  +  B²²' * (K̄ⁱⁿᵗΘ̅Θ̅ * B¹⁴  +  K̄ⁱⁿᵗΘ̅ * B²⁴)      
        
    BᵀKⁱⁿᵗB⁴⁴ =                              B¹⁴' * (K̄ⁱⁿᵗΘ̅ * B¹⁴  +  K̄ⁱⁿᵗΘ̅Θ̅ * B²⁴)  +  B²⁴' * (K̄ⁱⁿᵗΘ̅Θ̅ * B¹⁴  +  K̄ⁱⁿᵗΘ̅ * B²⁴)    

    Kⁱⁿᵗ¹¹ =  BᵀKⁱⁿᵗB¹¹   +    K̃¹¹
    Kⁱⁿᵗ¹² =  BᵀKⁱⁿᵗB¹²   +    K̃¹²
    Kⁱⁿᵗ¹³ =  -Kⁱⁿᵗ¹¹
    Kⁱⁿᵗ¹⁴ =  BᵀKⁱⁿᵗB¹⁴   +    K̃¹⁴

    Kⁱⁿᵗ²² =  BᵀKⁱⁿᵗB²²   +    K̃²²
    Kⁱⁿᵗ²³ =  -Kⁱⁿᵗ¹²'
    Kⁱⁿᵗ²⁴ =  BᵀKⁱⁿᵗB²⁴   +    K̃²⁴
    
    Kⁱⁿᵗ³³ =  Kⁱⁿᵗ¹¹
    Kⁱⁿᵗ³⁴ =  -Kⁱⁿᵗ¹⁴

    Kⁱⁿᵗ⁴⁴ =  BᵀKⁱⁿᵗB⁴⁴   +    K̃⁴⁴



    Kⁱⁿᵗ = hcat(vcat(Kⁱⁿᵗ¹¹, Kⁱⁿᵗ¹²', Kⁱⁿᵗ¹³', Kⁱⁿᵗ¹⁴'), vcat(Kⁱⁿᵗ¹², Kⁱⁿᵗ²², Kⁱⁿᵗ²³', Kⁱⁿᵗ²⁴'), vcat(Kⁱⁿᵗ¹³, Kⁱⁿᵗ²³, Kⁱⁿᵗ³³, Kⁱⁿᵗ³⁴'), vcat(Kⁱⁿᵗ¹⁴, Kⁱⁿᵗ²⁴, Kⁱⁿᵗ³⁴, Kⁱⁿᵗ⁴⁴))






    kinetic_energy = zero(T)

    Tᵏ¹ = zeros(Vec3{T})
    Tᵏ² = zeros(Vec3{T})
    Tᵏ³ = zeros(Vec3{T})
    Tᵏ⁴ = zeros(Vec3{T})


    M¹¹ = zeros(Mat33{T})
    M¹² = zeros(Mat33{T})
    M¹³ = zeros(Mat33{T})
    M¹⁴ = zeros(Mat33{T})
    M²¹ = zeros(Mat33{T})
    M²² = zeros(Mat33{T})
    M²³ = zeros(Mat33{T})
    M²⁴ = zeros(Mat33{T})
    M³¹ = zeros(Mat33{T})
    M³² = zeros(Mat33{T})
    M³³ = zeros(Mat33{T})
    M³⁴ = zeros(Mat33{T})
    M⁴¹ = zeros(Mat33{T})
    M⁴² = zeros(Mat33{T})
    M⁴³ = zeros(Mat33{T})
    M⁴⁴ = zeros(Mat33{T})

    Cᵏ¹¹ = zeros(Mat33{T})
    Cᵏ¹² = zeros(Mat33{T})
    Cᵏ¹³ = zeros(Mat33{T})
    Cᵏ¹⁴ = zeros(Mat33{T})
    Cᵏ²¹ = zeros(Mat33{T})
    Cᵏ²² = zeros(Mat33{T})
    Cᵏ²³ = zeros(Mat33{T})
    Cᵏ²⁴ = zeros(Mat33{T})
    Cᵏ³¹ = zeros(Mat33{T})
    Cᵏ³² = zeros(Mat33{T})
    Cᵏ³³ = zeros(Mat33{T})
    Cᵏ³⁴ = zeros(Mat33{T})
    Cᵏ⁴¹ = zeros(Mat33{T})
    Cᵏ⁴² = zeros(Mat33{T})
    Cᵏ⁴³ = zeros(Mat33{T})
    Cᵏ⁴⁴ = zeros(Mat33{T})




    contact_energy = zero(T)
        

    Tᶜ¹ = zeros(Vec3{T})
    Tᶜ² = zeros(Vec3{T})
    Tᶜ³ = zeros(Vec3{T})
    Tᶜ⁴ = zeros(Vec3{T})



    Kᶜ¹¹ = zeros(Mat33{T})
    Kᶜ¹² = zeros(Mat33{T})
    Kᶜ¹³ = zeros(Mat33{T})
    Kᶜ¹⁴ = zeros(Mat33{T})
    Kᶜ²¹ = zeros(Mat33{T})
    Kᶜ²² = zeros(Mat33{T})
    Kᶜ²³ = zeros(Mat33{T})
    Kᶜ²⁴ = zeros(Mat33{T})
    Kᶜ³¹ = zeros(Mat33{T})
    Kᶜ³² = zeros(Mat33{T})
    Kᶜ³³ = zeros(Mat33{T})
    Kᶜ³⁴ = zeros(Mat33{T})
    Kᶜ⁴¹ = zeros(Mat33{T})
    Kᶜ⁴² = zeros(Mat33{T})
    Kᶜ⁴³ = zeros(Mat33{T})
    Kᶜ⁴⁴ = zeros(Mat33{T})

    Cᶜ¹¹ = zeros(Mat33{T})
    Cᶜ¹² = zeros(Mat33{T})
    Cᶜ¹³ = zeros(Mat33{T})
    Cᶜ¹⁴ = zeros(Mat33{T})
    Cᶜ²¹ = zeros(Mat33{T})
    Cᶜ²² = zeros(Mat33{T})
    Cᶜ²³ = zeros(Mat33{T})
    Cᶜ²⁴ = zeros(Mat33{T})
    Cᶜ³¹ = zeros(Mat33{T})
    Cᶜ³² = zeros(Mat33{T})
    Cᶜ³³ = zeros(Mat33{T})
    Cᶜ³⁴ = zeros(Mat33{T})
    Cᶜ⁴¹ = zeros(Mat33{T})
    Cᶜ⁴² = zeros(Mat33{T})
    Cᶜ⁴³ = zeros(Mat33{T})
    Cᶜ⁴⁴ = zeros(Mat33{T})

        
    if isdynamic
        
        U̇₁ = Rₑ' * u̇₁
        U̇₂ = Rₑ' * u̇₂
        Ẇ₁ = Rₑ' * ẇ₁
        Ẇ₂ = Rₑ' * ẇ₂
        
        Ü₁ = Rₑ' * ü₁
        Ü₂ = Rₑ' * ü₂
        Ẅ₁ = Rₑ' * ẅ₁
        Ẅ₂ = Rₑ' * ẅ₂

        SU̇₁ = skew(U̇₁)
        SU̇₂ = skew(U̇₂)
        SẆ₁ = skew(Ẇ₁)
        SẆ₂ = skew(Ẇ₂)

        Ẇᵉ = Gᵀ¹ * U̇₁ + Gᵀ² * Ẇ₁ + Gᵀ³ * U̇₂ + Gᵀ⁴ * Ẇ₂
        SẆᵉ = skew(Ẇᵉ)

        rḋ = dot(r¹, u̇₁) + dot(r³, u̇₂)


        Gᵀ¹Rₑᵀ = Gᵀ¹ * Rₑ'
        Gᵀ²Rₑᵀ = Gᵀ² * Rₑ'
        Gᵀ⁴Rₑᵀ = Gᵀ⁴ * Rₑ'

        RₑG¹ = Rₑ * Gᵀ¹'
        RₑG² = Rₑ * Gᵀ²'
        RₑG⁴ = Rₑ * Gᵀ⁴'


        # cycle among the Gauss positions
        for iG in 1:nG

            zᴳ = zG[iG]
            ωᴳ = ωG[iG]

            ξ = l₀*(zᴳ+1)/2

            # Shape functions
            N₁ = 1-ξ/l₀
            N₂ = 1-N₁
            N₃ = ξ*(1-ξ/l₀)^2
            N₄ = -(1-ξ/l₀)*((ξ^2)/l₀)
            N₅ = (1-3*ξ/l₀)*(1-ξ/l₀)
            N₆ = (3*ξ/l₀-2)*(ξ/l₀)
            N₇ = N₃+N₄
            N₈ = N₅+N₆-1


            uᵗ = @SVector [0, N₃*Θ̅₁[3] + N₄*Θ̅₂[3], -N₃*Θ̅₁[2] + -N₄*Θ̅₂[2]]
            Θ̄  = @SVector [N₁*Θ̅₁[1] + N₂*Θ̅₂[1], N₅*Θ̅₁[2] + N₆*Θ̅₂[2], N₅*Θ̅₁[3] + N₆*Θ̅₂[3]]

            Suᵗ = skew(uᵗ)
            SΘ̄ = skew(Θ̄)

            R̄ = ID3 + SΘ̄

            Īᵨ = R̄*Jᵨ*R̄'

            N₇lₙ = N₇/lₙ
            N₇lₙ² = N₇lₙ/lₙ
            N₈lₙ = N₈/lₙ
            N₈lₙ² = N₈lₙ/lₙ

            P₁P¹ = @SMatrix [0 0 0; 0 N₇lₙ 0;0 0 N₇lₙ]
            P₁P² = @SMatrix [0 0 0; 0 0 N₃;0 -N₃ 0]
            P₁P³ = -P₁P¹
            P₁P⁴ = @SMatrix [0 0 0; 0 0 N₄;0 -N₄ 0]

            H₁¹ = N₁*ID3 + P₁P¹ - Suᵗ*Gᵀ¹
            H₁² =          P₁P² - Suᵗ*Gᵀ²
            H₁³ = N₂*ID3 + P₁P³ - Suᵗ*Gᵀ³
            H₁⁴ =          P₁P⁴ - Suᵗ*Gᵀ⁴

            H₂¹ = @SMatrix [0 0 0; 0  0 -N₈lₙ;0 N₈lₙ 0]
            H₂² = Diagonal(@SVector [N₁, N₅, N₅])
            H₂³ = -H₂¹
            H₂⁴ = Diagonal(@SVector [N₂, N₆, N₆])


            u̇ᵗ =  P₁P¹ * U̇₁ +  P₁P² * Ẇ₁ + P₁P³ * U̇₂ + P₁P⁴ * Ẇ₂

            Su̇ᵗ = skew(u̇ᵗ)
            
            N₇rḋ = N₇lₙ² * rḋ
            Ḣ₁¹ = Diagonal(@SVector [0, -N₇rḋ, -N₇rḋ]) - Su̇ᵗ * Gᵀ¹
            Ḣ₁² =                                      - Su̇ᵗ * Gᵀ²
            Ḣ₁⁴ =                                      - Su̇ᵗ * Gᵀ⁴

            N₈rḋ = N₈lₙ² * rḋ
            Ḣ₂¹ = @SMatrix [0 0 0; 0 0 N₈rḋ; 0 -N₈rḋ 0]

            h₁ = H₁¹ * U̇₁ + H₁² * Ẇ₁ + H₁³ * U̇₂ + H₁⁴ * Ẇ₂
            h₂ = H₂¹ * U̇₁ + H₂² * Ẇ₁ + H₂³ * U̇₂ + H₂⁴ * Ẇ₂
            Sh₁ = skew(h₁)
            Sh₂ = skew(h₂)

            C₁¹ = SẆᵉ * H₁¹ + Ḣ₁¹ - H₁¹ * SẆᵉ
            C₁² = SẆᵉ * H₁² + Ḣ₁² - H₁² * SẆᵉ
            C₁³ = -C₁¹
            C₁⁴ = SẆᵉ * H₁⁴ + Ḣ₁⁴ - H₁⁴ * SẆᵉ

            C₂¹ = SẆᵉ * H₂¹ + Ḣ₂¹ - H₂¹ * SẆᵉ
            C₂² = SẆᵉ * H₂²       - H₂² * SẆᵉ
            C₂³ = -C₂¹
            C₂⁴ = SẆᵉ * H₂⁴       - H₂⁴ * SẆᵉ

            H₁F₁ = H₁¹ * SU̇₁ + H₁² * SẆ₁ + H₁³ * SU̇₂ + H₁⁴ * SẆ₂
            H₂F₁ = H₂¹ * SU̇₁ + H₂² * SẆ₁ + H₂³ * SU̇₂ + H₂⁴ * SẆ₂

            A₁ḊrE¹ = @SMatrix [0 0 0; U̇₁[2]-U̇₂[2] 0 0; U̇₁[3]-U̇₂[3] 0 0]
            C₃¹ = -Sh₁*Gᵀ¹ + N₇lₙ²*A₁ḊrE¹ + SẆᵉ*P₁P¹ + H₁F₁*Gᵀ¹
            C₃² = -Sh₁*Gᵀ² +                  SẆᵉ*P₁P² + H₁F₁*Gᵀ²
            C₃⁴ = -Sh₁*Gᵀ⁴ +                  SẆᵉ*P₁P⁴ + H₁F₁*Gᵀ⁴

            A₂ḊrE¹ = @SMatrix [0 0 0; -U̇₁[3]+U̇₂[3] 0 0; U̇₁[2]-U̇₂[2] 0 0]
            C₄¹ = -Sh₂*Gᵀ¹ + N₈lₙ²*A₂ḊrE¹ + H₂F₁*Gᵀ¹
            C₄² = -Sh₂*Gᵀ²                  + H₂F₁*Gᵀ²
            C₄⁴ = -Sh₂*Gᵀ⁴                  + H₂F₁*Gᵀ⁴

            u̇₀ = Rₑ * h₁

            H₁Eᵀd̈ = H₁¹ * Ü₁ + H₁² * Ẅ₁ + H₁³ * Ü₂ + H₁⁴ * Ẅ₂
            C₁Eᵀḋ = C₁¹ * U̇₁ + C₁² * Ẇ₁ + C₁³ * U̇₂ + C₁⁴ * Ẇ₂
            Rₑᵀü₀ = H₁Eᵀd̈ + C₁Eᵀḋ

            Ẇ₀ = h₂
            ẇ₀ = Rₑ * Ẇ₀

            H₂Eᵀd̈ = H₂¹ * Ü₁ + H₂² * Ẅ₁ + H₂³ * Ü₂ + H₂⁴ * Ẅ₂
            C₂Eᵀḋ = C₂¹ * U̇₁ + C₂² * Ẇ₁ + C₂³ * U̇₂ + C₂⁴ * Ẇ₂
            Rₑᵀẅ₀ = H₂Eᵀd̈ + C₂Eᵀḋ

            SẆ₀ = skew(Ẇ₀)
            ĪᵨRₑᵀẅ₀ = Īᵨ*Rₑᵀẅ₀
            SẆ₀Īᵨ = SẆ₀*Īᵨ
            SẆ₀ĪᵨẆ₀ = SẆ₀Īᵨ*Ẇ₀
            ĪᵨRₑᵀẅ₀SẆ₀ĪᵨẆ₀ = ĪᵨRₑᵀẅ₀ + SẆ₀ĪᵨẆ₀
            AᵨH₁¹ᵀ = Aᵨ*H₁¹'
            AᵨH₁²ᵀ = Aᵨ*H₁²'
            AᵨH₁³ᵀ = Aᵨ*H₁³'
            AᵨH₁⁴ᵀ = Aᵨ*H₁⁴'

            Tᵏ¹G = ωᴳ * (AᵨH₁¹ᵀ*Rₑᵀü₀ + H₂¹'*ĪᵨRₑᵀẅ₀SẆ₀ĪᵨẆ₀)
            Tᵏ¹ += Tᵏ¹G
            Tᵏ² += ωᴳ * (AᵨH₁²ᵀ*Rₑᵀü₀ + H₂²'*ĪᵨRₑᵀẅ₀SẆ₀ĪᵨẆ₀)
            Tᵏ³ += -Tᵏ¹G + ωᴳ * Aᵨ * Rₑᵀü₀
            Tᵏ⁴ += ωᴳ * (AᵨH₁⁴ᵀ*Rₑᵀü₀ + H₂⁴'*ĪᵨRₑᵀẅ₀SẆ₀ĪᵨẆ₀)


            if damping>0
                Tᵈ¹G = ωᴳ * damping*(AᵨH₁¹ᵀ*h₁ + H₂¹'*Īᵨ*h₂)
                Tᵏ¹ += Tᵈ¹G
                Tᵏ² += ωᴳ * damping*(AᵨH₁²ᵀ*h₁ + H₂²'*Īᵨ*h₂)
                Tᵏ³ += -Tᵈ¹G + ωᴳ * Aᵨ * damping * h₁
                Tᵏ⁴ += ωᴳ * damping*(AᵨH₁⁴ᵀ*h₁ + H₂⁴'*Īᵨ*h₂)
            end


            M¹¹G = ωᴳ * (AᵨH₁¹ᵀ*H₁¹ + H₂¹'*Īᵨ*H₂¹)
            M¹²G = ωᴳ * (AᵨH₁¹ᵀ*H₁² + H₂¹'*Īᵨ*H₂²)
            M¹⁴G = ωᴳ * (AᵨH₁¹ᵀ*H₁⁴ + H₂¹'*Īᵨ*H₂⁴)

            M²²G = ωᴳ * (AᵨH₁²ᵀ*H₁² + H₂²'*Īᵨ*H₂²)
            M²⁴G = ωᴳ * (AᵨH₁²ᵀ*H₁⁴ + H₂²'*Īᵨ*H₂⁴)
            
            M⁴⁴G = ωᴳ * (AᵨH₁⁴ᵀ*H₁⁴ + H₂⁴'*Īᵨ*H₂⁴)


            M¹¹ += M¹¹G
            M¹² += M¹²G
            M¹³ += -M¹¹G + ωᴳ * AᵨH₁¹ᵀ
            M¹⁴ += M¹⁴G

            M²² += M²²G
            M²³ += -M¹²G' + ωᴳ * AᵨH₁²ᵀ
            M²⁴ += M²⁴G

            M³³ += M¹¹G + ωᴳ * Aᵨ * (ID3 - H₁¹ - H₁¹')
            M³⁴ += -M¹⁴G + ωᴳ * Aᵨ*H₁⁴
            
            M⁴⁴ += M⁴⁴G


            SẆ₀ĪᵨmSĪᵨẆ₀ = SẆ₀Īᵨ - skew(Īᵨ*Ẇ₀)

            AᵨC₁¹C₃¹ = Aᵨ*(C₁¹ + C₃¹)
            AᵨC₁¹C₃² = Aᵨ*(C₁² + C₃²)
            AᵨC₁¹C₃⁴ = Aᵨ*(C₁⁴ + C₃⁴)

            ĪᵨC₂¹C₄¹SẆ₀ĪᵨmSĪᵨẆ₀H₂¹ = Īᵨ*(C₂¹ + C₄¹) + SẆ₀ĪᵨmSĪᵨẆ₀ * H₂¹
            ĪᵨC₂²C₄¹SẆ₀ĪᵨmSĪᵨẆ₀H₂² = Īᵨ*(C₂² + C₄²) + SẆ₀ĪᵨmSĪᵨẆ₀ * H₂²
            ĪᵨC₂⁴C₄¹SẆ₀ĪᵨmSĪᵨẆ₀H₂⁴ = Īᵨ*(C₂⁴ + C₄⁴) + SẆ₀ĪᵨmSĪᵨẆ₀ * H₂⁴


            Cᵏ¹¹G = ωᴳ * ( H₁¹'*AᵨC₁¹C₃¹ + H₂¹'* ĪᵨC₂¹C₄¹SẆ₀ĪᵨmSĪᵨẆ₀H₂¹ ) 
            Cᵏ¹²G = ωᴳ * ( H₁¹'*AᵨC₁¹C₃² + H₂¹'* ĪᵨC₂²C₄¹SẆ₀ĪᵨmSĪᵨẆ₀H₂² ) 
            Cᵏ¹⁴G = ωᴳ * ( H₁¹'*AᵨC₁¹C₃⁴ + H₂¹'* ĪᵨC₂⁴C₄¹SẆ₀ĪᵨmSĪᵨẆ₀H₂⁴ ) 
 
            Cᵏ²¹G = ωᴳ * ( H₁²'*AᵨC₁¹C₃¹ + H₂²'* ĪᵨC₂¹C₄¹SẆ₀ĪᵨmSĪᵨẆ₀H₂¹ ) 
            Cᵏ²²G = ωᴳ * ( H₁²'*AᵨC₁¹C₃² + H₂²'* ĪᵨC₂²C₄¹SẆ₀ĪᵨmSĪᵨẆ₀H₂² ) 
            Cᵏ²⁴G = ωᴳ * ( H₁²'*AᵨC₁¹C₃⁴ + H₂²'* ĪᵨC₂⁴C₄¹SẆ₀ĪᵨmSĪᵨẆ₀H₂⁴ ) 
 
            Cᵏ⁴¹G = ωᴳ * ( H₁⁴'*AᵨC₁¹C₃¹ + H₂⁴'* ĪᵨC₂¹C₄¹SẆ₀ĪᵨmSĪᵨẆ₀H₂¹ ) 
            Cᵏ⁴²G = ωᴳ * ( H₁⁴'*AᵨC₁¹C₃² + H₂⁴'* ĪᵨC₂²C₄¹SẆ₀ĪᵨmSĪᵨẆ₀H₂² ) 
            Cᵏ⁴⁴G = ωᴳ * ( H₁⁴'*AᵨC₁¹C₃⁴ + H₂⁴'* ĪᵨC₂⁴C₄¹SẆ₀ĪᵨmSĪᵨẆ₀H₂⁴ ) 


            Cᵏ¹¹ += Cᵏ¹¹G
            Cᵏ¹² += Cᵏ¹²G
            Cᵏ¹⁴ += Cᵏ¹⁴G

            Cᵏ²¹ += Cᵏ²¹G
            Cᵏ²² += Cᵏ²²G
            Cᵏ²⁴ += Cᵏ²⁴G

            Cᵏ³¹ += -Cᵏ¹¹G + ωᴳ * AᵨC₁¹C₃¹
            Cᵏ³² += -Cᵏ¹²G + ωᴳ * AᵨC₁¹C₃²
            Cᵏ³⁴ += -Cᵏ¹⁴G + ωᴳ * AᵨC₁¹C₃⁴

            Cᵏ⁴¹ += Cᵏ⁴¹G
            Cᵏ⁴² += Cᵏ⁴²G
            Cᵏ⁴⁴ += Cᵏ⁴⁴G
            
            # kinetic energy
            Īᵨᵍ = Rₑ*Īᵨ*Rₑ'
            kinetic_energy += ωᴳ/2 * (Aᵨ*u̇₀'*u̇₀ + ẇ₀'*Īᵨᵍ*ẇ₀)

            if contactsearch

                xᴳ = N₁*x₁ + N₂*x₂ + Rₑ*uᵗ
                gausspoints.pos[(iB-1)*3 + iG] = xᴳ

                ḡₙ = sdf.r/4
            
                if incontact(xᴳ, sdf, ḡₙ)

                    gₙ, ∂gₙ∂x, ∂²gₙ∂x² = contact_gap(xᴳ, sdf)

                    @unpack kₙ, ηₙ, μ, εᵗ, kₜ, ηₜ, u̇ₛ = contactparams

                    pₙ, p′ₙ, Πₑ = regularize_gₙ(gₙ, ḡₙ)
                    ηₙ, η′ₙ = smoothstep(ηₙ, gₙ, ḡₙ)

                    u̇ₙ_mag = dot(u̇₀, ∂gₙ∂x)
                    u̇ₙ = u̇ₙ_mag*∂gₙ∂x
                    u̇ₜ = u̇₀ - u̇ₙ
                    u̇ₜ² = dot(u̇ₜ, u̇ₜ)
        
                    contact_energy += ωᴳ*kₙ*Πₑ
                    
                    𝓯ⁿ = kₙ * pₙ * ∂gₙ∂x - ηₙ * u̇ₙ

                    δₜ =  gausspoints.δₜ[(iB-1)*3 + iG]
                    δₜ = δₜ + Δt*norm(u̇ₜ)

                    μʳᵉᵍ = μ/sqrt(u̇ₜ²+εᵗ)
                    
                    if  norm(u̇ₜ) <= u̇ₛ
                        # if  gausspoints.status[(iB-1)*3 + iG] == 1
                        #     δₜ = μʳᵉᵍ*pₙ/kₜ - ηₜ/kₜ*sqrt(u̇ₜ²+εᵗ)
                        # end 
                        gausspoints.status[(iB-1)*3 + iG] = 2
                        # break
                    elseif norm(u̇ₜ) > u̇ₛ
                        gausspoints.status[(iB-1)*3 + iG] = 1
                        # break
                    end
                    # elseif gausspoints.status[(iB-1)*3 + iG] == 1 && norm(u̇ₜ) <= u̇ₛ

                    #     gausspoints.status[(iB-1)*3 + iG] = 2
                    #     # δₜ = μʳᵉᵍ*pₙ/kₜ - ηₜ/kₜ*sqrt(u̇ₜ²+εᵗ)
                    #     break
                        
                    # elseif gausspoints.status[(iB-1)*3 + iG] == 2
                        
                    #     𝓯ᵗ = - kₜ * δₜ * u̇ₜ/sqrt(u̇ₜ²+εᵗ) - ηₜ*u̇ₜ

                    #     if norm(𝓯ᵗ) > μstick*pₙ && norm(u̇ₜ) > u̇ₛ
                    #         gausspoints.status[(iB-1)*3 + iG] = 1
                    #     end 
                    #     break

                    # end
                    gausspoints.δₜ[(iB-1)*3 + iG] =  δₜ

                    𝓯ᵗ = zeros(Vec3)
                    if gausspoints.status[(iB-1)*3 + iG] == 1
                        𝓯ᵗ = - kₙ * pₙ * μʳᵉᵍ * u̇ₜ   
                    elseif gausspoints.status[(iB-1)*3 + iG] == 2
                        𝓯ᵗ = - kₜ * δₜ * u̇ₜ/sqrt(u̇ₜ²+εᵗ) - ηₜ*u̇ₜ
                    end

              
                    𝓯ᶜ = 𝓯ⁿ + 𝓯ᵗ

                    𝓕ᶜ = Rₑ' * 𝓯ᶜ

                    RₑH₁¹Rₑᵀ = Rₑ * H₁¹ * Rₑ'
                    RₑH₁²Rₑᵀ = Rₑ * H₁² * Rₑ'
                    RₑH₁³Rₑᵀ = Rₑ * H₁³ * Rₑ'
                    RₑH₁⁴Rₑᵀ = Rₑ * H₁⁴ * Rₑ'
        
                    Tᶜ¹ += ωᴳ * (RₑH₁¹Rₑᵀ' * 𝓯ᶜ)
                    Tᶜ² += ωᴳ * (RₑH₁²Rₑᵀ' * 𝓯ᶜ)
                    Tᶜ³ += ωᴳ * (RₑH₁³Rₑᵀ' * 𝓯ᶜ)
                    Tᶜ⁴ += ωᴳ * (RₑH₁⁴Rₑᵀ' * 𝓯ᶜ)
        
                    ŜH₁ᵀ𝓕ᶜ¹ = skew(H₁¹' * 𝓕ᶜ)
                    ŜH₁ᵀ𝓕ᶜ² = skew(H₁²' * 𝓕ᶜ)
                    ŜH₁ᵀ𝓕ᶜ³ = skew(H₁³' * 𝓕ᶜ)
                    ŜH₁ᵀ𝓕ᶜ⁴ = skew(H₁⁴' * 𝓕ᶜ)
        
                    RₑŜH₁ᵀ𝓕ᶜ¹ = Rₑ * ŜH₁ᵀ𝓕ᶜ¹
                    t₁¹¹ = -RₑŜH₁ᵀ𝓕ᶜ¹ * Gᵀ¹Rₑᵀ
                    t₁¹² = -RₑŜH₁ᵀ𝓕ᶜ¹ * Gᵀ²Rₑᵀ
                    t₁¹³ = -t₁¹¹
                    t₁¹⁴ = -RₑŜH₁ᵀ𝓕ᶜ¹ * Gᵀ⁴Rₑᵀ
        
                    RₑŜH₁ᵀ𝓕ᶜ² = Rₑ * ŜH₁ᵀ𝓕ᶜ²
                    t₁²¹ = -RₑŜH₁ᵀ𝓕ᶜ² * Gᵀ¹Rₑᵀ
                    t₁²² = -RₑŜH₁ᵀ𝓕ᶜ² * Gᵀ²Rₑᵀ
                    t₁²³ = -t₁²¹
                    t₁²⁴ = -RₑŜH₁ᵀ𝓕ᶜ² * Gᵀ⁴Rₑᵀ
        
                    RₑŜH₁ᵀ𝓕ᶜ³ = Rₑ * ŜH₁ᵀ𝓕ᶜ³
                    t₁³¹ = -RₑŜH₁ᵀ𝓕ᶜ³ * Gᵀ¹Rₑᵀ
                    t₁³² = -RₑŜH₁ᵀ𝓕ᶜ³ * Gᵀ²Rₑᵀ
                    t₁³³ = -t₁³¹
                    t₁³⁴ = -RₑŜH₁ᵀ𝓕ᶜ³ * Gᵀ⁴Rₑᵀ
        
                    RₑŜH₁ᵀ𝓕ᶜ⁴ = Rₑ * ŜH₁ᵀ𝓕ᶜ⁴
                    t₁⁴¹ = -RₑŜH₁ᵀ𝓕ᶜ⁴ * Gᵀ¹Rₑᵀ
                    t₁⁴² = -RₑŜH₁ᵀ𝓕ᶜ⁴ * Gᵀ²Rₑᵀ
                    t₁⁴³ = -t₁⁴¹
                    t₁⁴⁴ = -RₑŜH₁ᵀ𝓕ᶜ⁴ * Gᵀ⁴Rₑᵀ
        
        
                    v₁ = r³
                    A₁ᵀ𝓕ᶜr₁₁ = @SMatrix[0 0 0; 𝓕ᶜ[2]*v₁[1] 𝓕ᶜ[2]*v₁[2] 𝓕ᶜ[2]*v₁[3]; 𝓕ᶜ[3]*v₁[1] 𝓕ᶜ[3]*v₁[2] 𝓕ᶜ[3]*v₁[3]]
        
                    S𝓕ᶜ = skew(𝓕ᶜ)
        
                    S𝓕ᶜP₁P¹Rₑᵀ = S𝓕ᶜ * P₁P¹ * Rₑ'
                    S𝓕ᶜP₁P²Rₑᵀ = S𝓕ᶜ * P₁P² * Rₑ'
                    S𝓕ᶜP₁P⁴Rₑᵀ = S𝓕ᶜ * P₁P⁴ * Rₑ'
                    
                    t₂¹¹ = N₇lₙ² * Rₑ * A₁ᵀ𝓕ᶜr₁₁ - RₑG¹ * S𝓕ᶜP₁P¹Rₑᵀ
                    t₂¹² =                       - RₑG¹ * S𝓕ᶜP₁P²Rₑᵀ
                    t₂¹³ = -t₂¹¹
                    t₂¹⁴ =                       - RₑG¹ * S𝓕ᶜP₁P⁴Rₑᵀ
        
                    t₂²¹ =                       - RₑG² * S𝓕ᶜP₁P¹Rₑᵀ
                    t₂²² =                       - RₑG² * S𝓕ᶜP₁P²Rₑᵀ
                    t₂²³ = -t₂²¹
                    t₂²⁴ =                       - RₑG² * S𝓕ᶜP₁P⁴Rₑᵀ
        
                    t₂³¹ = t₂¹³
                    t₂³² = -t₂¹²
                    t₂³³ = t₂¹¹
                    t₂³⁴ = -t₂¹⁴
        
                    t₂⁴¹ =                       - RₑG⁴ * S𝓕ᶜP₁P¹Rₑᵀ
                    t₂⁴² =                       - RₑG⁴ * S𝓕ᶜP₁P²Rₑᵀ
                    t₂⁴³ = -t₂⁴¹
                    t₂⁴⁴ =                       - RₑG⁴ * S𝓕ᶜP₁P⁴Rₑᵀ
        
                    RₑH₁ᵀ¹S𝓕ᶜ = Rₑ * H₁¹' * S𝓕ᶜ
                    RₑH₁ᵀ²S𝓕ᶜ = Rₑ * H₁²' * S𝓕ᶜ
                    RₑH₁ᵀ³S𝓕ᶜ = Rₑ * H₁³' * S𝓕ᶜ
                    RₑH₁ᵀ⁴S𝓕ᶜ = Rₑ * H₁⁴' * S𝓕ᶜ
        
                    t₃¹¹ = RₑH₁ᵀ¹S𝓕ᶜ * Gᵀ¹Rₑᵀ
                    t₃¹² = RₑH₁ᵀ¹S𝓕ᶜ * Gᵀ²Rₑᵀ
                    t₃¹³ = -t₃¹¹
                    t₃¹⁴ = RₑH₁ᵀ¹S𝓕ᶜ * Gᵀ⁴Rₑᵀ
        
                    t₃²¹ = RₑH₁ᵀ²S𝓕ᶜ * Gᵀ¹Rₑᵀ
                    t₃²² = RₑH₁ᵀ²S𝓕ᶜ * Gᵀ²Rₑᵀ
                    t₃²³ = -t₃²¹
                    t₃²⁴ = RₑH₁ᵀ²S𝓕ᶜ * Gᵀ⁴Rₑᵀ
        
                    t₃³¹ = RₑH₁ᵀ³S𝓕ᶜ * Gᵀ¹Rₑᵀ
                    t₃³² = RₑH₁ᵀ³S𝓕ᶜ * Gᵀ²Rₑᵀ
                    t₃³³ = -t₃³¹
                    t₃³⁴ = RₑH₁ᵀ³S𝓕ᶜ * Gᵀ⁴Rₑᵀ
        
                    t₃⁴¹ = RₑH₁ᵀ⁴S𝓕ᶜ * Gᵀ¹Rₑᵀ
                    t₃⁴² = RₑH₁ᵀ⁴S𝓕ᶜ * Gᵀ²Rₑᵀ
                    t₃⁴³ = -t₃⁴¹
                    t₃⁴⁴ = RₑH₁ᵀ⁴S𝓕ᶜ * Gᵀ⁴Rₑᵀ

        
                    aux = dot(∂gₙ∂x, u̇₀) * ∂²gₙ∂x² + ∂gₙ∂x * (∂²gₙ∂x² * u̇₀)'
        
                    𝓐₁¹ =  aux * RₑH₁¹Rₑᵀ
                    𝓐₁² =  aux * RₑH₁²Rₑᵀ
                    𝓐₁³ =  aux * RₑH₁³Rₑᵀ
                    𝓐₁⁴ =  aux * RₑH₁⁴Rₑᵀ

        

                    RₑSh₁ = Rₑ * Sh₁
                    𝓐₂¹ = - RₑSh₁ * Gᵀ¹Rₑᵀ
                    𝓐₂² = - RₑSh₁ * Gᵀ²Rₑᵀ
                    𝓐₂⁴ = - RₑSh₁ * Gᵀ⁴Rₑᵀ

                    A₁Ḋr¹ = @SMatrix [0 0 0; v₁[1]*(U̇₁[2]-U̇₂[2]) v₁[2]*(U̇₁[2]-U̇₂[2]) v₁[3]*(U̇₁[2]-U̇₂[2]);v₁[1]*(U̇₁[3]-U̇₂[3]) v₁[2]*(U̇₁[3]-U̇₂[3]) v₁[3]*(U̇₁[3]-U̇₂[3])]
        
                    RₑSGᵀḊ = Rₑ * skew(Gᵀ¹ * U̇₁ + Gᵀ² * Ẇ₁ + Gᵀ³ * U̇₂ + Gᵀ⁴ * Ẇ₂)
                    𝓐₃¹ = N₇lₙ² * Rₑ * A₁Ḋr¹ + RₑSGᵀḊ * P₁P¹ * Rₑ'
                    𝓐₃² =                      RₑSGᵀḊ * P₁P² * Rₑ'
                    𝓐₃⁴ =                      RₑSGᵀḊ * P₁P⁴ * Rₑ'
        
                    RₑH₁SḊ = Rₑ*H₁F₁
                    𝓐₄¹ = RₑH₁SḊ * Gᵀ¹Rₑᵀ
                    𝓐₄² = RₑH₁SḊ * Gᵀ²Rₑᵀ
                    𝓐₄⁴ = RₑH₁SḊ * Gᵀ⁴Rₑᵀ
        
                    ∂u̇₀∂d¹ = 𝓐₂¹ + 𝓐₃¹ + 𝓐₄¹
                    ∂u̇₀∂d² = 𝓐₂² + 𝓐₃² + 𝓐₄²
                    ∂u̇₀∂d³ = -∂u̇₀∂d¹
                    ∂u̇₀∂d⁴ = 𝓐₂⁴ + 𝓐₃⁴ + 𝓐₄⁴


                    nn = ∂gₙ∂x*∂gₙ∂x'
        
                    ∂u̇ₙ∂d¹ = 𝓐₁¹ + nn * ∂u̇₀∂d¹
                    ∂u̇ₙ∂d² = 𝓐₁² + nn * ∂u̇₀∂d²
                    ∂u̇ₙ∂d³ = 𝓐₁³ + nn * ∂u̇₀∂d³
                    ∂u̇ₙ∂d⁴ = 𝓐₁⁴ + nn * ∂u̇₀∂d⁴

                    # u̇ₜ = u̇₀ - u̇ₙ
                    ∂u̇ₜ∂d¹ = ∂u̇₀∂d¹ - ∂u̇ₙ∂d¹
                    ∂u̇ₜ∂d² = ∂u̇₀∂d² - ∂u̇ₙ∂d²
                    ∂u̇ₜ∂d³ = ∂u̇₀∂d³ - ∂u̇ₙ∂d³
                    ∂u̇ₜ∂d⁴ = ∂u̇₀∂d⁴ - ∂u̇ₙ∂d⁴


                    ∂u̇ₙ∂ḋ¹ = nn * RₑH₁¹Rₑᵀ
                    ∂u̇ₙ∂ḋ² = nn * RₑH₁²Rₑᵀ
                    ∂u̇ₙ∂ḋ³ = nn * RₑH₁³Rₑᵀ
                    ∂u̇ₙ∂ḋ⁴ = nn * RₑH₁⁴Rₑᵀ

                    ∂u̇₀∂ḋ¹ = RₑH₁¹Rₑᵀ
                    ∂u̇₀∂ḋ² = RₑH₁²Rₑᵀ
                    ∂u̇₀∂ḋ³ = RₑH₁³Rₑᵀ
                    ∂u̇₀∂ḋ⁴ = RₑH₁⁴Rₑᵀ


                    # u̇ₜ = u̇₀ - u̇ₙ
                    ∂u̇ₜ∂ḋ¹ = ∂u̇₀∂ḋ¹ - ∂u̇ₙ∂ḋ¹
                    ∂u̇ₜ∂ḋ² = ∂u̇₀∂ḋ² - ∂u̇ₙ∂ḋ²
                    ∂u̇ₜ∂ḋ³ = ∂u̇₀∂ḋ³ - ∂u̇ₙ∂ḋ³
                    ∂u̇ₜ∂ḋ⁴ = ∂u̇₀∂ḋ⁴ - ∂u̇ₙ∂ḋ⁴


                    aux = kₙ*p′ₙ*nn + kₙ*pₙ*∂²gₙ∂x² - η′ₙ*u̇ₙ*∂gₙ∂x'
                    Kᶠⁿ¹ = aux * RₑH₁¹Rₑᵀ - ηₙ*∂u̇ₙ∂d¹
                    Kᶠⁿ² = aux * RₑH₁²Rₑᵀ - ηₙ*∂u̇ₙ∂d²
                    Kᶠⁿ³ = aux * RₑH₁³Rₑᵀ - ηₙ*∂u̇ₙ∂d³
                    Kᶠⁿ⁴ = aux * RₑH₁⁴Rₑᵀ - ηₙ*∂u̇ₙ∂d⁴

                    Cᶠⁿ¹ = - ηₙ * nn * ∂u̇₀∂ḋ¹
                    Cᶠⁿ² = - ηₙ * nn * ∂u̇₀∂ḋ²
                    Cᶠⁿ³ = - ηₙ * nn * ∂u̇₀∂ḋ³
                    Cᶠⁿ⁴ = - ηₙ * nn * ∂u̇₀∂ḋ⁴


                    Kᶠᵗ¹ = zeros(Mat33)
                    Kᶠᵗ² = zeros(Mat33)
                    Kᶠᵗ³ = zeros(Mat33)
                    Kᶠᵗ⁴ = zeros(Mat33)

                    Cᶠᵗ¹ = zeros(Mat33)
                    Cᶠᵗ² = zeros(Mat33)
                    Cᶠᵗ³ = zeros(Mat33)
                    Cᶠᵗ⁴ = zeros(Mat33)

                    if gausspoints.status[(iB-1)*3 + iG] == 1


                        aux1 = - μʳᵉᵍ * kₙ * p′ₙ * u̇ₜ * ∂gₙ∂x'
                        aux2 = - μʳᵉᵍ * kₙ * pₙ*(ID3 - 1/(u̇ₜ²+εᵗ)*u̇ₜ*u̇ₜ')
                        Kᶠᵗ¹ = aux1 * RₑH₁¹Rₑᵀ + aux2 * ∂u̇ₜ∂d¹
                        Kᶠᵗ² = aux1 * RₑH₁²Rₑᵀ + aux2 * ∂u̇ₜ∂d²
                        Kᶠᵗ³ = aux1 * RₑH₁³Rₑᵀ + aux2 * ∂u̇ₜ∂d³
                        Kᶠᵗ⁴ = aux1 * RₑH₁⁴Rₑᵀ + aux2 * ∂u̇ₜ∂d⁴

                        Cᶠᵗ¹ = aux2 * ∂u̇ₜ∂ḋ¹
                        Cᶠᵗ² = aux2 * ∂u̇ₜ∂ḋ²
                        Cᶠᵗ³ = aux2 * ∂u̇ₜ∂ḋ³
                        Cᶠᵗ⁴ = aux2 * ∂u̇ₜ∂ḋ⁴

                    elseif gausspoints.status[(iB-1)*3 + iG] == 2

                        aux = - kₜ * δₜ/sqrt(u̇ₜ²+εᵗ) * (ID3 - 1/(u̇ₜ²+εᵗ)*u̇ₜ*u̇ₜ' + ηₜ*ID3)
                        Kᶠᵗ¹ = aux * ∂u̇ₜ∂d¹
                        Kᶠᵗ² = aux * ∂u̇ₜ∂d²
                        Kᶠᵗ³ = aux * ∂u̇ₜ∂d³
                        Kᶠᵗ⁴ = aux * ∂u̇ₜ∂d⁴
    
                        Cᶠᵗ¹ = aux * ∂u̇ₜ∂ḋ¹
                        Cᶠᵗ² = aux * ∂u̇ₜ∂ḋ²
                        Cᶠᵗ³ = aux * ∂u̇ₜ∂ḋ³
                        Cᶠᵗ⁴ = aux * ∂u̇ₜ∂ḋ⁴

                    end 

                    Kᶠᶜ¹ = Kᶠⁿ¹ + Kᶠᵗ¹
                    Kᶠᶜ² = Kᶠⁿ² + Kᶠᵗ²
                    Kᶠᶜ³ = Kᶠⁿ³ + Kᶠᵗ³
                    Kᶠᶜ⁴ = Kᶠⁿ⁴ + Kᶠᵗ⁴


                    t₄¹¹ = RₑH₁¹Rₑᵀ' * Kᶠᶜ¹
                    t₄¹² = RₑH₁¹Rₑᵀ' * Kᶠᶜ²
                    t₄¹³ = RₑH₁¹Rₑᵀ' * Kᶠᶜ³
                    t₄¹⁴ = RₑH₁¹Rₑᵀ' * Kᶠᶜ⁴
        
                    t₄²¹ = RₑH₁²Rₑᵀ' * Kᶠᶜ¹
                    t₄²² = RₑH₁²Rₑᵀ' * Kᶠᶜ²
                    t₄²³ = RₑH₁²Rₑᵀ' * Kᶠᶜ³
                    t₄²⁴ = RₑH₁²Rₑᵀ' * Kᶠᶜ⁴
        
                    t₄³¹ = RₑH₁³Rₑᵀ' * Kᶠᶜ¹
                    t₄³² = RₑH₁³Rₑᵀ' * Kᶠᶜ²
                    t₄³³ = RₑH₁³Rₑᵀ' * Kᶠᶜ³
                    t₄³⁴ = RₑH₁³Rₑᵀ' * Kᶠᶜ⁴
        
                    t₄⁴¹ = RₑH₁⁴Rₑᵀ' * Kᶠᶜ¹
                    t₄⁴² = RₑH₁⁴Rₑᵀ' * Kᶠᶜ²
                    t₄⁴³ = RₑH₁⁴Rₑᵀ' * Kᶠᶜ³
                    t₄⁴⁴ = RₑH₁⁴Rₑᵀ' * Kᶠᶜ⁴
        
        
                    Cᶠᶜ¹ = Cᶠⁿ¹ + Cᶠᵗ¹
                    Cᶠᶜ² = Cᶠⁿ² + Cᶠᵗ²
                    Cᶠᶜ³ = Cᶠⁿ³ + Cᶠᵗ³
                    Cᶠᶜ⁴ = Cᶠⁿ⁴ + Cᶠᵗ⁴
        
        
                    Kᶜ¹¹ +=  ωᴳ * (t₁¹¹ + t₂¹¹ + t₃¹¹ + t₄¹¹)
                    Kᶜ¹² +=  ωᴳ * (t₁¹² + t₂¹² + t₃¹² + t₄¹²)
                    Kᶜ¹³ +=  ωᴳ * (t₁¹³ + t₂¹³ + t₃¹³ + t₄¹³)
                    Kᶜ¹⁴ +=  ωᴳ * (t₁¹⁴ + t₂¹⁴ + t₃¹⁴ + t₄¹⁴)
        
                    Kᶜ²¹ +=  ωᴳ * (t₁²¹ + t₂²¹ + t₃²¹ + t₄²¹)
                    Kᶜ²² +=  ωᴳ * (t₁²² + t₂²² + t₃²² + t₄²²)
                    Kᶜ²³ +=  ωᴳ * (t₁²³ + t₂²³ + t₃²³ + t₄²³)
                    Kᶜ²⁴ +=  ωᴳ * (t₁²⁴ + t₂²⁴ + t₃²⁴ + t₄²⁴)
        
                    Kᶜ³¹ +=  ωᴳ * (t₁³¹ + t₂³¹ + t₃³¹ + t₄³¹)
                    Kᶜ³² +=  ωᴳ * (t₁³² + t₂³² + t₃³² + t₄³²)
                    Kᶜ³³ +=  ωᴳ * (t₁³³ + t₂³³ + t₃³³ + t₄³³)
                    Kᶜ³⁴ +=  ωᴳ * (t₁³⁴ + t₂³⁴ + t₃³⁴ + t₄³⁴)
        
                    Kᶜ⁴¹ +=  ωᴳ * (t₁⁴¹ + t₂⁴¹ + t₃⁴¹ + t₄⁴¹)
                    Kᶜ⁴² +=  ωᴳ * (t₁⁴² + t₂⁴² + t₃⁴² + t₄⁴²)
                    Kᶜ⁴³ +=  ωᴳ * (t₁⁴³ + t₂⁴³ + t₃⁴³ + t₄⁴³)
                    Kᶜ⁴⁴ +=  ωᴳ * (t₁⁴⁴ + t₂⁴⁴ + t₃⁴⁴ + t₄⁴⁴)
        
                    Cᶜ¹¹ +=  ωᴳ * RₑH₁¹Rₑᵀ' * Cᶠᶜ¹
                    Cᶜ¹² +=  ωᴳ * RₑH₁¹Rₑᵀ' * Cᶠᶜ²
                    Cᶜ¹³ +=  ωᴳ * RₑH₁¹Rₑᵀ' * Cᶠᶜ³
                    Cᶜ¹⁴ +=  ωᴳ * RₑH₁¹Rₑᵀ' * Cᶠᶜ⁴
        
                    Cᶜ²¹ +=  ωᴳ * RₑH₁²Rₑᵀ' * Cᶠᶜ¹
                    Cᶜ²² +=  ωᴳ * RₑH₁²Rₑᵀ' * Cᶠᶜ²
                    Cᶜ²³ +=  ωᴳ * RₑH₁²Rₑᵀ' * Cᶠᶜ³
                    Cᶜ²⁴ +=  ωᴳ * RₑH₁²Rₑᵀ' * Cᶠᶜ⁴
        
                    Cᶜ³¹ +=  ωᴳ * RₑH₁³Rₑᵀ' * Cᶠᶜ¹
                    Cᶜ³² +=  ωᴳ * RₑH₁³Rₑᵀ' * Cᶠᶜ²
                    Cᶜ³³ +=  ωᴳ * RₑH₁³Rₑᵀ' * Cᶠᶜ³
                    Cᶜ³⁴ +=  ωᴳ * RₑH₁³Rₑᵀ' * Cᶠᶜ⁴
        
                    Cᶜ⁴¹ +=  ωᴳ * RₑH₁⁴Rₑᵀ' * Cᶠᶜ¹
                    Cᶜ⁴² +=  ωᴳ * RₑH₁⁴Rₑᵀ' * Cᶠᶜ²
                    Cᶜ⁴³ +=  ωᴳ * RₑH₁⁴Rₑᵀ' * Cᶠᶜ³
                    Cᶜ⁴⁴ +=  ωᴳ * RₑH₁⁴Rₑᵀ' * Cᶠᶜ⁴
                
                
                else
                    gausspoints.δₜ[(iB-1)*3 + iG] = 0 
                    gausspoints.status[(iB-1)*3 + iG] = 0 
                end


            end


            
        end

        l₀2 = l₀/2
        l₀2Rₑ = l₀2 * Rₑ


        Tᵏ¹ = l₀2Rₑ*Tᵏ¹
        Tᵏ² = l₀2Rₑ*Tᵏ²
        Tᵏ³ = l₀2Rₑ*Tᵏ³
        Tᵏ⁴ = l₀2Rₑ*Tᵏ⁴



        M¹¹ = l₀2Rₑ * M¹¹ * Rₑ'
        M¹² = l₀2Rₑ * M¹² * Rₑ'
        M¹³ = l₀2Rₑ * M¹³ * Rₑ'
        M¹⁴ = l₀2Rₑ * M¹⁴ * Rₑ'
        M²² = l₀2Rₑ * M²² * Rₑ'
        M²³ = l₀2Rₑ * M²³ * Rₑ'
        M²⁴ = l₀2Rₑ * M²⁴ * Rₑ'
        M³³ = l₀2Rₑ * M³³ * Rₑ'
        M³⁴ = l₀2Rₑ * M³⁴ * Rₑ'
        M⁴⁴ = l₀2Rₑ * M⁴⁴ * Rₑ'

        Cᵏ¹¹ = l₀2Rₑ * Cᵏ¹¹ * Rₑ' 
        Cᵏ¹² = l₀2Rₑ * Cᵏ¹² * Rₑ' 
        Cᵏ¹³ = -Cᵏ¹¹             
        Cᵏ¹⁴ = l₀2Rₑ * Cᵏ¹⁴ * Rₑ' 
        Cᵏ²¹ = l₀2Rₑ * Cᵏ²¹ * Rₑ' 
        Cᵏ²² = l₀2Rₑ * Cᵏ²² * Rₑ' 
        Cᵏ²³ = -Cᵏ²¹            
        Cᵏ²⁴ = l₀2Rₑ * Cᵏ²⁴ * Rₑ' 
        Cᵏ³¹ = l₀2Rₑ * Cᵏ³¹ * Rₑ'
        Cᵏ³² = l₀2Rₑ * Cᵏ³² * Rₑ' 
        Cᵏ³³ = -Cᵏ³¹              
        Cᵏ³⁴ = l₀2Rₑ * Cᵏ³⁴ * Rₑ' 
        Cᵏ⁴¹ = l₀2Rₑ * Cᵏ⁴¹ * Rₑ' 
        Cᵏ⁴² = l₀2Rₑ * Cᵏ⁴² * Rₑ' 
        Cᵏ⁴³ = -Cᵏ⁴¹
        Cᵏ⁴⁴ = l₀2Rₑ * Cᵏ⁴⁴ * Rₑ' 



        kinetic_energy = l₀2*kinetic_energy

        if damping>0
            Cᵏ¹¹ += damping*M¹¹
            Cᵏ¹² += damping*M¹²
            Cᵏ¹³ += damping*M¹³
            Cᵏ¹⁴ += damping*M¹⁴
            Cᵏ²¹ += damping*M¹²'
            Cᵏ²² += damping*M²²
            Cᵏ²³ += damping*M²³
            Cᵏ²⁴ += damping*M²⁴
            Cᵏ³¹ += damping*M¹³'
            Cᵏ³² += damping*M²³'
            Cᵏ³³ += damping*M³³
            Cᵏ³⁴ += damping*M³⁴
            Cᵏ⁴¹ += damping*M¹⁴'
            Cᵏ⁴² += damping*M²⁴'
            Cᵏ⁴³ += damping*M³⁴'
            Cᵏ⁴⁴ += damping*M⁴⁴
        end



        if exact
            
            Θ₁ = toangle(ΔR₁)
            Θ₂ = toangle(ΔR₂)
            
            Tₛ⁻¹Θ₁ = Tₛ⁻¹(Θ₁)
            Tₛ⁻¹Θ₂ = Tₛ⁻¹(Θ₂)

            M²¹ = M¹²' * Tₛ⁻¹Θ₁'
            M²² = M²²  * Tₛ⁻¹Θ₁' 
            M²³ = M²³  * Tₛ⁻¹Θ₁'
            M²⁴ = M²⁴  * Tₛ⁻¹Θ₁'
            M³¹ = M¹³'
            M³² = M²³'
            M⁴¹ = M¹⁴' * Tₛ⁻¹Θ₂'
            M⁴² = M²⁴' * Tₛ⁻¹Θ₂' 
            M⁴³ = M³⁴' * Tₛ⁻¹Θ₂'
            M⁴⁴ = M⁴⁴  * Tₛ⁻¹Θ₂'

            Cᵏ²¹ = Cᵏ²¹ * Tₛ⁻¹Θ₁'
            Cᵏ²² = Cᵏ²² * Tₛ⁻¹Θ₁' 
            Cᵏ²³ = Cᵏ²³ * Tₛ⁻¹Θ₁'
            Cᵏ²⁴ = Cᵏ²⁴ * Tₛ⁻¹Θ₁'
            Cᵏ⁴¹ = Cᵏ⁴¹ * Tₛ⁻¹Θ₂'
            Cᵏ⁴² = Cᵏ⁴² * Tₛ⁻¹Θ₂' 
            Cᵏ⁴³ = Cᵏ⁴³ * Tₛ⁻¹Θ₂'
            Cᵏ⁴⁴ = Cᵏ⁴⁴ * Tₛ⁻¹Θ₂'

        else

            M²¹ = M¹²'
            M³¹ = M¹³'
            M³² = M²³'
            M⁴¹ = M¹⁴' 
            M⁴² = M²⁴'
            M⁴³ = M³⁴'

        end




        if contactsearch


            Tᶜ¹ = l₀2*Tᶜ¹
            Tᶜ² = l₀2*Tᶜ²
            Tᶜ³ = l₀2*Tᶜ³
            Tᶜ⁴ = l₀2*Tᶜ⁴


            Cᶜ¹¹ = l₀2 * Cᶜ¹¹ 
            Cᶜ¹² = l₀2 * Cᶜ¹² 
            Cᶜ¹³ = l₀2 * Cᶜ¹³ 
            Cᶜ¹⁴ = l₀2 * Cᶜ¹⁴ 
            Cᶜ²¹ = l₀2 * Cᶜ²¹ 
            Cᶜ²² = l₀2 * Cᶜ²² 
            Cᶜ²³ = l₀2 * Cᶜ²³ 
            Cᶜ²⁴ = l₀2 * Cᶜ²⁴ 
            Cᶜ³¹ = l₀2 * Cᶜ³¹ 
            Cᶜ³² = l₀2 * Cᶜ³² 
            Cᶜ³³ = l₀2 * Cᶜ³³ 
            Cᶜ³⁴ = l₀2 * Cᶜ³⁴ 
            Cᶜ⁴¹ = l₀2 * Cᶜ⁴¹ 
            Cᶜ⁴² = l₀2 * Cᶜ⁴² 
            Cᶜ⁴³ = l₀2 * Cᶜ⁴³ 
            Cᶜ⁴⁴ = l₀2 * Cᶜ⁴⁴ 

            Kᶜ¹¹ = l₀2 * Kᶜ¹¹
            Kᶜ¹² = l₀2 * Kᶜ¹²
            Kᶜ¹³ = l₀2 * Kᶜ¹³
            Kᶜ¹⁴ = l₀2 * Kᶜ¹⁴
            Kᶜ²¹ = l₀2 * Kᶜ²¹
            Kᶜ²² = l₀2 * Kᶜ²²
            Kᶜ²³ = l₀2 * Kᶜ²³
            Kᶜ²⁴ = l₀2 * Kᶜ²⁴
            Kᶜ³¹ = l₀2 * Kᶜ³¹
            Kᶜ³² = l₀2 * Kᶜ³²
            Kᶜ³³ = l₀2 * Kᶜ³³
            Kᶜ³⁴ = l₀2 * Kᶜ³⁴
            Kᶜ⁴¹ = l₀2 * Kᶜ⁴¹
            Kᶜ⁴² = l₀2 * Kᶜ⁴²
            Kᶜ⁴³ = l₀2 * Kᶜ⁴³
            Kᶜ⁴⁴ = l₀2 * Kᶜ⁴⁴

            contact_energy = l₀2 * contact_energy


            if exact

                Cᶜ²¹ = Cᶜ²¹ * Tₛ⁻¹Θ₁'
                Cᶜ²² = Cᶜ²² * Tₛ⁻¹Θ₁' 
                Cᶜ²³ = Cᶜ²³ * Tₛ⁻¹Θ₁'
                Cᶜ²⁴ = Cᶜ²⁴ * Tₛ⁻¹Θ₁'
                Cᶜ⁴¹ = Cᶜ⁴¹ * Tₛ⁻¹Θ₂'
                Cᶜ⁴² = Cᶜ⁴² * Tₛ⁻¹Θ₂' 
                Cᶜ⁴³ = Cᶜ⁴³ * Tₛ⁻¹Θ₂'
                Cᶜ⁴⁴ = Cᶜ⁴⁴ * Tₛ⁻¹Θ₂'

            end


        end



    end


    Tᵏ = [Tᵏ¹; Tᵏ²; Tᵏ³; Tᵏ⁴]
    
    M = hcat(vcat(M¹¹, M²¹, M³¹, M⁴¹), vcat(M¹², M²², M³², M⁴²), vcat(M¹³, M²³, M³³, M⁴³), vcat(M¹⁴, M²⁴, M³⁴, M⁴⁴))
    Cᵏ = hcat(vcat(Cᵏ¹¹, Cᵏ²¹, Cᵏ³¹, Cᵏ⁴¹), vcat(Cᵏ¹², Cᵏ²², Cᵏ³², Cᵏ⁴²), vcat(Cᵏ¹³, Cᵏ²³, Cᵏ³³, Cᵏ⁴³), vcat(Cᵏ¹⁴, Cᵏ²⁴, Cᵏ³⁴, Cᵏ⁴⁴))

    
    Tᶜ = [Tᶜ¹; Tᶜ²; Tᶜ³; Tᶜ⁴]
    Kᶜ = hcat(vcat(Kᶜ¹¹, Kᶜ²¹, Kᶜ³¹, Kᶜ⁴¹), vcat(Kᶜ¹², Kᶜ²², Kᶜ³², Kᶜ⁴²), vcat(Kᶜ¹³, Kᶜ²³, Kᶜ³³, Kᶜ⁴³), vcat(Kᶜ¹⁴, Kᶜ²⁴, Kᶜ³⁴, Kᶜ⁴⁴))
    Cᶜ = hcat(vcat(Cᶜ¹¹, Cᶜ²¹, Cᶜ³¹, Cᶜ⁴¹), vcat(Cᶜ¹², Cᶜ²², Cᶜ³², Cᶜ⁴²), vcat(Cᶜ¹³, Cᶜ²³, Cᶜ³³, Cᶜ⁴³), vcat(Cᶜ¹⁴, Cᶜ²⁴, Cᶜ³⁴, Cᶜ⁴⁴))    

    return strain_energy, kinetic_energy, contact_energy, Tⁱⁿᵗ, Tᵏ, Tᶜ, Kⁱⁿᵗ, Kᶜ, M, Cᵏ, Cᶜ




    
end

function assemble!(conf, matrices, energy, params, Δt) 

    @unpack nodes, beams, colors, contact, sdf, gausspoints = conf
    
    # initialise the matrices associate to the whole structure
    fill!(matrices.K, 0)
    fill!(matrices.C, 0)
    fill!(matrices.M, 0)
    fill!(matrices.Tⁱⁿᵗ, 0)
    fill!(matrices.Tᵏ, 0)
    fill!(matrices.Tᶜ, 0)
    
    # initialise the energy values associate to the whole structure
    energy.strain_energy = 0
    energy.kinetic_energy = 0
    energy.contact_energy = 0
    energy.max_strain_energy = 0.0

    #Preparing constants
    gausspoints_info = (params.nᴳ, params.ωᴳ, params.zᴳ)

    for cidxs in colors

        @batch for idx in cidxs

            b = LazyRow(beams, idx)
            
            n1 = b.node1
            n2 = b.node2
            # information from node 1 and 2
            X₁, X₂ = nodes.X₀[n1], nodes.X₀[n2]
            u₁, u₂ = nodes.u[n1], nodes.u[n2]
            u̇₁, u̇₂ = nodes.u̇[n1], nodes.u̇[n2]
            ü₁, ü₂ = nodes.ü[n1], nodes.ü[n2]
            ẇ₁, ẇ₂ = nodes.ẇ[n1], nodes.ẇ[n2]
            ẅ₁, ẅ₂ = nodes.ẅ[n1], nodes.ẅ[n2]
            R₁, R₂ = nodes.R[n1], nodes.R[n2]
            ΔR₁, ΔR₂ = nodes.ΔR[n1], nodes.ΔR[n2]

            # Packing
            init = (b.ind, X₁, X₂, b.l₀, b.Rₑ⁰)      
            constants = (init, gausspoints_info, b.properties, contact, sdf, Δt, gausspoints)

            strain_energy, kinetic_energy, contact_energy, Tⁱⁿᵗ, Tᵏ, Tᶜ, Kⁱⁿᵗ, Kᶜ, M, Cᵏ, Cᶜ = compute(u₁, u₂, R₁, R₂, ΔR₁, ΔR₂, u̇₁, u̇₂, ẇ₁, ẇ₂, ü₁, ü₂, ẅ₁, ẅ₂, constants)
        
            K = Kⁱⁿᵗ - Kᶜ
            C = Cᵏ-(1+params.α)*Cᶜ

            #-----------------------
            # Assemble contributions

            idof1 = nodes.idof_6[n1]
            idof2 = nodes.idof_6[n2]
            
            idofs = vcat(idof1, idof2)

            energy.strain_energy +=  strain_energy
            energy.kinetic_energy += kinetic_energy
            energy.contact_energy += contact_energy
            
            energy.max_strain_energy = max(
                energy.max_strain_energy,
                strain_energy
            )


            matrices.Tᵏ[idofs] += Tᵏ
            matrices.Tⁱⁿᵗ[idofs] += Tⁱⁿᵗ
            matrices.Tᶜ[idofs] += Tᶜ

            matrices.K[idofs, idofs] += K
            matrices.C[idofs, idofs] += C
            matrices.M[idofs, idofs] += M
                                
        end

    end

end 