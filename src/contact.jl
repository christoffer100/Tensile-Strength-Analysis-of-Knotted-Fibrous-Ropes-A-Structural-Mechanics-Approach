abstract type SignedDistanceField end

struct Plane_z_SDF <: SignedDistanceField
    
    r::Float64  # beams radius
    z0::Float64 # plane position 
    
end 

struct Plane_y_SDF <: SignedDistanceField
    
    r::Float64  # beams radius
    y0::Float64 # plane position 
    
end 

struct Sphere_SDF <: SignedDistanceField
    
    r::Float64 # beams radius
    R::Float64 # sphere radius
    x0::Float64 # sphere position
    y0::Float64 # sphere position
    z0::Float64 # sphere position
    
end 

mutable struct Cylinder_SDF <: SignedDistanceField
    
    r::Float64 # beams radius
    R::Float64 # cylinder radius
    
end  

# Properties of a discrete SDF
struct Discrete_SDF{F} <: SignedDistanceField
    
    r::Float64
    sitp::F
    dom::NTuple{6, Float64}
    dx::Float64
    dy::Float64
    dz::Float64
    
end 

function Discrete_SDF(filename, radius, inside)
    
    # Read sdf from file
    npx, npy, npz, dx, dy, dz, dom, sdf = read_VTK_sdf(filename)
    if inside 
        field = reshape(sdf, (npx,npy,npz))
    else 
        field = reshape(-sdf, (npx,npy,npz))
    end 

    # Coordinates where the sdf values are taken 
    x0 = dom[1]
    y0 = dom[3]
    z0 = dom[5]

    xend = dom[2]
    yend = dom[4]
    zend = dom[6]

    x = range(x0; step=dx, stop=xend)
    y = range(y0; step=dy, stop=yend)
    z = range(z0; step=dz, stop=zend)

    # Build quadratic interpolation (quadratic should be enough in our case to get continuous gradient)
    itp = interpolate(field, BSpline(Quadratic(Reflect(OnCell()))))

    # Scale the interpolation on the defined coordinate grid
    sitp = scale(itp, x, y, z)  

    return Discrete_SDF{typeof(sitp)}(radius, sitp, dom, dx, dy, dz)  

end 

# Get gap, gradient and hession of sphere analytical SDF at point 
@inline function contact_gap(point, sdf::Sphere_SDF)
    
    aux = Vec3(point[1] - sdf.x0, point[2] - sdf.y0, point[3] - sdf.z0)
    
    norm_aux = norm(aux)
    invnorm = 1/norm_aux
    
    gₙ = norm_aux - sdf.R - sdf.r
    ∂gₙ∂x = invnorm * aux
    ∂²gₙ∂x² = invnorm*ID3 - (invnorm^3)*(aux*aux')
    
    return gₙ, ∂gₙ∂x, ∂²gₙ∂x²
    
end

@inline function incontact(point, sdf::Sphere_SDF, ḡₙ)
    
    aux = Vec3(point[1] - sdf.x0, point[2] - sdf.y0, point[3] - sdf.z0)
    norm_aux = norm(aux)
    gₙ = norm_aux - sdf.R - sdf.r
    return gₙ - sdf.r ≤ ḡₙ
    
end

# Get gap, gradient and hession of sphere analytical SDF at point 
@inline function contact_gap(point, sdf::Cylinder_SDF)
    
    aux = Vec3(point[1], point[2], 0)
    
    norm_aux = norm(aux)
    invnorm = 1/norm_aux
    
    gₙ = norm_aux - sdf.R + sdf.r
    ∂gₙ∂x = invnorm * aux
    ∂²gₙ∂x² = invnorm*ID3 + (invnorm^3)*(aux*aux')
    
    return -gₙ, -∂gₙ∂x, -∂²gₙ∂x²
    
end

@inline function incontact(point, sdf::Cylinder_SDF, ḡₙ)
    
    aux = Vec3(point[1], point[2], 0)
    norm_aux = norm(aux)
    gₙ = norm_aux - sdf.R + sdf.r
    
    return -gₙ - sdf.r ≤ ḡₙ
    
end


@inline function isinside(point, dom)

    l_x = point[1] - dom[1]  
    flag_lx = l_x>=0 && l_x<=(dom[2]-dom[1]) 
    l_y = point[2] - dom[3]
    flag_ly = l_y>=0 && l_y<=(dom[4]-dom[3])
    l_z = point[3] - dom[5] 
    flag_lz = l_z>=0 && l_z<=(dom[6]-dom[5])

    return flag_lx && flag_ly && flag_lz   

end

@inline function contact_gap(point, sdf::Discrete_SDF)

    sitp = sdf.sitp

    gₙ = sitp(point...)
    ∂gₙ∂x = Interpolations.gradient(sitp, point...)
    ∂²gₙ∂x² = Interpolations.hessian(sitp, point...)
    
    # Normalize
    nn = dot(∂gₙ∂x, ∂gₙ∂x)
    nmaginv = 1/sqrt(nn)
    n = ∂gₙ∂x*nmaginv 
    H = ∂²gₙ∂x²*nmaginv * (ID3 - (∂gₙ∂x*∂gₙ∂x')/nn)
    
    return gₙ - sdf.r, n, H
    
    
        
end 

@inline function incontact(point, sdf::Discrete_SDF, ḡₙ)

    sitp = sdf.sitp
    return isinside(point, sdf.dom) && (sitp(point...) - sdf.r) ≤ ḡₙ
        
end 

# Quadratically regulise penalty
@inline function regularize_gₙ(gₙ::T, radius) where T
    
    ḡₙ = radius/4
    p̄ₙ = ḡₙ/2
    
    pₙ = 0
    p′ₙ = 0
    Πₑ = 0

    if gₙ≤0

        pₙ = p̄ₙ - gₙ
        p′ₙ = -one(T)
        Πₑ = gₙ^2/2 - p̄ₙ*gₙ + (ḡₙ^2)/6

    elseif ḡₙ>=gₙ && gₙ>0

        aux = (ḡₙ-p̄ₙ)/(ḡₙ^2)
        pₙ =  aux*gₙ^2 - gₙ + p̄ₙ
        p′ₙ = 2*aux*gₙ - 1
        Πₑ = (ḡₙ-p̄ₙ)/(3*ḡₙ^2)*gₙ^3 - gₙ^2/2 + p̄ₙ*gₙ - ḡₙ^2/6

    end

    return pₙ, p′ₙ, Πₑ
    
end 

@inline function smoothstep(v::T, x, xᵘ) where T
    
    if x≤0

        y = v
        y′ = zero(T)

    else

        x = (x-xᵘ)/xᵘ
        y = x * x * (3 + 2 * x) * v
        y′ = 6/xᵘ * x * (1 + x) * v

    end
    
    return y, y′
    
end 


