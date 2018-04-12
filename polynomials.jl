using Polynomials
"""
    jacobi(p::Integer, α, β)
Compute the Legendre polynomial with parameters `α`, `β` of degree `p`
using the three term recursion [Karniadakis and Sherwin, Spectral/hp Element
Methods for CFD, Appendix A].
"""
function poly_jacobi(p::Integer, α, β, ::Type{T}=Float64, var=:x) where {T}
    a = Poly{T}([one(T)], var)
    b = Poly{T}([(α - β)/2, (2+α+β)/2], var)
    if p <= 0
        return a
    elseif p == 1
        return b
    end
    px = Poly{T}([zero(T), one(T)], var)
    for n in 2:p
        a1 = 2n*(n+α+β)*(2n-2+α+β)
        a2 = (2n-1+α+β)*(α+β)*(α-β)
        a3 = (2n-2+α+β)*(2n-1+α+β)*(2n+α+β)
        a4 = 2*(n-1+α)*(n-1+β)*(2n+α+β)
        a, b = b, ( (a2+a3*px)*b - a4*a ) / a1
    end
    b
end

function poly_legendre{T<:Number}(n, ::Type{T}=Float64, var=:x)
    return poly_jacobi(n,0.0,0.0,T,var)
end

function dubiner_base(j::Integer,z,w)
    S=2*w-1
    R=(2*(1+(2*z-1))/(1-S))-1
    t=-3/2+(1/2)*sqrt(1+8*j)
    n=((ceil(t)+1)*(ceil(t)+2))/2-j;
    k=2*n+1;
    m=ceil(t)-n;

    @assert ((w>=0)&&(z>=0))&&(1>=z+w) "evalued point not in domain"
    P=polyval(poly_jacobi(Int(n),0,0),R)*polyval(poly_jacobi(Int(m),k,0),S)*((1-S)/2)^n
    N=sqrt(2/((2*n+1)*(m+n+1)))
    return (2*P)/N
end

function poly_jacobi2(x, p::Integer, α, β)
    T = typeof( (2+α+β)*x / 2 )
    a = one(T)
    b = ((2+α+β)*x + α - β) / 2
    if p <= 0
        return a
    elseif p == 1
        return b
    end
    for n in 2:p
        a1 = 2n*(n+α+β)*(2n-2+α+β)
        a2 = (2n-1+α+β)*(α+β)*(α-β)
        a3 = (2n-2+α+β)*(2n-1+α+β)*(2n+α+β)
        a4 = 2*(n-1+α)*(n-1+β)*(2n+α+β)
        a, b = b, ( (a2+a3*x)*b - a4*a ) / a1
    end
    b
end

function poly_legendre2(x::T, n::Integer) where {T}
    return poly_jacobi2(x, n,zero(T),zero(T))
end

function dubiner_base2(j::Integer,z,w)
    S=2*w-1
    R=(2*(1+(2*z-1))/(1-S))-1
    t=-3/2+(1/2)*sqrt(1+8*j)
    n=((ceil(t)+1)*(ceil(t)+2))/2-j;
    k=2*n+1;
    m=ceil(t)-n;
    @assert ((w>=0)&&(z>=0))&&(1>=z+w) "evalued point not in domain"
    P=poly_jacobi2(R,Int(n),0,0)*poly_jacobi2(S, Int(m),k,0)*((1-S)/2)^n
    N=sqrt(2/((2*n+1)*(m+n+1)))
    return (2*P)/N
end

using PolynomialBases
function dubiner_base3(j::Integer,z,w)
    S=2*w-1
    R=(2*(1+(2*z-1))/(1-S))-1
    t=-3/2+(1/2)*sqrt(1+8*j)
    n=((ceil(t)+1)*(ceil(t)+2))/2-j;
    k=2*n+1;
    m=ceil(t)-n;
    @assert ((w>=0)&&(z>=0))&&(1>=z+w) "evalued point not in domain"
    P=jacobi(R,Int(n),0,0)*jacobi(S, Int(m),k,0)*((1-S)/2)^n
    N=sqrt(2/((2*n+1)*(m+n+1)))
    return (2*P)/N
end
@code_warntype poly_jacobi2(1.0,1,0.0,0.0)
