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
"""
eval_dubiner(z,w,j::Integer)
Compute the dubiner polynomial of degree `j` at point (z,w)
on the reference triangle ((0,0),(1,0),(0,1)) ?
"""
function eval_dubiner(z,w,j::Integer)
    S=2*w-1
    R=(2*(1+(2*z-1))/(1-S))-1
    t=-3/2+(1/2)*sqrt(1+8*j)
    n=((ceil(t)+1)*(ceil(t)+2))/2-j;
    k=2*n+1;
    m=ceil(t)-n;
    @assert ((w>=0)&&(z>=0))&&(1>=z+w) "point not in domain"
    P=jacobi(R,Int(n),0,0)*jacobi(S, Int(m),k,0)*((1-S)/2)^n
    N=sqrt(2/((2*n+1)*(m+n+1)))
    return (2*P)/N
end

"""
eval_dubiner_ex(z,w,j::Integer)
Compute the dubiner polynomial of degree `j` at point (z,w)
on the reference triangle ((0,0),(1,0),(0,1)) ?
Use exact functions when available
"""
function eval_dubiner_ex(r::T,s::T,j::Integer) where{T}
    if j == 0; return zero(T)
    elseif j == 1; return sqrt(2)*one(T)
    elseif j == 2; return 2*sqrt(3)*(2*r + s - 1)
    elseif j == 3; return 2*(3*s - 1)
    elseif j == 4; return sqrt(30)*(6*r^2 + 6*r*(s - 1) + s^2 - 2*s + 1)
    elseif j == 5; return 3*sqrt(2)*(5*s - 1)*(2*r + s - 1)
    elseif j == 6; return sqrt(6)*(10*s^2 - 8*s + 1)
    elseif j == 7; return 2*sqrt(14)*(2*r + s - 1)*(10*r^2 + 10*r*(s - 1) + s^2 - 2*s + 1)
    elseif j == 8; return 2*sqrt(10)*(7*s - 1)*(6*r^2 + 6*r*(s - 1) + s^2 - 2*s + 1)
    elseif j == 9; return 2*sqrt(6)*(21*s^2 - 12*s + 1)*(2*r + s - 1)
    elseif j == 10; return 2*sqrt(2)*(35*s^3 - 45*s^2 + 15*s - 1)
    elseif j == 11; return 3*sqrt(10)*(70*r^4 + 140*r^3*(s - 1) + 90*r^2*(s^2 - 2*s + 1) + 20*r*(s^3 - 3*s^2 + 3*s - 1) + s^4 - 4*s^3 + 6*s^2 - 4*s + 1)
    elseif j == 12; return sqrt(70)*(9*s - 1)*(2*r + s - 1)*(10*r^2 + 10*r*(s - 1) + s^2 - 2*s + 1)
    elseif j == 13; return 5*sqrt(2)*(36*s^2 - 16*s + 1)*(6*r^2 + 6*r*(s - 1) + s^2 - 2*s + 1)
    elseif j == 14; return sqrt(30)*(84*s^3 - 84*s^2 + 21*s - 1)*(2*r + s - 1)
    elseif j == 15; return sqrt(10)*(126*s^4 - 224*s^3 + 126*s^2 - 24*s + 1)
    else; return eval_dubiner(r,s,j)
    end
end

#test
r = 0.5; s=0.5
p = zeros(15,1)
p[1] = sqrt(2)
p[2] = 2*sqrt(3)*(2*r + s - 1)
p[3] = 2*(3*s - 1)
p[4] = sqrt(30)*(6*r^2 + 6*r*(s - 1) + s^2 - 2*s + 1)
p[5] = 3*sqrt(2)*(5*s - 1)*(2*r + s - 1)
p[6] = sqrt(6)*(10*s^2 - 8*s + 1)
p[7] = 2*sqrt(14)*(2*r + s - 1)*(10*r^2 + 10*r*(s - 1) + s^2 - 2*s + 1)
p[8] = 2*sqrt(10)*(7*s - 1)*(6*r^2 + 6*r*(s - 1) + s^2 - 2*s + 1)
p[9] = 2*sqrt(6)*(21*s^2 - 12*s + 1)*(2*r + s - 1)
p[10] = 2*sqrt(2)*(35*s^3 - 45*s^2 + 15*s - 1)
p[11] = 3*sqrt(10)*(70*r^4 + 140*r^3*(s - 1) + 90*r^2*(s^2 - 2*s + 1) + 20*r*(s^3 - 3*s^2 + 3*s - 1) + s^4 - 4*s^3 + 6*s^2 - 4*s + 1)
p[12] = sqrt(70)*(9*s - 1)*(2*r + s - 1)*(10*r^2 + 10*r*(s - 1) + s^2 - 2*s + 1)
p[13] = 5*sqrt(2)*(36*s^2 - 16*s + 1)*(6*r^2 + 6*r*(s - 1) + s^2 - 2*s + 1)
p[14] = sqrt(30)*(84*s^3 - 84*s^2 + 21*s - 1)*(2*r + s - 1)
p[15] = sqrt(10)*(126*s^4 - 224*s^3 + 126*s^2 - 24*s + 1)

for j = 1:15
    println(eval_dubiner(r,s,j) ≈ p[j])
end
