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
function eval_Djacobi(x,n::Integer,α, β)
    T = typeof((2+α+β)*x/2)
    a = zero(T)
    if n <= 0;return a;end
    return (n+α+β+1)/2*jacobi(x,n-1,α+1,β+1)
end
"""
eval_dubiner(z,w,j::Integer)
Compute the dubiner polynomial of degree `j` at point (x,y)
on the reference triangle ((0,0),(1,0),(0,1))
"""
function eval_dubiner(x,y,j::Integer)
    #check domain
    @assert ((y>=0)&&(x>=0))&&(1>=x+y) "point not in domain"
    # Map to reference square
    ξ=2*x/(1-y)-1
    η=2*y-1
    #Compute degrees
    t=-3/2+(1/2)*sqrt(1+8*j)
    n=((ceil(t)+1)*(ceil(t)+2))/2-j
    k=2*n+1
    m=ceil(t)-n
    #Compute Dubiner_nm(ξ, η)
    P=jacobi(ξ,Int(n),0,0)*jacobi(η, Int(m),k,0)*((1-η)/2)^n
    N=sqrt(2/((2*n+1)*(m+n+1)))
    return (2*P)/N
end

"""
eval_Ddubiner(z,w,j::Integer)
Compute the gradient of dubiner polynomial of degree `j` at point (x,y)
on the reference triangle ((0,0),(1,0),(0,1))
"""
function eval_Ddubiner(x,y,j::Integer)
    #check domain
    @assert ((y>=0)&&(x>=0))&&(1>=x+y) "point not in domain"
    # Map to reference square
    ξ=2*x/(1-y)-1
    η=2*y-1
    #Compute degrees
    t=-3/2+(1/2)*sqrt(1+8*j)
    n=Int(((ceil(t)+1)*(ceil(t)+2))/2-j)
    k=2*n+1
    m=Int(ceil(t)-n)
    #Compute ∇Dubiner_nm(ξ, η)
    Dφn = eval_Djacobi(ξ,n,0,0)
    φn = jacobi(ξ,n,0,0)
    φm = jacobi(η,m,k,0)
    Dφm = eval_Djacobi(η,m,k,0)
    Px=2/(1-η)*Dφn*φm*((1-η)/2)^n
    N=sqrt(2/((2*n+1)*(m+n+1)))
    Py=(2*x/(1-y)^2*Dφn*((1-η)/2)^n-n*((1-η)/2)^(n-1)*φn)*φm + 2*φn*((1-η)/2)^n*Dφm
    return [(4*Px)/N,(2*Py)/N]
end

"""
eval_dubiner_ex(z,w,j::Integer)
Compute the dubiner polynomial of degree `j` at point (r,s)
on the reference triangle ((0,0),(1,0),(0,1))
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

"""
eval_Ddubiner_ex(z,w,j::Integer)
Compute the gradient of dubiner polynomial of degree `j` at point (r,s)
on the reference triangle ((0,0),(1,0),(0,1))
Use exact functions when available
"""
function eval_Ddubiner_ex(r::T,s::T,j::Integer) where{T}
    if j == 0; return [zero(T),zero(T)]
    elseif j == 1; return [zero(T),zero(T)]
    elseif j == 2; return [4*sqrt(3), 2*sqrt(3)]*one(T)
    elseif j == 3; return [zero(T), one(T)*6]
    elseif j == 4; return [12*sqrt(30)*r + 6*sqrt(30)*(s - 1),
                            6*sqrt(30)*r + 2*sqrt(30)*s - 2*sqrt(30)]
    elseif j == 5; return [6*sqrt(2)*(5*s - 1),
                            30*sqrt(2)*r + 30*sqrt(2)*s - 18*sqrt(2)]
    elseif j == 6; return [zero(T),20*sqrt(6)*s - 8*sqrt(6)]
    elseif j == 7; return [24*sqrt(14)*(5*r^2 + 5*r*(s - 1) + s^2 - 2*s + 1),
                            60*sqrt(14)*r^2 + r*(48*sqrt(14)*s - 48*sqrt(14)) + 6*sqrt(14)*s^2 - 12*sqrt(14)*s + 6*sqrt(14)]
    elseif j == 8; return [12*sqrt(10)*(7*s - 1)*(2*r + s - 1),
                            84*sqrt(10)*r^2 + r*(168*sqrt(10)*s - 96*sqrt(10)) + 42*sqrt(10)*s^2 - 60*sqrt(10)*s + 18*sqrt(10)]
    elseif j == 9; return [4*sqrt(6)*(21*s^2 - 12*s + 1),
                            r*(168*sqrt(6)*s - 48*sqrt(6)) + 126*sqrt(6)*s^2 - 132*sqrt(6)*s + 26*sqrt(6)]
    elseif j == 10; return [zero(T),210*sqrt(2)*s^2 - 180*sqrt(2)*s + 30*sqrt(2)]
    elseif j == 11; return [840*sqrt(10)*r^3 + 1260*sqrt(10)*r^2*(s - 1) + 540*sqrt(10)*r*(s^2 - 2*s + 1) + 60*sqrt(10)*(s^3 - 3*s^2 + 3*s - 1),
                            420*sqrt(10)*r^3 + r^2*(540*sqrt(10)*s - 540*sqrt(10)) + r*(180*sqrt(10)*s^2 - 360*sqrt(10)*s + 180*sqrt(10)) + 12*sqrt(10)*s^3 - 36*sqrt(10)*s^2 + 36*sqrt(10)*s - 12*sqrt(10)]
    elseif j == 12; return [12*sqrt(70)*(9*s - 1)*(5*r^2 + 5*r*(s - 1) + s^2 - 2*s + 1),
                            180*sqrt(70)*r^3 + r^2*(540*sqrt(70)*s - 300*sqrt(70)) + r*(324*sqrt(70)*s^2 - 456*sqrt(70)*s + 132*sqrt(70)) + 36*sqrt(70)*s^3 - 84*sqrt(70)*s^2 + 60*sqrt(70)*s - 12*sqrt(70)]
    elseif j == 13; return [30*sqrt(2)*(36*s^2 - 16*s + 1)*(2*r + s - 1),
                            r^2*(2160*sqrt(2)*s - 480*sqrt(2)) + r*(3240*sqrt(2)*s^2 - 3120*sqrt(2)*s + 510*sqrt(2)) + 720*sqrt(2)*s^3 - 1320*sqrt(2)*s^2 + 690*sqrt(2)*s - 90*sqrt(2)]
    elseif j == 14; return [2*sqrt(30)*(84*s^3 - 84*s^2 + 21*s - 1),
                            r*(504*sqrt(30)*s^2 - 336*sqrt(30)*s + 42*sqrt(30)) + 336*sqrt(30)*s^3 - 504*sqrt(30)*s^2 + 210*sqrt(30)*s - 22*sqrt(30)]
    elseif j == 15; return [zero(T),504*sqrt(10)*s^3 - 672*sqrt(10)*s^2 + 252*sqrt(10)*s - 24*sqrt(10)]
    else; return eval_Ddubiner(r,s,j)
    end
end


#test
const rr = 0.5
const ss = 0.5
for j = 1:15
    println(eval_dubiner(rr,ss,j) ≈ eval_dubiner_ex(rr,ss,j))
    println(eval_Ddubiner(rr,ss,j) ≈ eval_Ddubiner_ex(rr,ss,j))
end
