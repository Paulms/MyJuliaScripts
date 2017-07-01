type LSFractal
  axiom::String
  rules::Dict{Char,String}
  angle::Real
  n::Int
  points::Matrix{Real}
end

function LSFractal(axiom, rules, angle=0.0,n=2;first_angle = nothing)
  axiom_ex = axiom
  if first_angle == nothing
    first_angle = angle
  end
  for i in 2:n
    axiom_ex = join([(c in keys(rules)) ? rules[c] : c  for c in axiom_ex])
  end
  points = zeros(1,3)
  points[1,3] = first_angle
  for action in axiom_ex
    if action in 'A':'Z'
      x = points[end,1] + cos(points[end,3]*(π/180))
      y = points[end,2] + sin(points[end,3]*(π/180))
      alfa = points[end,3]
      points = vcat(points,[x y alfa])
    else
      alfa = points[end,3]
      points[end, 3]=eval(parse("$action($alfa,$angle)"))
    end
  end
  LSFractal(axiom,rules,angle,n,points)
end
import GR: plot
function plot(fractal::LSFractal)
  l_x = minimum(fractal.points[:,1]); h_x = maximum(fractal.points[:,1])
  l_y = minimum(fractal.points[:,2]); h_y = maximum(fractal.points[:,2])
  w_x = h_x - l_x; w_y =h_y - l_y
  w_xy = max(w_x,w_y)
  setwindow(l_x, l_x + w_xy, l_y, l_y+w_xy)
  polyline(fractal.points[:,1],fractal.points[:,2])
end

using GR
inline("atom")
# Snowflake
axiom="F++F++F"
rules=Dict('F' => "F-F++F-F")
angle=60.0
n=5
snowflake2 = LSFractal(axiom,rules,angle,n)
plot(snowflake2)
GR.show()

#Sierpinsky triangle
axiom="R"
rules=Dict('L' => "R+L+R", 'R' => "L-R-L")
angle=-60.0
n=6
Striangle = LSFractal(axiom,rules,angle,n;first_angle = 0.0)
plot(Striangle)
GR.show()

#Koch Island
axiom="F-F-F-F"
rules=Dict('F' => "F-F+F+FF-F-F+F")
angle=90
n = 3
KochIsland = LSFractal(axiom,rules,angle,n)
plot(KochIsland)
GR.show()

#Koch Island
axiom="F-F-F-F"
rules=Dict('F' => "F-F+F+FF-F-F+F")
angle=90
n = 3
KochIsland = LSFractal(axiom,rules,angle,n)
plot(KochIsland)
GR.show()

# Hexagonal Gosper Curve
axiom="L"
rules=Dict('L'=>"L+R++R-L--LL-R+", 'R'=>"-L+RR++R+L--L-R")
angle=60
n=4
HGosper = LSFractal(axiom,rules,angle,n)
plot(HGosper)
GR.show()

# Cuadratic Gosper Curve
axiom="-R"
rules=Dict('L'=>"LL-R-R+L+L-R-RL+R+LLR-L+R+LL+R-LR-R-L+L+RR-",
           'R'=>"+LL-R-R+L+LR+L-RR-L-R+LRR-L-RL+L+R-R-L+L+RR")
angle=90
n=3
CGosper = LSFractal(axiom,rules,angle,n)
plot(CGosper)
GR.show()

#Sierpinsky Carpet
axiom="F"
rules=Dict('F'=>"F+F-F-F-G+F+F+F-F", 'G'=>"GGG")
angle=90
n=5
Carpet = LSFractal(axiom,rules,angle,n;first_angle = 45.0)
plot(Carpet)
GR.show()

#Sierpinsky Triangle 2
axiom="F-G-G"
rules=Dict('G'=>"GG", 'F'=>"F-G+F+G-F")
angle=120
n=6
STriangle2 = LSFractal(axiom,rules,angle,n)
plot(STriangle2)
GR.show()

#Sierpinsky Median Curve
axiom="L--F--L--F"
rules=Dict('L'=>"+R-F-R+", 'R'=>"-L+F+L-")
angle=45
n=8
SMedian = LSFractal(axiom,rules,angle,n)
plot(SMedian)
GR.show()
