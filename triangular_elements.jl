
using LinearAlgebra
using Symbolics
using SymPy

# Plane-stress (triangle):
# [ 1 x1 y1;
#   1 x2 y2; 
#   1 x3 y3 ]

M1 = [1. 0 0;
	1 20 10;
	1 0 10]
M2 = [1. 0 0;
       1 20 0;
       1 20 10]

@variables x y;
N1 = [1, x, y]' * inv(M1)
N2 = [1, x, y]' * inv(M2)
B1 = Symbolics.jacobian(N1, [x, y])'
B2 = Symbolics.jacobian(N2, [x, y])'

E = 30e6
V = 0.3
B1f = (x->Symbolics.symbolic_to_float(x)).(B1)
B2f = (x->Symbolics.symbolic_to_float(x)).(B2)
D = E/(1-V^2) * [1 V 0; V 1 0; 0 0 ((1-V)/2)]

# reshape for x,y
Bxy1 = zeros(3, 6)
Bxy1[1, 1:2:end] = B1f[1, :]
Bxy1[2, 2:2:end] = B1f[2, :]
Bxy1[3, 1:2:end] = B1f[2, :]
Bxy1[3, 2:2:end] = B1f[1, :]

Bxy2 = zeros(3, 6)
Bxy2[1, 1:2:end] = B2f[1, :]
Bxy2[2, 2:2:end] = B2f[2, :]
Bxy2[3, 1:2:end] = B2f[2, :]
Bxy2[3, 2:2:end] = B2f[1, :]

t = 1
A1 = det(M1)/2
A2 = det(M2)/2

Ke1 = Bxy1'*D*Bxy1 * A1 * t
Ke2 = Bxy2'*D*Bxy2 * A2 * t


# reshape for x,y
Nxy1 = zeros(Num, 2, 6)
Nxy1[1, 1:2:end] = N1
Nxy1[2, 2:2:end] = N1

Nxy2 = zeros(Num, 2, 6)
Nxy2[1, 1:2:end] = N2
Nxy2[2, 2:2:end] = N2

# traction/pressure pulling forward in x direction
T = [1000;0]
fs_int = Nxy2'*T  # integrand

e1 = Symbolics.sympy_integrate(fs_int[1], y)
e3 = Symbolics.sympy_integrate(fs_int[3], y)
e5 = Symbolics.sympy_integrate(fs_int[5], y)

f1 = substitute(e1, Dict(x=>20, y=>10))  - substitute(e1, Dict(x=>20, y=>0))
f3 = substitute(e3, Dict(x=>20, y=>10))  - substitute(e3, Dict(x=>20, y=>0))
f5 = substitute(e5, Dict(x=>20, y=>10))  - substitute(e5, Dict(x=>20, y=>0))

function assemble(Kg, Ke, connectivity; dofs_per_node=2)
    # --- 1. DOF map ---
    dof_map = [
        (node_id - 1) * dofs_per_node + dof_offset
        for node_id in connectivity
        for dof_offset in 1:dofs_per_node
    ]

    # --- 2. Assembly loops ---
    n = length(dof_map)
    for i in 1:n  # Local row index
        I = dof_map[i]       # Global row index
        for j in 1:n # Local col index
            J = dof_map[j]      # Global col index
            Kg[I, J] += Ke[i, j]
        end
    end
end

Kg = zeros(Float64, 8, 8)
elem1_nodes = [1, 3, 2]
elem2_nodes = [1, 4, 3]
assemble(Kg, Ke1, elem1_nodes)
assemble(Kg, Ke2, elem2_nodes)


fr  = [5000; 0; 5000; 0]
kr= Kg[5:end, 5:end]
dr = kr\fr

d1 = [0;0;dr[1];dr[2];0;0]
d2 = [0;0;dr[1];dr[2];dr[3];dr[4]]

s1 = D*Bxy1*d1
s2 = D*Bxy2*d2

# Output:
# julia> s1 = D*Bxy1*d1
# 3-element Vector{Float64}:
#  1004.8038430744599
#   301.44115292233795
#     2.401921537229731
# 
# julia> s2 = D*Bxy2*d2
# 3-element Vector{Float64}:
#  1103.6257577490567
#   630.847535170994
#    64.85188150520406

