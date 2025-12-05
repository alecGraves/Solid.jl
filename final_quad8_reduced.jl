include("./mesh.jl")
include("./plot.jl")
include("./shapes.jl")
include("./integration.jl")
include("./solve.jl")

using GLMakie
update_theme!(theme_dark())

## === 1. mesh ===
# concrete
c1 = trapezoid(
	60., 60., 
	100., 60., 
	100.0-60.0/tan(deg2rad(85.0)), 120.0,
	60.0+60.0/tan(deg2rad(72.0)), 120.0,
	6, 10)
c2 = rectangle(40., 45., 60., 60.,  4, 3)
c3 = rectangle(60., 45., 100., 60.,  6, 3)
c4 = rectangle(100., 45., 120., 60.,  4, 3)

# soil
s1 = rectangle(0., 45., 40., 60.,  5, 3)
s2 = rectangle(120., 45., 160., 60.,  5, 3)
s3 = rectangle(0., 0., 40., 45.,  5, 6)
s4 = rectangle(40., 0., 60.0, 45.,  4, 6)
s5 = rectangle(60., 0., 100.0, 45.,  6, 6)
s6 = rectangle(100., 0., 120.0, 45.,  4, 6)
s7 = rectangle(120., 0., 160.0, 45.,  5, 6)

concrete = vcat(c1, c2, c3, c4)
soil = vcat(s1, s2, s3, s4, s5, s6, s7)
elements = vcat(concrete, soil)

## Plot
fig, ax = plotelements(concrete, color=:gray)
plotelements(soil, color=:brown, start=size(concrete, 1), fig=fig, ax=ax)

# compute node indices
node_xs, node_ys, idx_map = deduplicate(elements[:,:,1][:], elements[:,:,2][:])  # removes overlapping
element_nodes = reshape(idx_map, (size(elements, 1), 4))  # (E x 4) node id's per element
nodes = 1:length(node_xs)  # node id's
# println(element_nodes[1, :])  # [7, 8, 2, 1]
# println(length(node_xs))  # 214
println(xys2node([0.], [0.], node_xs, node_ys)[1])

# add points to convert elements from linear quads (4-node) to cubic quads (8-node)
function addpoints(elements)
	# 4 - 3
	# 1 - 2
	elements2 = zeros(size(elements))
	# loop through each element and generate an intermediate point, adding it to the list
	for j in 1:size(elements, 1)
		# order: bottom, right, top, left
		for (i, f) in [(1, 2), (2, 3), (3, 4), (4, 1)]
			# x and y midpoint
			elements2[j, i, :] = 0.5 * (elements[j, i, :] + elements[j, f, :])
		end
	end

	return elements2
end

elements2 = addpoints(elements)
oct_xs, oct_ys, oct_idx_map = deduplicate(elements2[:, :, 1][:], elements2[:, :, 2][:])
oct_idx_map .+= length(node_xs)
oct_element_nodes = reshape(oct_idx_map, (size(elements2, 1), 4))

node_xs = vcat(node_xs, oct_xs)
node_ys = vcat(node_ys, oct_ys)
element_nodes = cat(element_nodes, oct_element_nodes, dims=2)
nodes = 1:length(node_xs)
println("octquad_elements: $(size(element_nodes))")
println("nodes: $(length(nodes))")

## Plot
scatter!(ax, node_xs, node_ys, color=:orange, markersize=4.0, alpha=0.9)
# plotnodes(node_xs, node_ys, color=:orange, ax=ax)

# === 1.5 Compute constrained nodes ===
constrained = zeros(Int64, length(nodes)*2)
println(length(nodes))
println(nodes[end])
println(size(constrained))

cons_y_s3 = reduce(hcat, lerp2d(0., 0., 40., 0., 9))  # 5 10    (n-1)*2+1
cons_y_s4 = reduce(hcat, lerp2d(40., 0., 60., 0., 7))  # 4 8
cons_y_s5 = reduce(hcat, lerp2d(60., 0., 100., 0., 11))  # 6 12
cons_y_s6 = reduce(hcat, lerp2d(100., 0., 120., 0., 7))  # 4 8
cons_y_s7 = reduce(hcat, lerp2d(120., 0., 160., 0., 9))  # 5 10
cons_y_pts = vcat(cons_y_s3, cons_y_s4, cons_y_s5, cons_y_s6, cons_y_s7)
cons_y_nodes = xys2node(cons_y_pts[:, 1], cons_y_pts[:, 2], node_xs, node_ys)
cons_y_nodes[cons_y_nodes .> length(nodes)] .= length(nodes)
constrained[2 .* cons_y_nodes] .= 1
# scatter!(ax, cons_y_pts[:, 1], cons_y_pts[:, 2], color=:red, marker=:utriangle, alpha=1)
scatter!(ax, node_xs[cons_y_nodes], node_ys[cons_y_nodes], color=:purple, marker=:utriangle, alpha=1)

cons_x_s1 = reduce(hcat, lerp2d(0., 45., 0., 60., 5))  # 3 6
cons_x_s2 = reduce(hcat, lerp2d(160., 45., 160., 60., 5))  # 3 6
cons_x_s3 = reduce(hcat, lerp2d(0., 0., 0., 45., 11))  # 6 12
cons_x_s7 = reduce(hcat, lerp2d(160., 0., 160., 45., 11))  # 6 12
cons_x_pts = vcat(cons_x_s1, cons_x_s2, cons_x_s3, cons_x_s7)
cons_x_nodes = xys2node(cons_x_pts[:, 1], cons_x_pts[:, 2], node_xs, node_ys)
cons_x_nodes[cons_x_nodes .> length(nodes)] .= length(nodes)
constrained[2 .* cons_x_nodes.-1] .= 1

# scatter!(ax, cons_x_pts[:, 1], cons_x_pts[:, 2], color=:blue, marker=:rtriangle, alpha=1)
scatter!(ax, node_xs[cons_x_nodes], node_ys[cons_x_nodes], color=:purple, marker=:rtriangle, alpha=1)


## === 2. Analyze forces ===

# global force vector
ns = 2
Fg = zeros(length(nodes)*ns)

# hydrostatic pressure
c1_xf, c1_yf = lerp2d(60., 60., 60.0+60.0/tan(deg2rad(72.0)), 120.0,  10)
c2_xf, c2_yf = lerp2d(40., 60., 60., 60.,  4)
s1_xf, s1_yf = lerp2d(0., 60., 40., 60.,  5)

c1_nf = xys2node(c1_xf, c1_yf, node_xs, node_ys)
c2_nf = xys2node(c2_xf, c2_yf, node_xs, node_ys)
s1_nf = xys2node(s1_xf, s1_yf, node_xs, node_ys)
hspnodes = sort(unique(vcat(c1_nf, c2_nf, s1_nf)))
# print(hspnodes) [1, 7, 13, 19, 25, 31, 37, 43, 49, 55, 56, 57, 58, 59, 60, 61, 62]

function hsp(x, norm; ρ=1000, g=9.81)
	return ρ*g*(115-x[2]) * ceil((115-x[2])/1e9) * -norm
end

function apply_hsp(Fg, element_nodes, node_xs, node_ys, hspnodes, hspfn)
	for e in eachslice(element_nodes, dims=1)
		fs = zeros(16)
		for (i, f) in [(1, 2); (2, 3); (3, 4); (4, 1)]
			if e[i] ∈ hspnodes && e[f] ∈ hspnodes  # bottom
				xi = [node_xs[e[i]], node_ys[e[i]]]
				xf = [node_xs[e[f]], node_ys[e[f]]]
				Ts = gauss_integrate_line(hspfn, xi, xf, x->shape_isoline_lagrange(x, 3), 3)
				fs[(i*2)-1:i*2] += Ts[1, :]
				fs[((i+4)*2)-1:(i+4)*2] += Ts[2, :]
				fs[(f*2)-1:f*2] += Ts[end, :]
				# println("fs: $fs")
			end
		end
		if sum(abs.(fs)) > 0
			gidx = hcat(2*e.-1, 2*e)'[:]
			Fg[gidx] += fs
		end
	end
end

apply_hsp(Fg, element_nodes, node_xs, node_ys, hspnodes, hsp)
plotforces(Fg, node_xs, node_ys, ax=ax, scale=4.0)

# gravity
cdensity = 2200
sdensity = 1800
function apply_gravity(Fg, element_nodes, node_xs, node_ys, density)
	for en in eachslice(element_nodes, dims=1)
		e = hcat(node_xs[en], node_ys[en])
		fs = gauss_integrate_bodyforce([0; -9.81*density], e, nshape_quad_cub, jacobian_nshape_quad_cub, 3)
		gidx = hcat(2*en.-1, 2*en)'[:]
		Fg[gidx] += fs
	end
end

nc = size(concrete, 1)
Fgrav = zeros(length(nodes)*ns)
apply_gravity(Fgrav, element_nodes[1:nc, :], node_xs, node_ys, cdensity)
apply_gravity(Fgrav, element_nodes[nc+1:end, :], node_xs, node_ys, sdensity)
# plotforces(Fgrav, node_xs, node_ys, ax=ax, scale=2.0, color=:red)
Fg += Fgrav


# constraints

## === 2. Analyze elements ===
Kg = zeros(length(nodes)*2, length(nodes)*2)

function assemble_elements(Kg, element_nodes, node_xs, node_ys, ν, E)
	D = E/(1-ν^2) * [1.0 ν 0; ν 1 0; 0 0 ((1.0-ν)/2.0)]
	for en in eachslice(element_nodes, dims=1)
		e = hcat(node_xs[en], node_ys[en])
		ke = gauss_integrate_element(D, e, nshape_quad_cub, jacobian_nshape_quad_cub, 3)
		gidx = hcat(2*en.-1, 2*en)'[:]
		Kg[gidx, gidx] += ke
	end
end

Ec = 10e9
νc = 0.2
Es = 10e6   # scale_factor=2.59388 for 10% deformation
			# scale_factor=2.590 for 10% deformation - 897 elements, 2848 nodes
			# scale_factor=5.27698 for 10%, no gravity
# Es = 320e6  # scale_factor=84.6210 for 10% deformation
			# scale_factor=167.357 for 10%, no gravity
νs = 0.3
assemble_elements(Kg, element_nodes[1:nc,:], node_xs, node_ys, νc, Ec)
assemble_elements(Kg, element_nodes[nc+1:end,:], node_xs, node_ys, νs, Es)

## === 3. Solve ===
d, fr = partition_solve(Kg, Fg, constrained)

desired_d = 0.1 * maximum(abs, vcat(node_xs, node_ys))
dmax = maximum(abs, d)
scale_factor = max(1,0, desired_d / dmax)
println("scale_factor=$scale_factor")

dx = d[1:2:end]
dy = d[2:2:end]

ndx = node_xs + dx * scale_factor
ndy = node_ys + dy*scale_factor
ed = [hcat(ndx[en], ndy[en]) for en in eachslice(element_nodes, dims=1)]
scatter!(ax, ndx, ndy, color=:purple, marker=:cross, alpha=0.4)
ed = permutedims(cat(ed..., dims=3), (3, 1, 2))



function compute_stress(d, element_nodes, node_xs, node_ys, ν, E)
	D = E/(1-ν^2) * [1.0 ν 0; ν 1 0; 0 0 ((1.0-ν)/2.0)]
	M=3
	xs = zeros(length(element_nodes), M*M)
	ys = zeros(length(element_nodes), M*M)
	vstress = zeros(length(element_nodes), M*M)
	for (i, en) in enumerate(eachslice(element_nodes, dims=1))
		e = hcat(node_xs[en], node_ys[en])
		d_e = reshape(hcat(d[2*en.-1], d[2*en])', (size(element_nodes, 2)*2, 1))
		# println("$(size(d_e)), $d_e")
		exs, eys, stress = compute_stress_intpts(D, d_e, e, nshape_quad_cub, jacobian_nshape_quad_cub, M)
		vonmises = sqrt.(stress[:, 1].^2 .+ (stress[:, 1].*stress[:, 2]) .+ stress[:, 2].^2 + 3. .* stress[:, 3].^2)
		xs[i, :] = exs
		ys[i, :] = eys
		vstress[i, :] = vonmises
	end
	return xs[:], ys[:], vstress[:]
end

xs, ys, vstress = compute_stress(d, element_nodes[1:nc, :], node_xs, node_ys, νc, Ec)

# plotelements(ed, color=:white, fig=fig, ax=ax, alpha=0.9, labels=false)
scatter!(ax, Point2f.(xs, ys), colormap=:buda, color=vstress, markersize=10, alpha=0.9)

scatter!(ax, Point2f(xs[argmax(vstress)], ys[argmax(vstress)]), color=:red, markersize=20, alpha=1.0)

mins = findall(x->(x>0.0 && x < 4*sum(vstress)/length(vstress)), vstress)
vsn = vstress[mins]
xsn = xs[mins]
ysn = ys[mins]


scatter!(ax, Point2f.(xsn, ysn), color=:cyan, markersize=20, alpha=0.9)

# wait(display(fig))

s6c4 = vcat(reduce(hcat, lerp2d(100., 0., 100., 45., (6-1)*2+1)),
	reduce(hcat, lerp2d(100., 45., 100., 60., (3-1)*2+1)))

s6c4_nodes = xys2node(s6c4[:, 1], s6c4[:, 2], node_xs, node_ys)
scatter!(ax, node_xs[s6c4_nodes], node_ys[s6c4_nodes], color=:purple, marker=:rtriangle, alpha=1)

println(dy[s6c4_nodes])
println(node_ys[s6c4_nodes])


using CairoMakie

# --- DATASET 1 ---
disp1 = [0.0, -0.63363945, -1.2351905, -1.8051275, -2.3437295, -2.8512566, -3.327507, -3.7725077, -4.1855927, -4.56711, -4.918353, -4.918353, -4.918902, -4.919468, -4.9203696, -4.921385]
pos1  = [0.0, 4.5, 9.0, 13.5, 18.0, 22.5, 27.0, 31.5, 36.0, 40.5, 45.0, 45.0, 48.75, 52.5, 56.25, 60.0]

# --- DATASET 2 ---
disp2 = Float32[0.0, -0.29227543, -0.57787365, -0.8568405, -1.1292192, -1.3950528, -1.6543789, -1.9072349, -2.1536467, -2.3936427, -2.62723, -2.854421, -3.0751967, -3.2895496, -3.4974327, -3.6988184, -3.893652, -4.081912, -4.2635937, -4.438741, -4.607473, -4.769952, -4.9264, -4.9264, -4.926645, -4.926869, -4.9270887, -4.9273186, -4.9275637, -4.927838, -4.9281354, -4.928509, -4.929101, -4.9297523]
pos2  = [0.0, 2.045454502105713, 4.090909004211426, 6.136363506317139, 8.181818008422852, 10.227272033691406, 12.272727012634277, 14.318181991577148, 16.363636016845703, 18.40909194946289, 20.454545974731445, 22.5, 24.545454025268555, 26.59090805053711, 28.636363983154297, 30.68181800842285, 32.727272033691406, 34.772727966308594, 36.818180084228516, 38.8636360168457, 40.90909194946289, 42.95454406738281, 45.0, 45.0, 46.5, 48.0, 49.5, 51.0, 52.5, 54.0, 55.5, 57.0, 58.5, 60.0]

disp3_strong = [0.0, -0.009308809, -0.018410176, -0.027306734, -0.03600117, -0.04449612, -0.052794326, -0.060898446, -0.06881128, -0.07653555, -0.084074154, -0.09142991, -0.098605946, -0.10560533, -0.11243183, -0.11908921, -0.12558289, -0.1319183, -0.13810457, -0.14415078, -0.15007219, -0.1558821, -0.16160594, -0.16160594, -0.16182062, -0.16202328, -0.16222474, -0.16243605, -0.16265982, -0.1629072, -0.16317254, -0.16349912, -0.1639986, -0.16454564]
pos3 = [0.0, 2.045454502105713, 4.090909004211426, 6.136363506317139, 8.181818008422852, 10.227272033691406, 12.272727012634277, 14.318181991577148, 16.363636016845703, 18.40909194946289, 20.454545974731445, 22.5, 24.545454025268555, 26.59090805053711, 28.636363983154297, 30.68181800842285, 32.727272033691406, 34.772727966308594, 36.818180084228516, 38.8636360168457, 40.90909194946289, 42.95454406738281, 45.0, 45.0, 46.5, 48.0, 49.5, 51.0, 52.5, 54.0, 55.5, 57.0, 58.5, 60.0]

# --- PLOTTING ---
f = Figure(size = (800, 600))
ax = Axis(f[1, 1], 
    title = "Displacement Comparison",
    xlabel = "Vertical Position (y)",
    ylabel = "Displacement (v)",
    xgridstyle = :dash, 
    ygridstyle = :dash
)

# Plot 1 (Blue Circles)
scatterlines!(ax, pos3, disp3_strong, 
    color = :blue, label = "Displacement (Stiff Soil well-graded gravel)", 
    markersize = 10
)

# Plot 2 (Red Triangles, Dashed Line)
scatterlines!(ax, pos2, disp2, 
    color = :red, label = "Displacement (Weak Soil)", 
    linestyle = :dash, marker = :utriangle, markersize = 10
)

axislegend(ax, position = :lb) # Lower Left (since data is negative)
display(f)
