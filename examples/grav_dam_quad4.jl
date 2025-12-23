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
println("quad_elements: $(size(element_nodes))")
println("nodes: $(length(nodes))")
# println(element_nodes[1, :])  # [7, 8, 2, 1]
# println(length(node_xs))  # 214
println(xys2node([0.], [0.], node_xs, node_ys)[1])

## Plot
plotnodes(node_xs, node_ys, color=:orange, ax=ax)

# === 1.5 Compute constrained nodes ===
constrained = zeros(Int64, length(nodes)*2)

cons_y_s3 = reduce(hcat, lerp2d(0., 0., 40., 0., 5))
cons_y_s4 = reduce(hcat, lerp2d(40., 0., 60., 0., 4))
cons_y_s5 = reduce(hcat, lerp2d(60., 0., 100., 0., 6))
cons_y_s6 = reduce(hcat, lerp2d(100., 0., 120., 0., 4))
cons_y_s7 = reduce(hcat, lerp2d(120., 0., 160., 0., 5))
cons_y_pts = vcat(cons_y_s3, cons_y_s4, cons_y_s5, cons_y_s6, cons_y_s7)
cons_y_nodes = xys2node(cons_y_pts[:, 1], cons_y_pts[:, 2], node_xs, node_ys)
constrained[2 .* cons_y_nodes] .= 1
plotpoints(cons_y_pts[:, 1], cons_y_pts[:, 2], color=:purple, shape=:utriangle, ax=ax)

cons_x_s1 = reduce(hcat, lerp2d(0., 45., 0., 60., 3))
cons_x_s2 = reduce(hcat, lerp2d(160., 45., 160., 60., 3))
cons_x_s3 = reduce(hcat, lerp2d(0., 0., 0., 45., 6))
cons_x_s7 = reduce(hcat, lerp2d(160., 0., 160., 45., 6))
cons_x_pts = vcat(cons_x_s1, cons_x_s2, cons_x_s3, cons_x_s7)
cons_x_nodes = xys2node(cons_x_pts[:, 1], cons_x_pts[:, 2], node_xs, node_ys)
constrained[2 .* cons_x_nodes.-1] .= 1
plotpoints(cons_x_pts[:, 1], cons_x_pts[:, 2], color=:blue, shape=:rtriangle, ax=ax)

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
		fs = zeros(8)
		for (i, f) in [(1, 2); (2, 3); (3, 4); (4, 1)]
			if e[i] ∈ hspnodes && e[f] ∈ hspnodes  # bottom
				xi = [node_xs[e[i]], node_ys[e[i]]]
				xf = [node_xs[e[f]], node_ys[e[f]]]
				Ts = gauss_integrate_line(hspfn, xi, xf, x->shape_isoline_lagrange(x, 2), 2)
				fs[(i*2)-1:i*2] += Ts[1, :]
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
		fs = gauss_integrate_bodyforce([0; -9.81*density], e, nshape_quad_lin, jacobian_nshape_quad_lin, 2)
		gidx = hcat(2*en.-1, 2*en)'[:]
		Fg[gidx] += fs
	end
end

nc = size(concrete, 1)
Fgrav = zeros(length(nodes)*ns)
apply_gravity(Fgrav, element_nodes[1:nc, :], node_xs, node_ys, cdensity)
apply_gravity(Fgrav, element_nodes[nc+1:end, :], node_xs, node_ys, sdensity)
plotforces(Fgrav, node_xs, node_ys, ax=ax, scale=2.0, color=:red)
Fg += Fgrav


# constraints

## === 2. Analyze elements ===
Kg = zeros(length(nodes)*2, length(nodes)*2)

function assemble_elements(Kg, element_nodes, node_xs, node_ys, ν, E)
	D = E/(1-ν^2) * [1.0 ν 0; ν 1 0; 0 0 ((1.0-ν)/2.0)]
	for en in eachslice(element_nodes, dims=1)
		e = hcat(node_xs[en], node_ys[en])
		ke = gauss_integrate_element(D, e, nshape_quad_lin, jacobian_nshape_quad_lin, 2)
		gidx = hcat(2*en.-1, 2*en)'[:]
		Kg[gidx, gidx] += ke
	end
end

Ec = 10e9
νc = 0.2
Es = 10e6
νs = 0.3
assemble_elements(Kg, element_nodes[1:nc,:], node_xs, node_ys, νc, Ec)
assemble_elements(Kg, element_nodes[nc+1:end,:], node_xs, node_ys, νs, Es)

## === 3. Solve ===
d, fr = partition_solve(Kg, Fg, constrained)

desired_d = 0.05 * maximum(abs, vcat(node_xs, node_ys))
dmax = maximum(abs, d)
scale_factor = max(1,0, desired_d / dmax)
println("scale_factor=$scale_factor")

dx = d[1:2:end]
dy = d[2:2:end]

ndx = node_xs + dx * scale_factor
ndy = node_ys + dy*scale_factor
ed = [hcat(ndx[en], ndy[en]) for en in eachslice(element_nodes, dims=1)]
scatter!(ax, ndx, ndy, color=:purple, marker=:cross, alpha=0.3)
ed = permutedims(cat(ed..., dims=3), (3, 1, 2))

plotelements(ed, color=:purple, fig=fig, ax=ax, alpha=0.9, labels=false)

function compute_stress(d, element_nodes, node_xs, node_ys, ν, E)
	D = E/(1-ν^2) * [1.0 ν 0; ν 1 0; 0 0 ((1.0-ν)/2.0)]
	xs = zeros(length(element_nodes), 2*2)
	ys = zeros(length(element_nodes), 2*2)
	vstress = zeros(length(element_nodes), 2*2)
	for (i, en) in enumerate(eachslice(element_nodes, dims=1))
		e = hcat(node_xs[en], node_ys[en])
		d_e = reshape(hcat(d[2*en.-1], d[2*en])', (8, 1))
		exs, eys, stress = compute_stress_intpts(D, d_e, e, nshape_quad_lin, jacobian_nshape_quad_lin, 2)
		vonmises = sqrt.(stress[:, 1].^2 .+ (stress[:, 1].*stress[:, 2]) .+ stress[:, 2].^2 + 3. .* stress[:, 3].^2)
		xs[i, :] = exs
		ys[i, :] = eys
		vstress[i, :] = vonmises
	end
	return xs[:], ys[:], vstress[:]
end

xs, ys, vstress = compute_stress(d, element_nodes[1:nc, :], node_xs, node_ys, νc, Ec)

scatter!(ax, Point2f.(xs, ys), colormap=:buda, color=vstress, markersize=10, alpha=0.6)

wait(display(fig))
