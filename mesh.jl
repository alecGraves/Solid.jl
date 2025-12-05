
function rectangle(xi, yi, xf, yf, nx, ny)
	#           . (xf,yf)
	#  . (xi,yi)
	xs = LinRange(xi, xf, nx)
	ys = LinRange(yf, yi, ny)

	elements = zeros(Float32, (nx-1)*(ny-1), 4, 2)

	for y = 1:(ny-1)
		for x = 1:(nx-1)
			# 4 - 3
			# 1 - 2

			elements[x+((nx-1)*(y-1)), :, :] = [
				xs[x] ys[y+1];
				xs[x+1] ys[y+1];
				xs[x+1] ys[y];
				xs[x] ys[y];
			]
		end
	end

	return elements
end

function lerp2d(xi, yi, xf, yf, n)
	xs = LinRange(xi, xf, n)
	ys = LinRange(yi, yf, n)
	return xs, ys
end

# Assume corners passed in as: bottom left, bottom right, top right, top left
#  4 - 3
#  1 - 2
function trapezoid(x1, y1, x2, y2, x3, y3, x4, y4, nx, ny)
	left_x, left_y = lerp2d(x4, y4, x1, y1, ny)
	right_x, right_y = lerp2d(x3, y3, x2, y2, ny)

	elements = zeros(Float32, (nx-1)*(ny-1), 4, 2)

	for i = 1:(ny-1)
		row_x, row_y = lerp2d(left_x[i], left_y[i], right_x[i], right_y[i], nx)
		row_xn, row_yn = lerp2d(left_x[i+1], left_y[i+1], right_x[i+1], right_y[i+1], nx)
		for x = 1:(nx-1)
			elements[x+((nx-1)*(i-1)), :, :] = [
				row_xn[x] row_yn[x];
				row_xn[x+1] row_yn[x+1];
				row_x[x+1] row_y[x+1]
				row_x[x] row_y[x];
			]
		end
	end
	return elements
end

# merges duplicate points in elements and returns mapping
function deduplicate(xs, ys, rounded=0.1, eps2=1e-6)
	xmin, xmax = extrema(xs)
	idxs = 1:length(xs)

	# sort left to right, top to bottom
	idxs_xsort = sort(idxs, by=i -> floor(xs[i]/rounded) + xmin + (- floor(ys[i]/rounded)) * (xmax-xmin))

	idxs_dd = zeros(Int64, size(idxs_xsort))  # deduplicated
	idxs_dd[1] = idxs_xsort[1]
	idxs_map = zeros(Int64, size(idxs))
	idxs_map[idxs_xsort[1]] = 1
	nd = 1

	eps2 = 1e-6;

	for i in idxs_xsort[2:end]
		if ((xs[i]-xs[idxs_dd[nd]])^2 + (ys[i] - ys[idxs_dd[nd]])^2 > eps2)
			nd += 1
			idxs_dd[nd] = i
		end
		idxs_map[i] = nd
	end
	return xs[idxs_dd[1:nd]], ys[idxs_dd[1:nd]], idxs_map
end

# return node (idx) closest to xy
function xys2node(xs::AbstractVector{<:Real},
                  ys::AbstractVector{<:Real},
                  node_xs::AbstractVector{<:Real},
                  node_ys::AbstractVector{<:Real})
    @assert length(xs) == length(ys) "xs and ys must have the same length"

    n  = length(xs)          # number of query points
    m  = length(node_xs)     # number of nodes
    @assert length(node_ys) == m "node_xs and node_ys must be the same length"

    idxs = Vector{Int}(undef, n)

    @inbounds for i in 1:n
        xi = xs[i]
        yi = ys[i]

        min_d   = typemax(Float64)   # initialize with +∞
        min_j   = 1                   # placeholder - will be updated

        @inbounds for j in 1:m
            dx = node_xs[j] - xi
            dy = node_ys[j] - yi
            d2 = dx*dx + dy*dy        # squared Euclidean distance

            if d2 < min_d
                min_d = d2
                min_j = j
            end
        end

        idxs[i] = min_j
    end

    return idxs
end


# function xys2node(xs, ys, node_xs, node_ys, rounded=0.1)
# 	xmin, xmax = extrema(node_xs)

# 	# we know that the node coordinates are already sorted left to right, top to bottom
# 	toval = (x_, y_) -> floor(x_/rounded) + xmin + (- floor(y_/rounded)) * (xmax-xmin)
# 	node_vals = floor.(node_xs.*(1.0/rounded)) .+ xmin .+ (- floor.(node_ys.*(1.0/rounded))) .* (xmax-xmin)

# 	idxs = zeros(Int64, length(xs))
# 	for i = 1:length(xs)
# 		idxs[i] = first(searchsorted(node_vals, toval(xs[i], ys[i])))
# 	end

# 	return idxs
# end

# function element_nodes(elements, nodes, node_xs, node_ys)
# 	xmin, xmax = extrema(elements[:,:,1])
# 	to_idxval = i -> floor(node_xs[i]/rounded) + xmin + (- floor(node_ys[i]/rounded)) * (xmax-xmin)
# 	idxs_xsort = searchsorted(nodes, by=)
# end

# function trapezoidal_nodes(xs, ys, round=0.1, max_dist=10)
# 	xmin, xmax = minmax(xs)
# 	ymin, ymax = minmax(ys)

# 	idxs = 1:length(xs)

# 	# sort left to right, top to bottom
# 	idxs_xsort = sort(idxs, i -> floor(xs[i]/eps) + floor((ymax - ys[i]) * (xmax-xmin))/eps)

# 	# sort top to bottom, left to right
# 	idxs_ysort = sort(idxs, i -> floor((ymax-ys[i])/eps) + floor(xs[i] * (ymax-ymin))/eps)

# 	elements = zeros(length(xs), 2, 4)
# 	nel = 0

# 	# loop through idxs, building quads
# 	#  1 - 2
# 	#  3 - 4
# 	for idx in idxs_xsort

# 	end
	
# end

