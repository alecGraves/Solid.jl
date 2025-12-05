
# TODO use static arrays

function shape_isoline_lin(ξ)
	return [0.5 * (1-ξ), 0.5 * (1+ξ)]
end

function shape_isoline_lagrange(ξ, n)
	nodes = range(-1.0, 1.0, length=n)    
	shape = ones(Float64, n)

	for j in 1:n
		for i in 1:n
			if i != j
				# Formula: L_j(ξ) = Π (ξ - ξ_i) / (ξ_j - ξ_i)
				shape[j] *= ((ξ - nodes[i]) / (nodes[j] - nodes[i]))
			end
		end
	end

	return shape
end

function nshape_quad_lin(ξ, η)
	N1 = 0.25 * (1 - ξ) * (1 - η)
	N2 = 0.25 * (1 + ξ) * (1 - η)
	N3 = 0.25 * (1 + ξ) * (1 + η)
	N4 = 0.25 * (1 - ξ) * (1 + η)
	return [N1 N2 N3 N4]
end

# using Symbolics
# @variables x, y
# quad = nshape_quad_lin(x, y)
# Symbolics.jacobian(quad, [x, y])'
function jacobian_nshape_quad_lin(x, y)
	return [-0.25(1 - y)   0.25(1 - y)  0.25(1 + y)  -0.25(1 + y);
 		-0.25(1 - x)  -0.25(1 + x)  0.25(1 + x)   0.25(1 - x)]
end


function nshape_quad_cub(ξ, η)
	N5 = 0.5 * (1 - ξ^2) * (1 - η)
	N6 = 0.5 * (1 + ξ) * (1 - η^2)
	N7 = 0.5 * (1 - ξ^2) * (1 + η)
	N8 = 0.5 * (1 - ξ) * (1 - η^2)
	N1 = 0.25 * (1 - ξ) * (1 - η) - 0.5 * (N8 + N5)
	N2 = 0.25 * (1 + ξ) * (1 - η) - 0.5 * (N5 + N6)
	N3 = 0.25 * (1 + ξ) * (1 + η) - 0.5 * (N6 + N7)
	N4 = 0.25 * (1 - ξ) * (1 + η) - 0.5 * (N7 + N8)
	return [N1 N2 N3 N4 N5 N6 N7 N8];
end

# using Symbolics
# @variables x, y
# quad = nshape_quad_cub(x, y)
# Symbolics.jacobian(quad, [x, y])'
function jacobian_nshape_quad_cub(x, y)
	return [
	 (-0.25(1 - y)-0.5(-x*(1 - y)-0.5(1 - (y^2))))  (-0.25(1 - x)-0.5(-0.5(1 - (x^2)) - (1 - x)*y));
	  (0.25(1 - y)-0.5(-x*(1 - y)+0.5(1 - (y^2)))) ( -0.25(1 + x)-0.5(-0.5(1 - (x^2)) - (1 + x)*y));
	 ( 0.25(1 + y)-0.5(-x*(1 + y)+0.5(1 - (y^2))))    (0.25(1 + x)-0.5(0.5(1 - (x^2)) - (1 + x)*y));
	 (-0.25(1 + y)-0.5(-x*(1 + y)-0.5(1 - (y^2))))   ( 0.25(1 - x)-0.5(0.5(1 - (x^2)) - (1 - x)*y));
					  (-x*(1 - y))                                (-0.5(1 - (x^2)));
				      (0.5(1 - (y^2)))                                     (-(1 + x)*y);
					  (-x*(1 + y))                                ( 0.5(1 - (x^2)));
				     (-0.5(1 - (y^2)))                                     (-(1 - x)*y)]'
end