using LinearAlgebra

GAUSS_M = [1, 2, 3, 4, 5]

GAUSS_W = [
	2.0 0 0 0 0;
	1.0 1.0 0 0 0;
	0.5555555555555556 0.8888888888888888 0.5555555555555556 0 9;
	0.3478548451374538 0.6521451548625461 0.6521451548625461 0.3478548451374538 0;
	0.23692688505618908 0.47862867049936647 0.5688888888888889 0.47862867049936647 0.23692688505618908;
]

GAUSS_X = [
	0.0 0 0 0 0;
	-0.5773502691896258 0.5773502691896258 0 0 0;
	-0.7745966692414834 0.0 0.7745966692414834 0 0;
	0.3478548451374538 0.6521451548625461 0.6521451548625461 0.3478548451374538 0;
	-0.906179845938664 -0.5384693101056831 0.0 0.5384693101056831 0.906179845938664;
]

function gauss_x(ξ, xi, xf)
	return (xi + xf) / 2.0 + ξ * (xf - xi) / 2.0
end

function jacobian_1d(xi, xf)
	return norm(xf - xi) / 2.0
end

function gauss_integrate_1d(f, xi, xf, M)
	ws = GAUSS_W[M, 1:M]
	ξs = GAUSS_X[M, 1:M]
	int = zeros(size(xi))

	for (ξ, w) in zip(ws, ξs)
		int += w * f(gauss_x(ξ, xi, xf))
	end

	return int * jacobian_1d(xi, xf)
end

# takes a function f(x, normal), endpoints, isoparametric shape function, and order M
# returns line integral at shape function points.
function gauss_integrate_line(f, xi, xf, shapefn, M)
	d = size(shapefn(0.), 1)

	edge_vec = xf-xi
	det_J = jacobian_1d(xi, xf)
	normal = [edge_vec[2], -edge_vec[1]] / (2*det_J)  # -90 degree 'z' rotation

	ws = GAUSS_W[M, 1:M]
	ξs = GAUSS_X[M, 1:M]
	int = zeros(d, 2)

	for (w, ξ) in zip(ws, ξs)
		val = w * f(gauss_x(ξ, xi, xf), normal) * det_J
		# shapefn(ξ) is (n). val' is (1x2).
		# Result is (n x 2)
		int += shapefn(ξ) * val'
	end
	return int
end

# general integration of 2d volume with arbitrary function
function gauss_integrate_2d(f, e, dshapefn, M)
	ws = GAUSS_W[M, 1:M]
	ξs = GAUSS_X[M, 1:M]

	int = zeros(size(e, 1) * 2)

	for (wξ, ξ) in zip(ws, ξs)
		for (wη, η) in zip(ws, ξs)
			det_J = det(dshapefn(ξ, η) * e)
			int += wξ * wη * f(ξ, η) * det_J
		end
	end
	return int
end

# body force integration
function gauss_integrate_bodyforce(x, e, nodeshape, dshapefn, M)
	d = size(e, 1)  # 4, 8, etc.

	ws = GAUSS_W[M, 1:M]
	ξs = GAUSS_X[M, 1:M]
	int = zeros(d * 2)
	for (wξ, ξ) in zip(ws, ξs)
		for (wη, η) in zip(ws, ξs)
			ns = nodeshape(ξ, η)
			N = zeros(2, d*2)
			N[1, 1:2:end] = ns
			N[2, 2:2:end] = ns

			det_J = det(dshapefn(ξ, η) * e)
			int += wξ * wη * N' * x * det_J
		end
	end
	return int
end

function gauss_integrate_element(D, e, nodeshape, dshapefn, M)
	d = size(e, 1)  # 4, 8, etc.

	B1 = [1. 0 0 0;
		0 0 0 1
		0 1 1 0];

	ws = GAUSS_W[M, 1:M]
	ξs = GAUSS_X[M, 1:M]
	int = zeros(2*d, 2*d)
	for (wξ, ξ) in zip(ws, ξs)
		for (wη, η) in zip(ws, ξs)
			ns = nodeshape(ξ, η)
			N = zeros(2, d*2)
			N[1, 1:2:end] = ns
			N[2, 2:2:end] = ns

			dN = dshapefn(ξ, η)
			J =  dN * e
			det_J = det(J)

			Γ = inv(J)
			B2 = zeros(4, 4)
			B2[1:2, 1:2] = Γ
			B2[3:end, 3:end] = Γ

			B3 = zeros(4, d*2)
			B3[1:2, 1:2:end] = dN
			B3[3:4, 2:2:end] = dN

			B = B1*B2*B3

			int += wξ * wη * B' * D * B *det_J
		end
	end
	return int
end

function compute_stress_intpts(D, d, e, nodeshape, dshapefn, M)
	dim = size(e, 1)  # 4, 8, etc.

	ws = GAUSS_W[M, 1:M]
	ξs = GAUSS_X[M, 1:M]
	B1 = [1. 0 0 0;
		0 0 0 1
		0 1 1 0];

	xs = zeros(M*M)
	ys = zeros(M*M)
	stresses = zeros(M*M, 3)
	i = 1;

	for (wξ, ξ) in zip(ws, ξs)
		for (wη, η) in zip(ws, ξs)
			# Compute stress at integraiont point
			dN = dshapefn(ξ, η)
			J =  dN * e
			Γ = inv(J)
			B2 = zeros(4, 4)
			B2[1:2, 1:2] = Γ
			B2[3:end, 3:end] = Γ
			B3 = zeros(4, dim*2)
			B3[1:2, 1:2:end] = dN
			B3[3:4, 2:2:end] = dN
			B = B1*B2*B3
			stresses[i, :] = D * B * d

			# Compute physical coordinates of integration point
			ns = nodeshape(ξ, η)
			N = zeros(2, dim*2)
			N[1, 1:2:end] = ns
			N[2, 2:2:end] = ns
			xy = N * e'[:];
			xs[i] = xy[1]
			ys[i] = xy[2]
			i += 1
		end
	end
	return xs, ys, stresses
end
