
function partition_solve(Kg, fg, constrained, displacements=nothing)
	if displacements == nothing  # optional constrained displacements.
		displacements = zeros(Float32, length(constrained))
	end

	# pre-scale the matrix to improve numerical stability
	# unscale = maximum(abs, Kg.nzval)  # sparse, use nzval to prevent accessing zeros
	unscale = maximum(abs, Kg)
	scale = 1.0/unscale
	K = Kg * scale
	forces = fg * scale

	# add force induced by constrained displacements
	notfree = (constrained .> 0 )
	free = (.!notfree)

	# compute load caused by constraints
	forces_offset = K[:, notfree] * displacements[notfree]

	K_f = K[free, free]
	f_f = forces[free] - forces_offset[free]

	# Solve
	d_f = K_f \ f_f
	displacements[free] .= d_f

	# Compute reactions
	reactions = K * displacements
	reactions[free] .= f_f
	reactions *= unscale

	return displacements, reactions
end

