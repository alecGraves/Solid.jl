
function plotelements(elements; color=:blue, start=0, alpha=0.5, labels=true, fig=nothing, ax=nothing)
	if ax == nothing
		fig=Figure()
		ax = Axis(fig[1, 1], aspect = DataAspect(), xlabel = "x", ylabel = "y")
	end
	xs =  elements[:, :, 1][:]
	ys = elements[:, :, 2][:]
	for i = 1:size(elements, 1)
		e = elements[i, 1:4, :]
		idx = [1, 2, 3, 4, 1]
		exs = e[idx, 1][:]
		eys = e[idx, 2][:]
		lines!(ax, exs, eys, color=color, alpha=alpha)
		if labels
			text!(sum(exs[1:4])/4, sum(eys[1:4])/4, text="$(i+start)", align = (:center, :center), color=color, alpha=alpha)
		end
	end
	scatter!(ax, xs, ys, color=color, alpha=alpha)
	return fig, ax
end


function plotnodes(xs, ys; color=:blue, fig=nothing, ax=nothing)
	if ax == nothing
		fig=Figure()
		ax = Axis(fig[1, 1], aspect = DataAspect(), xlabel = "x", ylabel = "y")
	end

	for i = 1:length(xs)
		text!(xs[i], ys[i], text="$i", color=color)
	end
end

# :rtriangle
function plotpoints(xs, ys; color=:blue, shape=:utriangle, fig=nothing, ax=nothing)
	if ax == nothing
		fig=Figure()
		ax = Axis(fig[1, 1], aspect = DataAspect(), xlabel = "x", ylabel = "y")
	end

	scatter!(ax, xs, ys, color=color, marker=shape, alpha=0.7)
end

function plotforces(fs, xs, ys; color=:blue, scale=1.0, fig=nothing, ax=nothing)
	if ax == nothing
		fig=Figure()
		ax = Axis(fig[1, 1], aspect = DataAspect(), xlabel = "x", ylabel = "y")
	end

	fs *= scale / max(1e-6, maximum(abs.(fs)))
	for i = 1:length(xs)
		if sum(abs.(fs[2*i-1:2*i])) > 0
			arrows2d!(ax, Point2f(xs[i], ys[i]), Vec2f(fs[2*i-1], fs[2*i]), color=color)
		end
	end
end

