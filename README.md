# Solid.jl

A Julia CAD and FEM code library for Julia.


There are several key steps to the finite element method:

1. Preprocessing:
	1. Element selection
	2. Meshing / discretization
2. Analysis:
	1. Element formulation, conversion to weak form, numerical integration
	2. Assembly
	3. Solution
3. Post-processing:
	1. Quantify targets (stress, strain, etc.) based on the solution (displacemnet, heat, etc.)
	2. Visualization

so
optimization:

1. Create a static array and use that for matrix construction (200ms -> 10ms)
2. Convert that static Kg into a sparse matrix (cost 10ms)
3. Solve using the sparse matrix (750ms -> 10ms)
4. use static arrays, hunt and eliminate numerous small matrix heap allocations during element matrix construction, force generation, and post-processing.

try gpu parallelizing things too?
- this is useful for larger solves and matrix construction
- under 10k nodes, probably better to just use cpu, or try to parallelize matrix construction.
- Need to test how useful it is for nodes vs solving, construction methods, etc.
- Element matrix and force integration using it would probably be pretty quick.


## TODO
- [ ] 2D elements
- [ ] 3D elements
- [ ] 2D heat
- [ ] 3D heat

Biggest win: 
- automatic mesh creation and alignment, optimizing element selection.

