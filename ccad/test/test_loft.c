/*  test_loft.c
 *
 *  Simple test program that uses the CCAD API to create a box
 *  and export it as an STL file.
 *
 *  The program will several box model files.
 */
#include "ccad.h"
#include <stdio.h>

int main(void) {
	double r1 = 1.0, h_frustum = 4.0;
	Shape circle = ccad_circle(0.0, 0.0, r1);
	Shape square = ccad_rectangle(0.0, 0.0, 1.0, 1.0);
	Shape square_z = ccad_translate(square, 0.0, 0.0, h_frustum);
	Shape square2_z = ccad_translate(square, 0.0, 0.0, -h_frustum);
	const Shape loft_shapes[3] = {square2_z, circle, square_z};
	Shape loft = ccad_loft(loft_shapes, 3, 0, 0);
	Shape loft_simple = ccad_simplify_to_solid(loft);

	if(!loft_simple) {
		fprintf(stderr, "Error: could not create box.\n");
		return 1;
	}

	if(ccad_write_stl(loft_simple, "loft.step", -1, 0.001) != 0) {
		fprintf(stderr, "Error: could not write STL file.\n");
		return 1;
	}
	printf("Box written successfully to 'loft.stl'.\n");

	/* Clean up */
	// ccad_free(circle);
	// ccad_free(square);
	// ccad_free(square_z);
	// ccad_free(square2_z);
	// ccad_free(loft);

	return 0;
}
