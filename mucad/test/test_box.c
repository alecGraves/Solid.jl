/*  test_box.c
 *
 *  Simple test program that uses the MUCAD API to create a box
 *  and export it as an STL file.
 *
 *  The program will several box model files.
 */

#include "mucad.h"
#include <stdio.h>

int main(void) {
	/* Create a unit box from (0,0,0) to (1,1,1) */
	Shape box = mucad_box(0.0, 0.0, 0.0,
			      1.0, 1.0, 1.0);

	if(!box) {
		fprintf(stderr, "Error: could not create box.\n");
		return 1;
	}

	/* Write the shape to an STL file */
	if(mucad_write_stl(box, "box.stl", -1, 0.1) != 0) {
		fprintf(stderr, "Error: could not write STL file.\n");
		return 1;
	}
	printf("Box written successfully to 'box.stl'.\n");

	if(mucad_write_step(box, "box.step", -1) != 0) {
		fprintf(stderr, "Error: could not write STEP file.\n");
		return 1;
	}
	printf("Box written successfully to 'box.step'.\n");

	if(mucad_write_obj(box, "box.obj", -1, 0.1) != 0) {
		fprintf(stderr, "Error: could not write OBJ file.\n");
		return 1;
	}
	printf("Box written successfully to 'box.obj'.\n");

	if(mucad_write_iges(box, "box.iges", -1) != 0) {
		fprintf(stderr, "Error: could not write IGES file.\n");
		return 1;
	}
	printf("Box written successfully to 'box.iges'.\n");

	/* Clean up */
	// mucad_free(box);

	return 0;
}
