/*  test_loft.c
 *
 *  Simple test program that uses the MUCAD API to create a box
 *  and export it as an STL file.
 *
 *  The program will several box model files.
 */
#include <stdio.h>
#include "mucad.h"

int main(void)
{
    double r1 = 1.0, h_frustum = 4.0;
    Shape circle = mucad_circle(0.0, 0.0, r1);
    Shape square = mucad_rectangle(0.0, 0.0, 1.0, 1.0);
    Shape square_z = mucad_translate(square, 0.0, 0.0, h_frustum);
    Shape square2_z = mucad_translate(square, 0.0, 0.0, -h_frustum);
    const Shape loft_shapes[3] = {square2_z, circle, square_z};
    Shape loft = mucad_loft(loft_shapes, 3, 1, 0);

    if (!loft) {
        fprintf(stderr, "Error: could not create box.\n");
        return 1;
    }

    if (mucad_write_stl(loft, "loft.stl", 0.001) != 0) {
        fprintf(stderr, "Error: could not write STL file.\n");
        return 1;
    }
    printf("Box written successfully to 'loft.stl'.\n");

    /* Clean up */
    // mucad_free(circle);
    // mucad_free(square);
    // mucad_free(square_z);
    // mucad_free(square2_z);
    // mucad_free(loft);

    return 0;
}
