/*  test_ccad_api.c
 *
 *  A self-contained test routine that verifies the core CCAD API.
 *  Compile with:
 *      gcc -o test_ccad_api test_ccad_api.c -lccad -lm
 *
 *  The function does not write any files - it only checks the
 *  geometry produced by the library.
 */

#include "ccad.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

#define EPS 1e-5 /* tolerance for floating-point comparisons */
#define M_PI 3.14159265358979323846264338327950288419
/* ------------------------------------------------------------------ */
/*  Test routine - call from main() or a test harness.                */
/* ------------------------------------------------------------------ */
void test_ccad_api(void) {
	double cx, cy, cz;
	int rc;

	/* -------------------------------------------------------------- */
	/* 1. Box - unit cube */
	Shape box = ccad_box(0.0, 0.0, 0.0,
			      1.0, 1.0, 1.0);
	assert(box != NULL);
	double vol_box = ccad_volume(box);
	assert(fabs(vol_box - 1.0) < EPS);
	rc = ccad_centroid(box, &cx, &cy, &cz);
	assert(rc == 0);
	assert(fabs(cx - 0.5) < EPS);
	assert(fabs(cy - 0.5) < EPS);
	assert(fabs(cz - 0.5) < EPS);
	ccad_free(box);

	/* -------------------------------------------------------------- */
	/* 2. Sphere - radius 2 */
	double r = 2.0;
	Shape sphere = ccad_sphere(0.0, 0.0, 0.0, r);
	assert(sphere != NULL);
	double vol_sphere = ccad_volume(sphere);
	double expected_vol_sphere = (4.0 / 3.0) * M_PI * r * r * r;
	assert(fabs(vol_sphere - expected_vol_sphere) < EPS);
	rc = ccad_centroid(sphere, &cx, &cy, &cz);
	assert(rc == 0);
	assert(fabs(cx) < EPS);
	assert(fabs(cy) < EPS);
	assert(fabs(cz) < EPS);
	ccad_free(sphere);

	/* -------------------------------------------------------------- */
	/* 3. Translate - move a unit box by (3,-2,5) */
	Shape t_box = ccad_box(0.0, 0.0, 0.0,
				1.0, 1.0, 1.0);
	Shape t_box_translated = ccad_translate(t_box, 3.0, -2.0, 5.0);
	assert(t_box_translated != NULL);
	rc = ccad_centroid(t_box_translated, &cx, &cy, &cz);
	assert(rc == 0);
	assert(fabs(cx - (0.5 + 3.0)) < EPS);
	assert(fabs(cy - (0.5 - 2.0)) < EPS);
	assert(fabs(cz - (0.5 + 5.0)) < EPS);
	double vol_t_box = ccad_volume(t_box_translated);
	assert(fabs(vol_t_box - 1.0) < EPS);
	ccad_free(t_box);
	ccad_free(t_box_translated);

	/* -------------------------------------------------------------- */
	/* 4. Rotate - 90° around the Z-axis (about the origin) */
	Shape r_box = ccad_box(1.0, 0.0, 0.0,
				2.0, 1.0, 1.0);
	Shape r_box_rotated = ccad_rotate(r_box,
					   0.0, 0.0, 0.0, /* axis point */
					   0.0, 0.0, 1.0, /* axis direction (z) */
					   M_PI / 2.0);	  /* 90° */
	assert(r_box_rotated != NULL);
	rc = ccad_centroid(r_box_rotated, &cx, &cy, &cz);
	assert(rc == 0);
	double orig_cx = 1.5, orig_cy = 0.5, orig_cz = 0.5;
	double expected_cx = -orig_cy; /* rotation around Z */
	double expected_cy = orig_cx;
	double expected_cz = orig_cz;
	assert(fabs(cx - expected_cx) < EPS);
	assert(fabs(cy - expected_cy) < EPS);
	assert(fabs(cz - expected_cz) < EPS);
	double vol_r_box = ccad_volume(r_box_rotated);
	assert(fabs(vol_r_box - 1.0) < EPS);
	ccad_free(r_box);
	ccad_free(r_box_rotated);

	/* -------------------------------------------------------------- */
	/* 5. Scale - scale a unit box by (2,3,4) about the origin */
	Shape s_box = ccad_box(0.0, 0.0, 0.0,
				1.0, 1.0, 1.0);
	Shape s_box_scaled = ccad_scale(s_box, 2.0, 3.0, 4.0);
	assert(s_box_scaled != NULL);
	rc = ccad_centroid(s_box_scaled, &cx, &cy, &cz);
	assert(rc == 0);
	double expected_cx_s = 0.5 * 2.0;
	double expected_cy_s = 0.5 * 3.0;
	double expected_cz_s = 0.5 * 4.0;
	assert(fabs(cx - expected_cx_s) < EPS);
	assert(fabs(cy - expected_cy_s) < EPS);
	assert(fabs(cz - expected_cz_s) < EPS);
	double vol_s_box = ccad_volume(s_box_scaled);
	double expected_vol_s_box = 1.0 * 2.0 * 3.0 * 4.0;
	assert(fabs(vol_s_box - expected_vol_s_box) < EPS);
	ccad_free(s_box);
	ccad_free(s_box_scaled);

	/* -------------------------------------------------------------- */
	/* 6. Revolve - rectangle (0,0)-(r,h) around the Y-axis -> cylinder */
	double r_cyl = 1.5;
	double h_cyl = 3.0;
	Shape rect = ccad_rectangle(0.0, 0.0, r_cyl, h_cyl);
	Shape cyl = ccad_revolve(rect,
				  0.0, 0.0, 0.0, /* axis point */
				  0.0, 1.0, 0.0, /* axis direction (Y) */
				  2.0 * M_PI);	 /* full revolution */
	assert(cyl != NULL);
	rc = ccad_centroid(cyl, &cx, &cy, &cz);
	assert(rc == 0);
	assert(fabs(cx - 0.0) < EPS);
	assert(fabs(cy - h_cyl / 2.0) < EPS);
	assert(fabs(cz - 0.0) < EPS);
	double vol_cyl = ccad_volume(cyl);
	double expected_vol_cyl = M_PI * r_cyl * r_cyl * h_cyl;
	assert(fabs(vol_cyl - expected_vol_cyl) < EPS);
	ccad_free(rect);
	ccad_free(cyl);

	/* -------------------------------------------------------------- */
	/* 7. Loft - frustum from two circles at z=0 and z=h */
	double r1 = 1.0, r2 = 2.0, h_frustum = 4.0;
	Shape circle1 = ccad_circle(0.0, 0.0, r1);
	Shape circle2 = ccad_circle(0.0, 0.0, r2);
	Shape circle2_z = ccad_translate(circle2, 0.0, 0.0, h_frustum);
	const Shape loft_shapes[2] = {circle1, circle2_z};
	Shape frustum = ccad_loft(loft_shapes, 2, 1, 1);
	assert(frustum != NULL);
	rc = ccad_centroid(frustum, &cx, &cy, &cz);
	assert(rc == 0);
	double vol_frustum = M_PI * h_frustum * (r1 * r1 + r1 * r2 + r2 * r2) / 3.0;
	double expected_cz2 = h_frustum * (r1 * r1 + 2.0 * r1 * r2 + 3.0 * r2 * r2)
			/ (4.0 * (r1 * r1 + r1 * r2 + r2 * r2));
	assert(fabs(cz - expected_cz2) < EPS);
	assert(fabs(cx) < EPS);
	assert(fabs(cy) < EPS);
	double vol_f = ccad_volume(frustum);
	assert(fabs(vol_f - vol_frustum) < EPS);
	ccad_free(circle1);
	ccad_free(circle2);
	ccad_free(circle2_z);
	ccad_free(frustum);

	/* -------------------------------------------------------------- */
	/* 8. Union - two overlapping unit boxes */
	Shape box1 = ccad_box(0.0, 0.0, 0.0,
			       1.0, 1.0, 1.0);
	Shape box2 = ccad_box(0.5, 0.0, 0.0,
			       1.5, 1.0, 1.0);
	Shape union_box = ccad_union(box1, box2);
	assert(union_box != NULL);
	double vol_box1 = ccad_volume(box1);
	double vol_box2 = ccad_volume(box2);
	double vol_union = ccad_volume(union_box);
	/* Union volume must lie between the largest component and the sum */
	assert(vol_union >= (vol_box1 > vol_box2 ? vol_box1 : vol_box2) - EPS);
	assert(vol_union <= vol_box1 + vol_box2 + EPS);
	ccad_free(box1);
	ccad_free(box2);
	ccad_free(union_box);

	/* -------------------------------------------------------------- */
	/* 9. Difference - subtract a smaller box from a larger one */
	Shape big_box = ccad_box(0.0, 0.0, 0.0,
				  2.0, 2.0, 2.0);
	Shape small_box = ccad_box(0.5, 0.5, 0.5,
				    1.5, 1.5, 1.5);
	Shape diff_box = ccad_difference(big_box, small_box);
	assert(diff_box != NULL);
	double vol_big = ccad_volume(big_box);
	double vol_small = ccad_volume(small_box);
	double vol_diff = ccad_volume(diff_box);
	assert(fabs(vol_diff - (vol_big - vol_small)) < EPS);
	ccad_free(big_box);
	ccad_free(small_box);
	ccad_free(diff_box);

	printf("All tests passed.\n");
}

int main(void) {
	test_ccad_api(); /* run the checks */
	return 0;	  /* success */
}
