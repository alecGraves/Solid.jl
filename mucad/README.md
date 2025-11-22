# ccad 

mucad is a small C-ABI-library computer aided design (CAD) package.

It was created to enable me to easily make better CAD files using a parametric constructive solid geometry (CSG) paradigm.

I also hope to use it in FEM and topology optimization.

## Usage

You can use mucad by creating primitives and 
```C
#include "mucad.h"
#include <stdio.h>

int main(void) {
	/* Create a unit box from (0,0,0) to (1,1,1) */
	Shape box = mucad_box(0.0, 0.0, 0.0,
			      1.0, 1.0, 1.0);

	/* Write the shape to a STEP file */
	if(mucad_write_step(box, "box.step", -1, 0.1) != 0) {
		fprintf(stderr, "Error: could not write STEP file.\n");
		return 1;
	}
}
```


## Installation
mucad is a wrapper around Open CASCADE, like CadQuery. But smaller/faster/etc.

## Why "gigacad"?
mucad was taken, as was picocad, femtocad, microcad, megacad, etc.

When I searched for gigacad, I did not find any CAD prorgams, so here we are.

## TODO
- [ ] static lib of Open CASCADE
