/* Texture functions for cs580 GzLib	*/
#include "stdafx.h"
#include "stdio.h"
#include "Gz.h"

GzColor *image = NULL;
int xs, ys;
int reset = 1;

/* Image texture function */
int tex_fun(float u, float v, float w, GzColor color)
{
	if (w == 99) {
		color[0] = 0.7;
		color[1] = 0.7;
		color[2] = 0.7;
		return 0;
	}

	unsigned char pixel[3];
	unsigned char dummy;
	char foo[8];
	int i, j;
	FILE *fd;

	if (reset)
	{ /* open and load texture file */
		fd = fopen("texture", "rb");
		if (fd == NULL)
		{
			fprintf(stderr, "texture file not found\n");
			exit(-1);
		}
		fscanf(fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
		image = (GzColor *)malloc(sizeof(GzColor) * (xs + 1) * (ys + 1));
		if (image == NULL)
		{
			fprintf(stderr, "malloc for texture image failed\n");
			exit(-1);
		}

		for (i = 0; i < xs * ys; i++)
		{ /* create array of GzColor values */
			fread(pixel, sizeof(pixel), 1, fd);
			image[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
			image[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
			image[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
		}

		reset = 0; /* init is done */
		fclose(fd);
	}

	/* bounds-test u,v to make sure nothing will overflow image array bounds */
	/* determine texture cell corner values and perform bilinear interpolation */
	/* set color to interpolated GzColor value and return */

	u = max(min(1, u), 0);
	v = max(min(1, v), 0);

	float resx = u * (xs - 1);
	float resy = v * (ys - 1);

	int lowerx, upperx, lowery, uppery;

	lowerx = (int)resx;
	upperx = (int)ceil(resx);
	lowery = (int)resy;
	uppery = (int)ceil(resy);

	float s = resx - lowerx;
	float t = resy - lowery;

	GzColor cA, cB, cC, cD;

	for (int i = 0; i < 3; i++) {
		cA[i] = image[xs * lowery + lowerx][i];
		cB[i] = image[xs * lowery + upperx][i];
		cC[i] = image[xs * uppery + upperx][i];
		cD[i] = image[xs * uppery + lowerx][i];
	}

	for (int i = 0; i < 3; i++)
		color[i] = s * t * cC[i] + (1 - s) * t * cD[i] + s * (1 - t) * cB[i] + (1 - s) * (1 - t) * cA[i];
	return GZ_SUCCESS;
}

/* Procedural texture function */
int ptex_fun(float u, float v, float w, GzColor color)
{

	if (w == 99) {
		color[0] = u * (0.7 - 0.1) + 0.1;
		color[1] = u * (0.9 - 0.1) + 0.1;
		color[2] = u * (0.7 - 0.1) + 0.1;
		return 0;
	}
	else {
		color[0] = 0.1;
		color[1] = 0.1;
		color[2] = 0.1;
		return 0;
	}

	return GZ_SUCCESS;
}

/* Free texture memory */
int GzFreeTexture()
{
	if (image != NULL)
		free(image);
	return GZ_SUCCESS;
}
