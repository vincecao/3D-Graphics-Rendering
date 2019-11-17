/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"
#include "math.h"

GzColor	*image=NULL;
int xs, ys;
int reset = 1;

/* Image texture function */
int tex_fun(float u, float v, GzColor color)
{
  unsigned char		pixel[3];
  unsigned char     dummy;
  char  		foo[8];
  int   		i, j;
  FILE			*fd;

  if (reset) {          /* open and load texture file */
    fd = fopen ("texture", "rb");
    if (fd == NULL) {
      fprintf (stderr, "texture file not found\n");
      exit(-1);
    }
    fscanf (fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
    image = (GzColor*)malloc(sizeof(GzColor)*(xs+1)*(ys+1));
    if (image == NULL) {
      fprintf (stderr, "malloc for texture image failed\n");
      exit(-1);
    }

    for (i = 0; i < xs*ys; i++) {	/* create array of GzColor values */
      fread(pixel, sizeof(pixel), 1, fd);
      image[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
      image[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
      image[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
      }

    reset = 0;          /* init is done */
	fclose(fd);
  }


/* bounds-test u,v to make sure nothing will overflow image array bounds */
/* determine texture cell corner values and perform bilinear interpolation */
/* set color to interpolated GzColor value and return */
  if (u<0){u=0;}
  else if(u>1){u=1;}
  if (v<0){v=0;}
  else if(v>1){v=1;}

  float s,t;
  int xl,xh,yl,yh;
  int A,B,C,D;
  u=u*(xs-1);
  v=v*(ys-1);
  xl=floor(u);
  xh = ceil(u);

  yl=floor(v);
  yh=ceil(v);  

  s=u-floor(u);
  t=v-floor(v);
  C = yh*xs + xh;
  D = yh*xs + xl;
  A = yl*xs + xl;
  B = yl*xs + xh;


  for (int i=0;i<3;i++){
    color[i]= s*t* image[C][i]+ (1-s)*t *image[D][i]+ s*(1-t)*image[B][i]+ (1-s)*(1-t)*image[A][i]; 
  }
	return GZ_SUCCESS;
}

/* Procedural texture function */
int ptex_fun(float u, float v, GzColor color)
{
  //chess board
	int N = 7;
	bool flag;
	u = u * N;
	v = v * N;
	flag=((int)(round(u))%N + (int)(round(v)) %N) % 2;
	if (flag) {
		color[0] = 0;
		color[1] = 1;
		color[2] = 1;
	}
	else {
		color[0] = 1;
		color[1] = 0;
		color[2] = 1;
	}
  


	return GZ_SUCCESS;
}

/* Free texture memory */
int GzFreeTexture()
{
	if(image!=NULL)
		free(image);
	return GZ_SUCCESS;
}

