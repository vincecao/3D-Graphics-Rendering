/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
#include <vector>

#define PI (float) 3.14159265358979323846
/***********************************************/
/* HW1 methods: copy here the methods from HW1 */

GzRender::GzRender(int xRes, int yRes)
{
	/* HW1.1 create a framebuffer for MS Windows display:
	 -- set display resolution
	 -- allocate memory for framebuffer : 3 bytes(b, g, r) x width x height
	 -- allocate memory for pixel buffer
	 */

	const int pixelSize = xRes * yRes;
	const int frameSize = pixelSize * 3;
	framebuffer = new char[frameSize];
	pixelbuffer = new GzPixel[pixelSize];
	xres = xRes;
	yres = yRes;
	// init
	m_camera.position[0] = DEFAULT_IM_X;
	m_camera.position[1] = DEFAULT_IM_Y;
	m_camera.position[2] = DEFAULT_IM_Z;
	m_camera.lookat[0] = 0;
	m_camera.lookat[1] = 0;
	m_camera.lookat[2] = 0;
	m_camera.worldup[0] = 0;
	m_camera.worldup[1] = 1;
	m_camera.worldup[2] = 0;
	m_camera.FOV = DEFAULT_FOV;
	matlevel = 0;
	numlights = 0;

}

GzRender::~GzRender()
{
	/* HW1.2 clean up, free buffer memory */
	delete[] framebuffer;
	delete[] pixelbuffer;

}

int GzRender::GzDefault()
{
	/* HW1.3 set pixel buffer to some default values - start a new frame */
	int loc;
	for (int i = 0; i < xres; i++) {
		for (int j = 0; j < yres; j++) {
			loc = ARRAY(i, j);
			// 32767
			GzPut(i, j, 128 << 4, 128 << 4, 128 << 4, 1, INT_MAX);//r,g,b,a,z(max of the short)

		}
	}

	return GZ_SUCCESS;
}


int GzRender::GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
	/* HW1.4 write pixel values into the buffer */

	if (i >= 0 && j >= 0 && i < xres && j < yres) {
		// boundary condition
		if (b < 0) b = 0;
		if (g < 0) g = 0;
		if (r < 0) r = 0;
		if (b > 4095) b = 4095;
		if (g > 4095) g = 4095;
		if (r > 4095) r = 4095;

		// put value into the target pixel
		int loc = ARRAY(i, j);
		// what is the order for the value in framegbuffer?
		pixelbuffer[loc].red = r;
		pixelbuffer[loc].green = g;
		pixelbuffer[loc].blue = b;
		pixelbuffer[loc].alpha = a;
		pixelbuffer[loc].z = z;

	}
	return GZ_SUCCESS;
}


int GzRender::GzGet(int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z)
{
	/* HW1.5 retrieve a pixel information from the pixel buffer */
	if (i < xres && i >= 0 && j < yres &&j >= 0) {
		int loc = ARRAY(i, j);
		*b = pixelbuffer[loc].blue;
		*g = pixelbuffer[loc].green;
		*r = pixelbuffer[loc].red;
		*a = pixelbuffer[loc].alpha;
		*z = pixelbuffer[loc].z;
	}

	return GZ_SUCCESS;
}


int GzRender::GzFlushDisplay2File(FILE* outfile)
{
	/* HW1.6 write image to ppm file -- "P6 %d %d 255\r" */
		// from pixel buffer to PPM

	(void)fprintf(outfile, "P6\n%d %d\n255\n", xres, yres);
	int loc;
	char r, g, b;
	for (int i = 0; i < xres; i++) {
		for (int j = 0; j < yres; j++) {
			loc = ARRAY(j, i);
			r = pixelbuffer[loc].red >> 4;
			g = pixelbuffer[loc].green >> 4;
			b = pixelbuffer[loc].blue >> 4;
			(void)fprintf(outfile, "%c%c%c", r, g, b);

		}
	}

	//the outfile closes at main script
	return GZ_SUCCESS;


}

int GzRender::GzFlushDisplay2FrameBuffer()
{
	/* HW1.7 write pixels to framebuffer:
		- put the pixels into the frame buffer
		- CAUTION: when storing the pixels into the frame buffer, the order is blue, green, and red
		- NOT red, green, and blue !!!
	*/
	int loc;
	for (int i = 0; i < xres; i++) {
		for (int j = 0; j < yres; j++) {
			loc = ARRAY(i, j);
			framebuffer[loc * 3] = pixelbuffer[loc].blue >> 4;
			framebuffer[loc * 3 + 1] = pixelbuffer[loc].green >> 4;
			framebuffer[loc * 3 + 2] = pixelbuffer[loc].red >> 4;

		}
	}

	return GZ_SUCCESS;
}


/***********************************************/
/* HW2 methods: implement from here */

int GzRender::GzPutAttribute(int numAttributes, GzToken	*nameList, GzPointer *valueList)
{
	/* HW 2.1
	-- Set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
	-- In later homeworks set shaders, interpolaters, texture maps, and lights
	*/

	GzPointer value;
	for (int i = 0; i < numAttributes; i++) {
		value = valueList[i];
		if (nameList[i] == GZ_RGB_COLOR) {
			GzColor* p;
			p = (GzColor*)value;
			for (int i = 0; i < 3; i++) {
				flatcolor[i] = (*p)[i];
			}
		}
		if (nameList[i] == GZ_INTERPOLATE) {
			int* p;
			p = (int*)value;
			interp_mode = *p;
		}
		if (nameList[i] == GZ_DIRECTIONAL_LIGHT) {
			GzLight* p;
			p = (GzLight*)value;
			for (int i = 0; i < 3; i++) {
				lights[numlights].direction[i] = (*p).direction[i];
				lights[numlights].color[i] = (*p).color[i];
			}
			numlights++;
		}
		if (nameList[i] == GZ_AMBIENT_LIGHT) {
			GzLight* p;
			p = (GzLight*)value;
			for (int i = 0; i < 3; i++) {
				ambientlight.direction[i] = (*p).direction[i];
				ambientlight.color[i] = (*p).color[i];
			}
		}
		if (nameList[i] == GZ_AMBIENT_COEFFICIENT) {
			GzColor* p;
			p = (GzColor*)value;
			for (int i = 0; i < 3; i++) {
				Ka[i] = (*p)[i];
			}
		}
		if (nameList[i] == GZ_SPECULAR_COEFFICIENT) {
			GzColor* p;
			p = (GzColor*)value;
			for (int i = 0; i < 3; i++) {
				Ks[i] = (*p)[i];
			}
		}

		if (nameList[i] == GZ_DIFFUSE_COEFFICIENT) {
			GzColor* p;
			p = (GzColor*)value;
			for (int i = 0; i < 3; i++) {
				Kd[i] = (*p)[i];
			}
		}
		if (nameList[i] == GZ_DISTRIBUTION_COEFFICIENT) {
			float* p;
			p = (float*)value;
			spec = *p;
		}
		if (nameList[i] == GZ_TEXTURE_MAP) {
			tex_fun = (GzTexture)value;
		}
	}
	return GZ_SUCCESS;
}

float CalSlope(float x1, float x2, float y1, float y2, bool* flag) {
	float slope;
	*flag = false;
	if (y1 == y2) {
		*flag = true;
		slope = D3D10_FLOAT32_MAX;
	}
	else {
		slope = (x1 - x2) / (y1 - y2);
	}
	return slope;
}

float CalLoc(float start_y, float slope, float start_x, float end_x, bool flag) {
	float end_y;
	if (flag) {
		end_y = start_y;
	}
	else {
		end_y = start_y + slope * (end_x - start_x);
	}
	return end_y;
}

void convertCoordinates(GzCoord* x, GzMatrix matrix) {
	float coord[4], changed[4];
	for (int i = 0; i < 3; i++) {
		coord[i] = (*x)[i];
	}
	coord[3] = 1;
	// vector multiplication
	for (int i = 0; i < 4; i++) {
		changed[i] = 0;
		for (int j = 0; j < 4; j++) {
			changed[i] += matrix[i][j] * coord[j];
		}
	}
	for (int i = 0; i < 3; i++) {
		(*x)[i] = changed[i] / changed[3];
	}

}
void vecProduct(GzCoord x, GzCoord y, GzCoord result) {
	for (int i = 0; i < 3; i++) {
		result[i] = x[i] * y[i];
	}

}

void crossProduct(GzCoord y, GzCoord z, GzCoord* result) {

	(*result)[0] = y[1] * z[2] - y[2] * z[1];
	(*result)[1] = y[2] * z[0] - y[0] * z[2];
	(*result)[2] = y[0] * z[1] - y[1] * z[0];

}
float dotproduct(GzCoord x, GzCoord y) { return (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]); }
float norm2(GzCoord x) {
	float var;
	var = pow(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2), 0.5);
	if (var == 0) {
		var = 1;
	}
	return var;
}
int signcheck(float var) {
	if (var > 0) {
		return 1;
	}
	else if (var < 0) {
		return -1;
	}
	else {
		return 0;
	}
}

void computeColor(GzCoord N, GzLight* lights, int numlights, GzColor La, GzColor Ka, GzColor Ks, GzColor Kd, float spec, GzColor* C) {
	// image space E(0,0,-1)
	// R = 2*dotproduct(N,L)*N - L
	// computeColor(*n1,KaLa,lights,numlights,Ks,Kd,spec,&c1);  

	GzColor R, E, sum_s, sum_d, KaLa;
	float temp, dot, RE, NL;
	int flag = 0;

	for (int i = 0; i < 3; i++) {
		KaLa[i] = Ka[i] * La[i];
	}
	//norm N
	temp = norm2(N);
	for (int i = 0; i < 3; i++) {
		N[i] = N[i] / temp;
	}
	for (int i = 0; i < 3; i++) {
		E[i] = 0;
		sum_s[i] = 0;
		sum_d[i] = 0;
	}
	E[2] = -1;


	for (int l = 0; l < numlights; l++) {
		// cal NL
		flag = 0;
		temp = dotproduct(N, E);
		flag += signcheck(temp);
		temp = dotproduct(N, lights[l].direction);
		flag += signcheck(temp);
		for (int i = 0; i < 3; ++i)
		{
			R[i] = 2 * temp*N[i] - lights[l].direction[i];
		}
		if (flag == 0) {
			// skip it
			NL = 0;
		}
		else if (flag == 2) {
			// cal
			NL = temp;
		}
		else if (flag == -2) {
			// FLIP
			NL = -temp;
		}

		RE = dotproduct(R, E);
		if (RE > 1) { RE = 1; }
		if (RE < 0) { RE = 0; }
		RE = pow(RE, spec);
		for (int i = 0; i < 3; i++) {
			sum_d[i] += lights[l].color[i] * NL;
			sum_s[i] += lights[l].color[i] * RE;
		}
	}

	for (int i = 0; i < 3; i++) {
		(*C)[i] = Ks[i] * sum_s[i] + Kd[i] * sum_d[i] + KaLa[i];
	}
}
void computeColorPart(GzCoord N, GzLight* lights, int numlights, GzColor La, float spec, GzColor* C) {
	// image space E(0,0,-1)
	// R = 2*dotproduct(N,L)*N - L

	GzColor R, E, sum_s, sum_d, KaLa;
	float temp, dot, RE, NL;
	int flag = 0;

	//norm N
	temp = norm2(N);
	for (int i = 0; i < 3; i++) {
		N[i] = N[i] / temp;
	}

	for (int i = 0; i < 3; i++) {
		E[i] = 0;
		sum_s[i] = 0;
		sum_d[i] = 0;
	}
	E[2] = -1;


	for (int l = 0; l < numlights; l++) {
		// cal NL
		flag = 0;
		temp = dotproduct(N, E);
		flag += signcheck(temp);
		temp = dotproduct(N, lights[l].direction);
		flag += signcheck(temp);
		//RE
		for (int i = 0; i < 3; ++i)
		{
			R[i] = 2 * temp*N[i] - lights[l].direction[i];
		}
		RE = pow(dotproduct(R, E), spec);

		if (flag == 0) {
			// skip it
			NL = 0;
		}
		else if (flag == 2) {
			// cal
			NL = temp;

		}
		else if (flag == -2) {
			// FLIP
			NL = -temp;
		}

		RE = dotproduct(R, E);
		if (RE > 1) { RE = 1; }
		if (RE < 0) { RE = 0; }
		RE = pow(RE, spec);
		for (int i = 0; i < 3; i++) {
			sum_d[i] += lights[l].color[i] * NL;
			sum_s[i] += lights[l].color[i] * RE;
		}
	}
	for (int i = 0; i < 3; i++) {
		(*C)[i] = sum_s[i] + sum_d[i] + La[i];
	}
}
void calPlaneParams(GzCoord v1, GzCoord v2, GzCoord v3, double* A, double* B, double* C, double* D) {
	// calPlaneParams(p1,p2,p3,&A1,&B1,&C1,&D1);
	GzCoord norm;
	GzCoord e1, e2;
	for (int i = 0; i < 3; i++) {
		e1[i] = v1[i] - v2[i];
		e2[i] = v2[i] - v3[i];
	}
	crossProduct(e1, e2, &norm);
	*A = norm[0];
	*B = norm[1];
	*C = norm[2];
	*D = -norm[0] * v1[0] - norm[1] * v1[1] - norm[2] * v1[2];
}
float calZ(float x, float y, double A, double B, double C, double D) {
	return (-D - A * x - B * y) / C;
}
int GzRender::GzPutTriangle(int	numParts, GzToken *nameList, GzPointer *valueList)
/* numParts - how many names and values */
{
	/* HW 2.2
	-- Pass in a triangle description with tokens and values corresponding to
		  GZ_NULL_TOKEN:		do nothing - no values
		  GZ_POSITION:		3 vert positions in model space
	-- Invoke the rastrizer/scanline framework
	-- Return error code
	*/

	typedef	struct {
		float dx_dy;
		float dz_dy;
		bool fxy;
		bool fzy;
	} EdgeSlope;
	typedef GzCoord   Coord3[3];
	typedef GzTextureIndex   Texture3[3];
	bool fzx = false;
	int sorted_y[3] = { 0,1,2 };
	float deltaX, deltaY, slope_dz_dx;
	float start_x, start_z, end_z, end_x, y, x;
	EdgeSlope slope12, slope13, slope23;//12,13,23, dx/dy,dz/dy
	int temp;
	int first_vertix;
	double var, z_var, z_var1, z_var2, z_var3;
	GzCoord *v1, *v2, *v3;
	GzCoord *n1, *n2, *n3;
	GzTextureIndex *t1, *t2, *t3;
	GzTextureIndex t;
	GzColor c1, c2, c3;
	GzColor C, Ct;
	GzIntensity r, g, b, a;
	GzDepth z_pre, z;
	double A1, B1, C1, D1;
	double A2, B2, C2, D2;
	double A3, B3, C3, D3;
	double Au, Bu, Cu, Du;
	double Av, Bv, Cv, Dv;
	double Az, Bz, Cz, Dz;
	GzCoord p1, p2, p3, N;
	Coord3* loc;
	Texture3* tex;


	loc = (Coord3*)*valueList;
	//interp_mode = GZ_COLOR;




	for (int i = 0; i < 3; i++) {
		convertCoordinates((GzCoord*)(*loc)[i], Ximage[matlevel - 1]);
	}


	// last questiona:check vertice


	// step 1 check the edge relationship find top y
	for (int i = 0; i < 3 - 1; i++) {
		for (int j = 0; j < 3 - 1 - i; j++) {
			if ((*loc)[sorted_y[j]][Y] > (*loc)[sorted_y[j + 1]][Y]) {
				temp = sorted_y[j];
				sorted_y[j] = sorted_y[j + 1];
				sorted_y[j + 1] = temp;
			}
		}
	}

	for (int i = 0; i < 3; i++) {
		if (sorted_y[i] == 0) {
			first_vertix = i;
			break;
		}
	}

	v1 = (*loc + sorted_y[0]);
	v2 = (*loc + sorted_y[1]);
	v3 = (*loc + sorted_y[2]);

	// skip if any z value >0 => calculate
	if ((*v1)[Z] >= 0 || (*v2)[Z] >= 0 || (*v3)[Z] >= 0) {

		loc = (Coord3*)valueList[1];// the loc of norm

		n1 = (*loc + sorted_y[0]);
		n2 = (*loc + sorted_y[1]);
		n3 = (*loc + sorted_y[2]);

		// convert to image space
		convertCoordinates(n1, Xnorm[matlevel - 1]);
		convertCoordinates(n2, Xnorm[matlevel - 1]);
		convertCoordinates(n3, Xnorm[matlevel - 1]);

		tex = (Texture3*)valueList[2];
		t1 = (*tex + sorted_y[0]);
		t2 = (*tex + sorted_y[1]);
		t3 = (*tex + sorted_y[2]);



		// calculate the light 
		for (int i = 0; i < 3; i++) {
			p1[i] = (*v1)[i];
			p2[i] = (*v2)[i];
			p3[i] = (*v3)[i];
		}
		calPlaneParams(p1, p2, p3, &Az, &Bz, &Cz, &Dz);

		// only for normal
		if (interp_mode == GZ_NORMALS) {
			p1[2] = (*n1)[0];
			p2[2] = (*n2)[0];
			p3[2] = (*n3)[0];
			calPlaneParams(p1, p2, p3, &A1, &B1, &C1, &D1);
			p1[2] = (*n1)[1];
			p2[2] = (*n2)[1];
			p3[2] = (*n3)[1];
			calPlaneParams(p1, p2, p3, &A2, &B2, &C2, &D2);
			p1[2] = (*n1)[2];
			p2[2] = (*n2)[2];
			p3[2] = (*n3)[2];
			calPlaneParams(p1, p2, p3, &A3, &B3, &C3, &D3);
		}
		else if (interp_mode == GZ_COLOR)
		{
			computeColorPart(*n1, lights, numlights, ambientlight.color, spec, &c1);
			computeColorPart(*n2, lights, numlights, ambientlight.color, spec, &c2);
			computeColorPart(*n3, lights, numlights, ambientlight.color, spec, &c3);
			p1[2] = c1[0];
			p2[2] = c2[0];
			p3[2] = c3[0];
			calPlaneParams(p1, p2, p3, &A1, &B1, &C1, &D1);
			p1[2] = c1[1];
			p2[2] = c2[1];
			p3[2] = c3[1];
			calPlaneParams(p1, p2, p3, &A2, &B2, &C2, &D2);
			p1[2] = c1[2];
			p2[2] = c2[2];
			p3[2] = c3[2];
			calPlaneParams(p1, p2, p3, &A3, &B3, &C3, &D3);
		}
		// get u,v parameter
		z_var1 = ((double)(*v1)[Z]) / (double)(INT_MAX - (*v1)[Z]);
		z_var2 = ((double)(*v2)[Z]) / (double)(INT_MAX - (*v2)[Z]);
		z_var3 = ((double)(*v3)[Z]) / (double)(INT_MAX - (*v3)[Z]);

		/*p1[2] = (*t1)[U];
		p2[2] = (*t2)[U] ;
		p3[2] = (*t3)[U] ;*/

		p1[2] = ((double)(*t1)[U]) / (z_var1 + 1);
		p2[2] = ((double)(*t2)[U]) / (z_var2 + 1);
		p3[2] = ((double)(*t3)[U]) / (z_var3 + 1);
		calPlaneParams(p1, p2, p3, &Au, &Bu, &Cu, &Du);


		/*p1[2] = (*t1)[V];
		p2[2] = (*t2)[V];
		p3[2] = (*t3)[V];*/

		p1[2] = ((double)(*t1)[V]) / (z_var1 + 1);
		p2[2] = ((double)(*t2)[V]) / (z_var2 + 1);
		p3[2] = ((double)(*t3)[V]) / (z_var3 + 1);
		calPlaneParams(p1, p2, p3, &Av, &Bv, &Cv, &Dv);

		// skip if any z value >0 => calculate
				////		// step 3 computing slope
		slope12.dx_dy = CalSlope((*v1)[X], (*v2)[X], (*v1)[Y], (*v2)[Y], &slope12.fxy);
		slope12.dz_dy = CalSlope((*v1)[Z], (*v2)[Z], (*v1)[Y], (*v2)[Y], &slope12.fzy);
		slope13.dx_dy = CalSlope((*v1)[X], (*v3)[X], (*v1)[Y], (*v3)[Y], &slope13.fxy);
		slope13.dz_dy = CalSlope((*v1)[Z], (*v3)[Z], (*v1)[Y], (*v3)[Y], &slope13.fzy);
		slope23.dx_dy = CalSlope((*v3)[X], (*v2)[X], (*v3)[Y], (*v2)[Y], &slope23.fxy);
		slope23.dz_dy = CalSlope((*v3)[Z], (*v2)[Z], (*v3)[Y], (*v2)[Y], &slope23.fzy);


		// step 4 fill in value
		// part 1
		for (int y = ceil((*v1)[Y]); y < (*v2)[Y]; y++) {

			start_x = CalLoc((*v1)[X], slope12.dx_dy, (*v1)[Y], y, slope12.fxy);
			start_z = CalLoc((*v1)[Z], slope12.dz_dy, (*v1)[Y], y, slope12.fzy);

			end_x = CalLoc((*v1)[X], slope13.dx_dy, (*v1)[Y], y, slope13.fxy);
			end_z = CalLoc((*v1)[Z], slope13.dz_dy, (*v1)[Y], y, slope13.fzy);
			if (end_x < start_x) {
				var = end_x;
				end_x = start_x;
				start_x = var;

				var = end_z;
				end_z = start_z;
				start_z = var;
			}
			slope_dz_dx = CalSlope(end_z, start_z, end_x, start_x, &fzx);

			for (int x = ceil(start_x); x < end_x; x++) {

				z = calZ(x, y, Az, Bz, Cz, Dz);
				//z = CalLoc(start_z, slope_dz_dx, start_x, x,fzx);
					// check the state
				GzGet(x, y, &r, &g, &b, &a, &z_pre);

				if (z_pre > z && z >= 0) {
					//  change u,v to the projective domain
					t[U] = calZ(x, y, Au, Bu, Cu, Du);
					t[V] = calZ(x, y, Av, Bv, Cv, Dv);
					z_var = ((double)z) / (double)(INT_MAX - z);
					t[U] = ((double)t[U])*(z_var + 1.0);
					t[V] = ((double)t[V])*(z_var + 1.0);
					(*tex_fun)(t[U], t[V], Ct);


					if (interp_mode == GZ_FLAT) {
						GzPut(x, y, ctoi(flatcolor[RED]), ctoi(flatcolor[GREEN]), ctoi(flatcolor[BLUE]), 1, z);
					}
					else if (interp_mode == GZ_COLOR) {
						C[0] = calZ(x, y, A1, B1, C1, D1);
						C[1] = calZ(x, y, A2, B2, C2, D2);
						C[2] = calZ(x, y, A3, B3, C3, D3);

						for (int i = 0; i < 3; i++) {
							C[i] = Ct[i] * C[i];
							if (C[i] > 1)C[i] = 1;
							if (C[i] < 0)C[i] = 0;
						}
						GzPut(x, y, ctoi(C[RED]), ctoi(C[GREEN]), ctoi(C[BLUE]), 1, z);
					}
					else if (interp_mode == GZ_NORMALS) {
						// the domain might have some problem!!
						N[0] = calZ(x, y, A1, B1, C1, D1);
						N[1] = calZ(x, y, A2, B2, C2, D2);
						N[2] = calZ(x, y, A3, B3, C3, D3);
						//computeColor(N, lights, numlights, ambientlight.color, Ka, Ks, Kd, spec, &C);
						computeColor(N, lights, numlights, ambientlight.color, Ct, Ks, Ct, spec, &C);
						for (int i = 0; i < 3; i++) {
							if (C[i] > 1)C[i] = 1;
							if (C[i] < 0)C[i] = 0;
						}
						GzPut(x, y, ctoi(C[RED]), ctoi(C[GREEN]), ctoi(C[BLUE]), 1, z);
					}
				}
			}
		}

		// part2
		for (int y = ceil((*v2)[Y]); y < (*v3)[Y]; y++) {
			start_x = CalLoc((*v2)[X], slope23.dx_dy, (*v2)[Y], y, slope23.fxy);
			start_z = CalLoc((*v2)[Z], slope23.dz_dy, (*v2)[Y], y, slope23.fzy);
			// not right
			end_x = CalLoc((*v1)[X], slope13.dx_dy, (*v1)[Y], y, slope13.fxy);
			end_z = CalLoc((*v1)[Z], slope13.dz_dy, (*v1)[Y], y, slope13.fzy);

			if (end_x < start_x) {
				var = end_x;
				end_x = start_x;
				start_x = var;

				var = end_z;
				end_z = start_z;
				start_z = var;
			}

			slope_dz_dx = CalSlope(end_z, start_z, end_x, start_x, &fzx);

			for (int x = ceil(start_x); x < end_x; x++) {
				z = calZ(x, y, Az, Bz, Cz, Dz);
				//z= CalLoc(start_z, slope_dz_dx, start_x, x,fzx);
					// check the state

				GzGet(x, y, &r, &g, &b, &a, &z_pre);

				if (z_pre > z && z >= 0) {
					// todo copy paste the previous method
					t[U] = calZ(x, y, Au, Bu, Cu, Du);
					t[V] = calZ(x, y, Av, Bv, Cv, Dv);
					z_var = ((double)z) / ((double)(INT_MAX - z));
					t[U] = t[U] * (z_var + 1.0);
					t[V] = t[V] * (z_var + 1.0);
					(*tex_fun)(t[U], t[V], Ct);
					if (interp_mode == GZ_FLAT) {
						GzPut(x, y, ctoi(flatcolor[RED]), ctoi(flatcolor[GREEN]), ctoi(flatcolor[BLUE]), 1, z);
					}
					else if (interp_mode == GZ_COLOR) {
						C[0] = calZ(x, y, A1, B1, C1, D1);
						C[1] = calZ(x, y, A2, B2, C2, D2);
						C[2] = calZ(x, y, A3, B3, C3, D3);

						for (int i = 0; i < 3; i++) {
							C[i] = Ct[i] * C[i];
							if (C[i] > 1)C[i] = 1;
							if (C[i] < 0)C[i] = 0;
						}
						GzPut(x, y, ctoi(C[RED]), ctoi(C[GREEN]), ctoi(C[BLUE]), 1, z);
					}
					else if (interp_mode == GZ_NORMALS) {
						// the domain might have some problem!!
						N[0] = calZ(x, y, A1, B1, C1, D1);
						N[1] = calZ(x, y, A2, B2, C2, D2);
						N[2] = calZ(x, y, A3, B3, C3, D3);
						//computeColor(N, lights, numlights, ambientlight.color, Ka, Ks, Kd, spec, &C);
						computeColor(N, lights, numlights, ambientlight.color, Ct, Ks, Ct, spec, &C);
						for (int i = 0; i < 3; i++) {
							if (C[i] > 1)C[i] = 1;
							if (C[i] < 0)C[i] = 0;
						}
						GzPut(x, y, ctoi(C[RED]), ctoi(C[GREEN]), ctoi(C[BLUE]), 1, z);
					}
				}
			}

		}

	}
	return GZ_SUCCESS;
}

float** GzRender::GzAddGrass(int numParts, GzToken *nameList, GzPointer *valueList)
{
	typedef	struct {
		float dx_dy;
		float dz_dy;
		bool fxy;
		bool fzy;
	} EdgeSlope;
	typedef GzCoord   Coord3[3];
	typedef GzTextureIndex   Texture3[3];
	bool fzx = false;
	int sorted_y[3] = { 0,1,2 };
	float deltaX, deltaY, slope_dz_dx;
	float start_x, start_z, end_z, end_x, y, x;
	EdgeSlope slope12, slope13, slope23;//12,13,23, dx/dy,dz/dy
	int temp;
	int first_vertix;
	double var, z_var, z_var1, z_var2, z_var3;
	GzCoord *v1, *v2, *v3;
	GzCoord top, base1, base2, normal;
	GzCoord *n1, *n2, *n3;
	GzTextureIndex t1, t2, t3;
	GzTextureIndex t;
	GzColor c1, c2, c3;
	GzColor C, Ct;
	GzIntensity r, g, b, a;
	GzDepth z_pre, z;
	double A1, B1, C1, D1;
	double A2, B2, C2, D2;
	double A3, B3, C3, D3;
	double Au, Bu, Cu, Du;
	double Av, Bv, Cv, Dv;
	double Az, Bz, Cz, Dz;
	GzCoord p1, p2, p3, N;
	Coord3* loc;
	Texture3* tex;
	float heigth;
	float width;

	float **result = (float **)malloc(3 * sizeof(float *));
	for (int i = 0; i < 3; i++)
		result[i] = (float *)malloc(8 * sizeof(float));

	/*float *resultOne = new float[8];
	float *resultTwo = new float[8];
	float *resultThree = new float[8];*/

	loc = (Coord3*)*valueList;

	// step 1 check the edge relationship find top y
	for (int i = 0; i < 3 - 1; i++) {
		for (int j = 0; j < 3 - 1 - i; j++) {
			if ((*loc)[sorted_y[j]][Y] > (*loc)[sorted_y[j + 1]][Y]) {
				temp = sorted_y[j];
				sorted_y[j] = sorted_y[j + 1];
				sorted_y[j + 1] = temp;
			}
		}
	}

	v1 = (*loc + sorted_y[0]);
	v2 = (*loc + sorted_y[1]);
	v3 = (*loc + sorted_y[2]);

	// the loc of norm
	loc = (Coord3*)valueList[1];

	n1 = (*loc + sorted_y[0]);
	n2 = (*loc + sorted_y[1]);
	n3 = (*loc + sorted_y[2]);

	// calculate the light 
	for (int i = 0; i < 3; i++) {
		p1[i] = (*v1)[i];
		p2[i] = (*v2)[i];
		p3[i] = (*v3)[i];
	}
	calPlaneParams(p1, p2, p3, &Az, &Bz, &Cz, &Dz);

	// only for normal
	p1[2] = (*n1)[0];
	p2[2] = (*n2)[0];
	p3[2] = (*n3)[0];
	calPlaneParams(p1, p2, p3, &A1, &B1, &C1, &D1);
	p1[2] = (*n1)[1];
	p2[2] = (*n2)[1];
	p3[2] = (*n3)[1];
	calPlaneParams(p1, p2, p3, &A2, &B2, &C2, &D2);
	p1[2] = (*n1)[2];
	p2[2] = (*n2)[2];
	p3[2] = (*n3)[2];
	calPlaneParams(p1, p2, p3, &A3, &B3, &C3, &D3);

	// skip if any z value >0 => calculate
	// step 3 computing slope
	slope12.dx_dy = CalSlope((*v1)[X], (*v2)[X], (*v1)[Y], (*v2)[Y], &slope12.fxy);
	slope12.dz_dy = CalSlope((*v1)[Z], (*v2)[Z], (*v1)[Y], (*v2)[Y], &slope12.fzy);
	slope13.dx_dy = CalSlope((*v1)[X], (*v3)[X], (*v1)[Y], (*v3)[Y], &slope13.fxy);
	slope13.dz_dy = CalSlope((*v1)[Z], (*v3)[Z], (*v1)[Y], (*v3)[Y], &slope13.fzy);
	slope23.dx_dy = CalSlope((*v3)[X], (*v2)[X], (*v3)[Y], (*v2)[Y], &slope23.fxy);
	slope23.dz_dy = CalSlope((*v3)[Z], (*v2)[Z], (*v3)[Y], (*v2)[Y], &slope23.fzy);

	// step 4 fill in value
	// part 1
	for (int y = ceil((*v1)[Y]); y < (*v2)[Y]; y++) {

		start_x = CalLoc((*v1)[X], slope12.dx_dy, (*v1)[Y], y, slope12.fxy);
		start_z = CalLoc((*v1)[Z], slope12.dz_dy, (*v1)[Y], y, slope12.fzy);

		end_x = CalLoc((*v1)[X], slope13.dx_dy, (*v1)[Y], y, slope13.fxy);
		end_z = CalLoc((*v1)[Z], slope13.dz_dy, (*v1)[Y], y, slope13.fzy);
		if (end_x < start_x) {
			var = end_x;
			end_x = start_x;
			start_x = var;

			var = end_z;
			end_z = start_z;
			start_z = var;
		}
		slope_dz_dx = CalSlope(end_z, start_z, end_x, start_x, &fzx);

		for (int x = ceil(start_x); x < end_x; x++) {

			z = calZ(x, y, Az, Bz, Cz, Dz);

			//check the state

			GzGet(x, y, &r, &g, &b, &a, &z_pre);

			if (z_pre > z && z >= 0) {
				
				//x,y,z & norm
				N[0] = calZ(x, y, A1, B1, C1, D1);
				N[1] = calZ(x, y, A2, B2, C2, D2);
				N[2] = calZ(x, y, A3, B3, C3, D3);
				top[0] = N[0] + x;
				top[1] = N[1] + y;
				top[2] = N[2] + z;
				z = calZ(x, y - 0.5, Az, Bz, Cz, Dz);
				base1[0] = x;
				base1[1] = y - 0.5;
				base1[2] = z;
				z = calZ(x, y + 0.5, Az, Bz, Cz, Dz);
				base2[0] = x;
				base2[1] = y + 0.5;
				base2[2] = z;
				calPlaneParams(top, base1, base2, &A1, &B1, &C1, &D1);
				normal[0] = A1;
				normal[1] = B1;
				normal[2] = C1;
				t1[0] = 0.5;
				t1[1] = 0;
				t2[0] = -0.5;
				t2[1] = 0;
				t3[0] = 0;
				t3[1] = 1;

			}
		}
	}

	// part2
	for (int y = ceil((*v2)[Y]); y < (*v3)[Y]; y++) {
		start_x = CalLoc((*v2)[X], slope23.dx_dy, (*v2)[Y], y, slope23.fxy);
		start_z = CalLoc((*v2)[Z], slope23.dz_dy, (*v2)[Y], y, slope23.fzy);
		// not right
		end_x = CalLoc((*v1)[X], slope13.dx_dy, (*v1)[Y], y, slope13.fxy);
		end_z = CalLoc((*v1)[Z], slope13.dz_dy, (*v1)[Y], y, slope13.fzy);

		if (end_x < start_x) {
			var = end_x;
			end_x = start_x;
			start_x = var;

			var = end_z;
			end_z = start_z;
			start_z = var;
		}

		slope_dz_dx = CalSlope(end_z, start_z, end_x, start_x, &fzx);

		for (int x = ceil(start_x); x < end_x; x++) {
			z = calZ(x, y, Az, Bz, Cz, Dz);
			// check the state

			GzGet(x, y, &r, &g, &b, &a, &z_pre);

			if (z_pre > z && z >= 0) {
				//x,y,z & norm
				N[0] = calZ(x, y, A1, B1, C1, D1);
				N[1] = calZ(x, y, A2, B2, C2, D2);
				N[2] = calZ(x, y, A3, B3, C3, D3);
				top[0] = N[0] + x;
				top[1] = N[1] + y;
				top[2] = N[2] + z;
				z = calZ(x, y - 0.5, Az, Bz, Cz, Dz);
				base1[0] = x;
				base1[1] = y - 0.5;
				base1[2] = z;
				z = calZ(x, y + 0.5, Az, Bz, Cz, Dz);
				base2[0] = x;
				base2[1] = y + 0.5;
				base2[2] = z;
				calPlaneParams(top, base1, base2, &A1, &B1, &C1, &D1);
				normal[0] = A1;
				normal[1] = B1;
				normal[2] = C1;
				t1[0] = 0.5;
				t1[1] = 0;
				t2[0] = -0.5;
				t2[1] = 0;
				t3[0] = 0;
				t3[1] = 1;

				
			}
		}

	}

	result[0][0] = top[0];
	result[0][1] = top[1];
	result[0][2] = top[2];
	result[0][3] = normal[0];
	result[0][4] = normal[1];
	result[0][5] = normal[2];
	result[0][6] = t1[0];
	result[0][7] = t1[1];

	result[1][0] = base1[0];
	result[1][1] = base1[1];
	result[1][2] = base1[2];
	result[1][3] = normal[0];
	result[1][4] = normal[1];
	result[1][5] = normal[2];
	result[1][6] = t2[0];
	result[1][7] = t2[1];

	result[2][0] = base2[0];
	result[2][1] = base2[1];
	result[2][2] = base2[2];
	result[2][3] = normal[0];
	result[2][4] = normal[1];
	result[2][5] = normal[2];
	result[2][6] = t3[0];
	result[2][7] = t3[1];

	return result;
}

int GzRender::GzRotXMat(float degree, GzMatrix mat)
{
	/* HW 3.1
	// Create rotate matrix : rotate along x axis
	// Pass back the matrix using mat value
	*/
	float radians = degree * PI / 180;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i == j) {
				mat[i][j] = 1;
			}
			else {
				mat[i][j] = 0;
			}
		}
	}
	mat[1][1] = cos(radians);
	mat[1][2] = -sin(radians);
	mat[2][1] = sin(radians);
	mat[2][2] = cos(radians);
	return GZ_SUCCESS;
}

int GzRender::GzRotYMat(float degree, GzMatrix mat)
{
	/* HW 3.2
	// Create rotate matrix : rotate along y axis
	// Pass back the matrix using mat value
	*/
	float radians = degree * PI / 180;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i == j) {
				mat[i][j] = 1;
			}
			else {
				mat[i][j] = 0;
			}
		}
	}
	mat[0][0] = cos(radians);
	mat[0][2] = sin(radians);
	mat[2][0] = -sin(radians);
	mat[2][2] = cos(radians);

	return GZ_SUCCESS;
}

int GzRender::GzRotZMat(float degree, GzMatrix mat)
{
	/* HW 3.3
	// Create rotate matrix : rotate along z axis
	// Pass back the matrix using mat value
	*/
	float radians = degree * PI / 180;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i == j) {
				mat[i][j] = 1;
			}
			else {
				mat[i][j] = 0;
			}
		}
	}
	mat[0][0] = cos(radians);
	mat[0][1] = -sin(radians);
	mat[1][0] = sin(radians);
	mat[1][1] = cos(radians);

	return GZ_SUCCESS;
}

int GzRender::GzTrxMat(GzCoord translate, GzMatrix mat)
{
	/* HW 3.4
	// Create translation matrix
	// Pass back the matrix using mat value
	*/
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i == j) {
				mat[i][j] = 1;
			}
			else {
				mat[i][j] = 0;
			}
		}
	}

	for (int i = 0; i < 3; i++) {
		mat[i][3] = translate[i];
	}

	return GZ_SUCCESS;
}


int GzRender::GzScaleMat(GzCoord scale, GzMatrix mat)
{
	/* HW 3.5
	// Create scaling matrix
	// Pass back the matrix using mat value
	*/
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			mat[i][j] = 0;
		}
	}
	mat[0][0] = scale[0];
	mat[1][1] = scale[1];
	mat[2][2] = scale[2];
	mat[3][3] = 1;

	return GZ_SUCCESS;
}




int GzRender::GzBeginRender()
{
	/* HW 3.7
	- setup for start of each frame - init frame buffer color,alpha,z
	- compute Xiw and projection xform Xpi from camera definition
	- init Ximage - put Xsp at base of stack, push on Xpi and Xiw
	- now stack contains Xsw and app can push model Xforms when needed
	*/
	GzMatrix*        Xiw = &m_camera.Xiw;  		  /* xform from world to image space */
	GzMatrix*        Xpi = &m_camera.Xpi;
	GzCoord x, y, z, CL, up_prime;
	float temp;

	// set matlevel =0
	matlevel = 0;
	// cal Xpi
	float d_reciprocal = tan(m_camera.FOV*PI / 180 / 2);
	//init
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			(*Xpi)[i][j] = 0;
			(*Xiw)[i][j] = 0;
		}
	}
	for (int i = 0; i < 4; i++) {
		(*Xpi)[i][i] = 1;
	}
	(*Xpi)[2][2] = d_reciprocal;
	(*Xpi)[3][2] = d_reciprocal;
	// cal Xiw
	(*Xiw)[3][3] = 1;
	for (int i = 0; i < 3; i++) {
		CL[i] = m_camera.lookat[i] - m_camera.position[i];
	}
	//vecProduct(m_camera.position,m_camera.lookat,(float*)&CL);
	temp = norm2(CL);
	for (int i = 0; i < 3; i++) {
		z[i] = CL[i] / temp;
	}
	temp = dotproduct(m_camera.worldup, z);
	for (int i = 0; i < 3; i++) {
		up_prime[i] = m_camera.worldup[i] - temp * z[i];
	}
	temp = norm2(up_prime);
	for (int i = 0; i < 3; i++) {
		y[i] = up_prime[i] / temp;
	}

	crossProduct(y, z, &x);

	for (int i = 0; i < 3; i++) {
		(*Xiw)[0][i] = x[i];
		(*Xiw)[1][i] = y[i];
		(*Xiw)[2][i] = z[i];
	}
	(*Xiw)[0][3] = -dotproduct(x, m_camera.position);
	(*Xiw)[1][3] = -dotproduct(y, m_camera.position);
	(*Xiw)[2][3] = -dotproduct(z, m_camera.position);


	// init to 0
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			Xsp[i][j] = 0;
		}
	}
	Xsp[0][0] = xres / 2;
	Xsp[0][3] = xres / 2;
	Xsp[1][3] = yres / 2;
	Xsp[1][1] = -yres / 2;
	Xsp[2][2] = (double)INT_MAX;
	Xsp[3][3] = 1;


	GzPushMatrix(Xsp);

	GzPushMatrix(*Xpi);
	GzPushMatrix(*Xiw);

	return GZ_SUCCESS;
}

int GzRender::GzPutCamera(GzCamera camera)
{
	/* HW 3.8
	/*- overwrite renderer camera structure with new camera definition
	*/
	GzCoord*        pos = &m_camera.position;
	GzCoord*        lookat = &m_camera.lookat;
	GzCoord*        worldup = &m_camera.worldup;
	GzMatrix*		Xpi = &m_camera.Xpi;
	GzMatrix*		Xiw = &m_camera.Xiw;
	for (int i = 0; i < 3; i++) {
		(*pos)[i] = camera.position[i];
		(*lookat)[i] = camera.lookat[i];
		(*worldup)[i] = camera.worldup[i];
	}
	m_camera.FOV = camera.FOV;
	//set Xiw Xpi
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			(*Xpi)[i][j] = camera.Xpi[i][j];
			(*Xiw)[i][j] = camera.Xiw[i][j];
		}
	}



	return GZ_SUCCESS;
}
void mulMatrix(GzMatrix x, GzMatrix y, GzMatrix* result) {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			(*result)[i][j] = 0;
			for (int k = 0; k < 4; k++) {
				(*result)[i][j] += x[i][k] * y[k][j];
			}
		}
	}
}
int GzRender::GzPushMatrix(GzMatrix	matrix)
{
	/* HW 3.9
	- push a matrix onto the Ximage stack
	- check for stack overflow
	*/
	GzMatrix* var = &(Ximage[matlevel]);
	GzMatrix* nvar = &(Xnorm[matlevel]);
	float norm;
	// for Ximage
	if (matlevel) {
		mulMatrix(Ximage[matlevel - 1], matrix, var);
	}
	else {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				(*var)[i][j] = matrix[i][j];
			}
		}
	}
	//for Xnorm
	if (matlevel < 2) {
		// Xsp and Xpi identity matrix
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				if (i == j) {
					(*nvar)[i][j] = 1;
				}
				else {
					(*nvar)[i][j] = 0;
				}
			}
		}
	}
	else {

		for (int i = 0; i < 3; i++) {
			matrix[i][3] = 0;
		}
		// norm the norm array
		for (int i = 0; i < 3; i++) {
			// cal norm
			norm = 0;
			for (int j = 0; j < 3; j++) {
				norm += pow(matrix[i][j], 2);
			}
			// norm
			norm = sqrt(norm);
			for (int j = 0; j < 3; j++) {
				matrix[i][j] = matrix[i][j] / norm;
			}

		}
		mulMatrix(Xnorm[matlevel - 1], matrix, nvar);

	}
	matlevel++;

	return GZ_SUCCESS;
}

int GzRender::GzPopMatrix()
{
	/* HW 3.10
	- pop a matrix off the Ximage stack
	- check for stack underflow
	*/
	if (matlevel > 0) { matlevel--; }

	return GZ_SUCCESS;
}






