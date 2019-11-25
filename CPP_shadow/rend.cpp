/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"

#define PI (float) 3.14159265358979323846
#include <iostream>
using namespace std;

void swap(float&, float&);
void sort_base_y(GzCoord, GzCoord, GzCoord);
int getcolor(int, float, GzLight*, GzLight, GzColor, GzColor, GzColor, GzCoord, GzColor, int, GzColor);

float CalSlope(float x1, float x2, float y1, float y2, bool* flag);
float CalLoc(float start_y, float slope, float start_x, float end_x, bool flag);
void calPlaneParams(GzCoord v1, GzCoord v2, GzCoord v3, double* A, double* B, double* C, double* D, double* flag);
void vecProduct(GzCoord x, GzCoord y, GzCoord result);
void crossProduct(GzCoord y, GzCoord z, GzCoord* result);
float dotproduct(GzCoord x, GzCoord y);
float norm2(GzCoord x);
int signcheck(float var);
float calZ(float x, float y, double A, double B, double C, double D, double flag);

void setNormal(GzCoord V1, GzCoord V2, GzCoord V3, GzCoord N1, GzCoord N2, GzCoord N3, double** normalTablePram);
void setZ(GzCoord V1, GzCoord V2, GzCoord V3, GzCoord N1, GzCoord N2, GzCoord N3, double* zTablePram);
float getNormalX(float x, float y, float z, double** normalTablePram);
float getNormalY(float x, float y, float z, double** normalTablePram);
float getNormalZ(float x, float y, float z, double** normalTablePram);
float getZ(float x, float y, float z, double* zTablePram);
float getMidX(GzCoord V1, GzCoord V2, GzCoord V3);
float getMidY(GzCoord V1, GzCoord V2, GzCoord V3);
float getMidZ(GzCoord V1, GzCoord V2, GzCoord V3);

double *zTablePram;
double **normalTablePram;
float **result;
float *z_map;

float rand_angle() {
	return  static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / PI));
}
float rand_height() {
	return  static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 0.1)) + 0.4;
}
float rand_base() {
	return  static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 0.02)) + 0.01;
}
float rand_sign() {
	float a = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 1));
	if (a > 0.5) {
		return 1;
	}
	else {
		return -1;
	}
}
float rand_curve() {
	return  static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 0.1)) + 0.1;
}

int GzRender::GzRotXMat(float degree, GzMatrix mat)
{
	if (degree == NULL || mat == NULL) return GZ_FAILURE;
	double radian = (degree / 180) * PI;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			mat[i][j] = 0;
		}
	}
	mat[0][0] = 1;

	mat[1][1] = cosf(radian);
	mat[1][2] = -sinf(radian);

	mat[2][1] = sinf(radian);
	mat[2][2] = cosf(radian);

	mat[3][3] = 1;
	return GZ_SUCCESS;
}

int GzRender::GzRotYMat(float degree, GzMatrix mat)
{
	if (degree == NULL || mat == NULL) return GZ_FAILURE;
	double radian = (degree / 180) * PI;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			mat[i][j] = 0;
		}
	}
	mat[0][0] = cosf(radian);
	mat[0][2] = sinf(radian);

	mat[1][1] = 1;

	mat[2][0] = -sinf(radian);
	mat[2][2] = cosf(radian);

	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzRender::GzRotZMat(float degree, GzMatrix mat)
{
	if (degree == NULL || mat == NULL) return GZ_FAILURE;
	double radian = (degree / 180) * PI;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			mat[i][j] = 0;
		}
	}
	mat[0][0] = cosf(radian);
	mat[0][1] = -sinf(radian);

	mat[1][0] = sinf(radian);
	mat[1][1] = cosf(radian);

	mat[2][2] = 1;

	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzRender::GzTrxMat(GzCoord translate, GzMatrix mat)
{
	if (translate == NULL || mat == NULL) return GZ_FAILURE;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			mat[i][j] = 0;
		}
	}
	mat[0][0] = 1;
	mat[1][1] = 1;
	mat[2][2] = 1;

	mat[0][3] = translate[0];
	mat[1][3] = translate[1];
	mat[2][3] = translate[2];
	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzRender::GzScaleMat(GzCoord scale, GzMatrix mat)
{
	if (scale == NULL || mat == NULL) return GZ_FAILURE;
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


GzRender::GzRender(int xRes, int yRes)
{
	xres = xRes;
	yres = yRes;

	int resolution = 0;
	resolution = xres * yres;

	framebuffer = (char*)malloc(sizeof(GzPixel) * xres * yres);
	pixelbuffer = (GzPixel*)malloc(sizeof(GzPixel) * xres * yres);

	matlevel = -1;
	lightmatlevel = -1;
	numlights = 0;

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			Xsp[i][j] = 0;
			m_camera.Xiw[i][j] = 0;
			m_camera.Xpi[i][j] = 0;
		}
	}

	Xsp[0][0] = (float)(xres / 2.0);
	Xsp[0][3] = (float)(xres / 2.0);
	Xsp[1][1] = -((float)yres / 2.0);
	Xsp[1][3] = (float)(yres / 2.0);
	Xsp[2][2] = (float)MAXINT;
	Xsp[3][3] = 1.0;

	m_camera.FOV = DEFAULT_FOV;

	m_camera.lookat[0] = 0;
	m_camera.lookat[1] = 0;
	m_camera.lookat[2] = 0;

	m_camera.position[0] = DEFAULT_IM_X;
	m_camera.position[1] = DEFAULT_IM_Y;
	m_camera.position[2] = DEFAULT_IM_Z;

	m_camera.worldup[0] = 0;
	m_camera.worldup[1] = 1;
	m_camera.worldup[2] = 0;

	//shade

	//interp_mode = GZ_RGB_COLOR;

	numlights = 0;

	GzColor ka = DEFAULT_AMBIENT;
	GzColor kd = DEFAULT_DIFFUSE;
	GzColor ks = DEFAULT_SPECULAR;

	Ka[RED] = ka[RED];
	Ka[GREEN] = ka[GREEN];
	Ka[BLUE] = ka[BLUE];

	Kd[RED] = kd[RED];
	Kd[GREEN] = kd[GREEN];
	Kd[BLUE] = kd[BLUE];

	Ks[RED] = ks[RED];
	Ks[GREEN] = ks[GREEN];
	Ks[BLUE] = ks[BLUE];

	spec = DEFAULT_SPEC;


	//shade
	z_map = (float *)malloc(sizeof(float) * xres * yres);

	//grass
	zTablePram = (double *)malloc(5 * sizeof(double));
	normalTablePram = (double **)malloc(3 * sizeof(double*));
	for (int i = 0; i < 3; i++)
		normalTablePram[i] = (double *)malloc(5 * sizeof(double));

	result = (float **)malloc(3 * 5 * sizeof(float *));
	for (int i = 0; i < 3 * 5; i++)
		result[i] = (float *)malloc(8 * sizeof(float));

	// introduce Xlw
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			Xlw[i][j] = 0;
	}
	GzCoord cl;
	for (int i = 0; i < 3; i++)
		cl[i] = m_camera.lookat[i] - m_camera.position[i];

	GzCoord newZ;
	float sq = sqrt(cl[0] * cl[0] + cl[1] * cl[1] + cl[2] * cl[2]);
	for (int i = 0; i < 3; i++)
		newZ[i] = cl[i] / sq;

	//dz = up * z
	float dZ = m_camera.worldup[0] * newZ[0] + m_camera.worldup[1] * newZ[1] + m_camera.worldup[2] * newZ[2];

	GzCoord newup;//up' = up - dz*z

	for (int i = 0; i < 3; i++)
		newup[i] = m_camera.worldup[i] - dZ * newZ[i];

	GzCoord newY;
	float sq2 = sqrt(newup[0] * newup[0] + newup[1] * newup[1] + newup[2] * newup[2]);
	for (int i = 0; i < 3; i++)
		newY[i] = newup[i] / sq2;


	GzCoord newX;
	newX[0] = newY[1] * newZ[2] - newY[2] * newZ[1];
	newX[1] = newY[2] * newZ[0] - newY[0] * newZ[2];
	newX[2] = newY[0] * newZ[1] - newY[1] * newZ[0];

	for (int j = 0; j < 3; j++) {
		for (int i = 0; i < 4; i++) {
			if (j == 0) {
				if (i == 3) {
					m_camera.Xiw[j][i] = -1.0f * (newX[0] * m_camera.position[0] + newX[1] * m_camera.position[1] + newX[2] * m_camera.position[2]);
				}
				else
					m_camera.Xiw[j][i] = newX[i];
			}
			else if (j == 1) {
				if (i == 3) {
					m_camera.Xiw[j][i] = -1.0f * (newY[0] * m_camera.position[0] + newY[1] * m_camera.position[1] + newY[2] * m_camera.position[2]);
				}
				else
					m_camera.Xiw[j][i] = newY[i];
			}
			else if (j == 2) {
				if (i == 3) {
					m_camera.Xiw[j][i] = -1.0f * (newZ[0] * m_camera.position[0] + newZ[1] * m_camera.position[1] + newZ[2] * m_camera.position[2]);
				}
				else
					m_camera.Xiw[j][i] = newZ[i];
			}
		}
	}

	m_camera.Xiw[3][3] = 1.0;

}

GzRender::~GzRender()
{
	free(framebuffer);
	free(pixelbuffer);
	free(z_map);
	free(zTablePram);
	free(normalTablePram);
	free(result);
}

int GzRender::GzDefault()
{
	/* HW1.3 set pixel buffer to some default values - start a new frame */

	for (int i = 0; i < xres * yres; i++)
	{
		pixelbuffer[i] = { 3200, 3500, 3700, 1, MAXINT };
		z_map[i] = MAXINT;
		framebuffer[3 * i] = (char)3200;
		framebuffer[3 * i + 1] = (char)3500;
		framebuffer[3 * i + 2] = (char)3700;
	}

	return GZ_SUCCESS;
}

int GzRender::GzBeginRender()
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			m_camera.Xpi[i][j] = 0;
	}

	float di = tanf((m_camera.FOV * PI / 180.0) / 2.0);
	m_camera.Xpi[0][0] = 1.0;
	m_camera.Xpi[1][1] = 1.0;
	m_camera.Xpi[2][2] = di;
	m_camera.Xpi[3][2] = di;
	m_camera.Xpi[3][3] = 1.0;
	//iw
	GzCoord cl;
	for (int i = 0; i < 3; i++)
		cl[i] = m_camera.lookat[i] - m_camera.position[i];

	GzCoord newZ;
	float sq = sqrt(cl[0] * cl[0] + cl[1] * cl[1] + cl[2] * cl[2]);
	for (int i = 0; i < 3; i++)
		newZ[i] = cl[i] / sq;

	//dz = up * z
	float dZ = m_camera.worldup[0] * newZ[0] + m_camera.worldup[1] * newZ[1] + m_camera.worldup[2] * newZ[2];

	GzCoord newup;//up' = up - dz*z

	for (int i = 0; i < 3; i++)
		newup[i] = m_camera.worldup[i] - dZ * newZ[i];

	GzCoord newY;
	float sq2 = sqrt(newup[0] * newup[0] + newup[1] * newup[1] + newup[2] * newup[2]);
	for (int i = 0; i < 3; i++)
		newY[i] = newup[i] / sq2;


	GzCoord newX;
	newX[0] = newY[1] * newZ[2] - newY[2] * newZ[1];
	newX[1] = newY[2] * newZ[0] - newY[0] * newZ[2];
	newX[2] = newY[0] * newZ[1] - newY[1] * newZ[0];

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++) {
			m_camera.Xiw[i][j] = 0;
			m_camera.Xwi[i][j] = 0;
		}
	}

	for (int j = 0; j < 3; j++) {
		for (int i = 0; i < 4; i++) {
			if (j == 0) {
				if (i == 3) {
					m_camera.Xiw[j][i] = -1.0f * (newX[0] * m_camera.position[0] + newX[1] * m_camera.position[1] + newX[2] * m_camera.position[2]);
				}
				else {
					m_camera.Xiw[j][i] = newX[i];
					m_camera.Xwi[i][j] = newX[i];
				}

			}
			else if (j == 1) {
				if (i == 3) {
					m_camera.Xiw[j][i] = -1.0f * (newY[0] * m_camera.position[0] + newY[1] * m_camera.position[1] + newY[2] * m_camera.position[2]);
				}
				else {
					m_camera.Xiw[j][i] = newY[i];
					m_camera.Xwi[i][j] = newY[i];
				}

			}
			else if (j == 2) {
				if (i == 3) {
					m_camera.Xiw[j][i] = -1.0f * (newZ[0] * m_camera.position[0] + newZ[1] * m_camera.position[1] + newZ[2] * m_camera.position[2]);
				}
				else {
					m_camera.Xiw[j][i] = newZ[i];
					m_camera.Xwi[i][j] = newZ[i];
				}
			}
		}
	}

	m_camera.Xiw[3][3] = 1.0;
	m_camera.Xwi[3][3] = 1.0;

	GzPushMatrix(Xsp);
	GzPushMatrix(m_camera.Xpi);
	GzPushMatrix(m_camera.Xiw);

	//change light to world space
	for (int ii = 0; ii < 3; ii++) {
		m_light.lookat[ii] = 0;
		for (int jj = 0; jj < 3; jj++) {
			m_light.lookat[ii] += m_camera.Xwi[ii][jj] * lights[0].direction[jj];
		}
	}
	m_light.lookat[0] = 7.8;
	m_light.lookat[1] = 0.7;
	m_light.lookat[2] = 6.5;


	// add light stack
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			m_light.Xpi[i][j] = 0;
	}

	di = tanf((m_light.FOV * PI / 180.0) / 2.0);
	m_light.Xpi[0][0] = 1.0;
	m_light.Xpi[1][1] = 1.0;
	m_light.Xpi[2][2] = di;
	m_light.Xpi[3][2] = di;
	m_light.Xpi[3][3] = 1.0;
	//iw

	for (int i = 0; i < 3; i++)
		cl[i] = m_light.lookat[i] - m_light.position[i];


	sq = sqrt(cl[0] * cl[0] + cl[1] * cl[1] + cl[2] * cl[2]);
	for (int i = 0; i < 3; i++)
		newZ[i] = cl[i] / sq;

	//dz = up * z
	dZ = m_light.worldup[0] * newZ[0] + m_light.worldup[1] * newZ[1] + m_light.worldup[2] * newZ[2];



	for (int i = 0; i < 3; i++)
		newup[i] = m_light.worldup[i] - dZ * newZ[i];


	sq2 = sqrt(newup[0] * newup[0] + newup[1] * newup[1] + newup[2] * newup[2]);
	for (int i = 0; i < 3; i++)
		newY[i] = newup[i] / sq2;



	newX[0] = newY[1] * newZ[2] - newY[2] * newZ[1];
	newX[1] = newY[2] * newZ[0] - newY[0] * newZ[2];
	newX[2] = newY[0] * newZ[1] - newY[1] * newZ[0];

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++) {
			m_light.Xiw[i][j] = 0;
		}
	}

	for (int j = 0; j < 3; j++) {
		for (int i = 0; i < 4; i++) {
			if (j == 0) {
				if (i == 3) {
					m_light.Xiw[j][i] = -1.0f * (newX[0] * m_light.position[0] + newX[1] * m_light.position[1] + newX[2] * m_light.position[2]);
				}
				else {
					m_light.Xiw[j][i] = newX[i];
				}

			}
			else if (j == 1) {
				if (i == 3) {
					m_light.Xiw[j][i] = -1.0f * (newY[0] * m_light.position[0] + newY[1] * m_light.position[1] + newY[2] * m_light.position[2]);
				}
				else {
					m_light.Xiw[j][i] = newY[i];
				}

			}
			else if (j == 2) {
				if (i == 3) {
					m_light.Xiw[j][i] = -1.0f * (newZ[0] * m_light.position[0] + newZ[1] * m_light.position[1] + newZ[2] * m_light.position[2]);
				}
				else {
					m_light.Xiw[j][i] = newZ[i];
				}
			}
		}
	}

	m_light.Xiw[3][3] = 1.0;
	GzPushLightMatrix(Xsp);
	GzPushLightMatrix(m_light.Xpi);
	GzPushLightMatrix(m_light.Xiw);

	return GZ_SUCCESS;
}

int GzRender::GzPutCamera(GzCamera camera)
{
	/* HW 3.8
	/*- overwrite renderer camera structure with new camera definition
	*/

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			m_camera.Xiw[i][j] = camera.Xiw[i][j];
			m_camera.Xpi[i][j] = camera.Xpi[i][j];
		}
	}

	m_camera.position[X] = camera.position[X];
	m_camera.position[Y] = camera.position[Y];
	m_camera.position[Z] = camera.position[Z];
	m_camera.lookat[X] = camera.lookat[X];
	m_camera.lookat[Y] = camera.lookat[Y];
	m_camera.lookat[Z] = camera.lookat[Z];
	m_camera.worldup[X] = camera.worldup[X];
	m_camera.worldup[Y] = camera.worldup[Y];
	m_camera.worldup[Z] = camera.worldup[Z];
	m_camera.FOV = camera.FOV;

	return GZ_SUCCESS;
}

int GzRender::GzPutLight(GzCamera light)
{
	/* HW 3.8
	/*- overwrite renderer camera structure with new camera definition
	*/

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			light.Xiw[i][j] = light.Xiw[i][j];
			light.Xpi[i][j] = light.Xpi[i][j];
		}
	}

	m_light.position[X] = light.position[X];
	m_light.position[Y] = light.position[Y];
	m_light.position[Z] = light.position[Z];
	m_light.worldup[X] = light.worldup[X];
	m_light.worldup[Y] = light.worldup[Y];
	m_light.worldup[Z] = light.worldup[Z];
	m_light.FOV = light.FOV;

	return GZ_SUCCESS;
}
int GzRender::GzPushLightMatrix(GzMatrix matrix)
{
	if (lightmatlevel < MATLEVELS)
	{
		//identity_mat
		if (lightmatlevel != -1)
		{
			//multiply
			float temp = 0;
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					temp = 0;
					for (int k = 0; k < 4; k++) {
						temp = temp + XLight[lightmatlevel][i][k] * matrix[k][j];
					}
					XLight[lightmatlevel + 1][i][j] = temp;
				}
			}
		}
		else
		{
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					XLight[0][i][j] = matrix[i][j];
				}
			}
		}
	}
	lightmatlevel++;
	return GZ_SUCCESS;
}
int GzRender::GzPushMatrix(GzMatrix	matrix)
{
	/* HW 3.9
	- push a matrix onto the Ximage stack
	- check for stack overflow
	*/

	if (matlevel >= MATLEVELS || matrix == NULL) return GZ_FAILURE;
	if (matlevel < MATLEVELS)
	{
		//identity_mat
		GzMatrix identity_mat;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				if (i == j) {
					identity_mat[i][j] = 1;
				}
				else
					identity_mat[i][j] = 0;
			}
		}


		if (matlevel != -1)
		{
			//multiply
			float temp = 0;
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					temp = 0;
					for (int k = 0; k < 4; k++) {
						temp = temp + Ximage[matlevel][i][k] * matrix[k][j];
					}
					Ximage[matlevel + 1][i][j] = temp;
				}
			}
		}
		else
		{
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					Ximage[0][i][j] = matrix[i][j];
				}
			}
		}

		if (matlevel == -1 || matlevel == 0) {

			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					Xnorm[matlevel + 1][i][j] = (float)identity_mat[i][j];
				}
			}
		}
		else {
			GzMatrix uni_res;
			GzMatrix uni_res_xnorm;
			uni_res[0][0] = matrix[0][0];
			uni_res[0][1] = matrix[0][1];
			uni_res[0][2] = matrix[0][2];
			uni_res[0][3] = 0;
			uni_res[1][0] = matrix[1][0];
			uni_res[1][1] = matrix[1][1];
			uni_res[1][2] = matrix[1][2];
			uni_res[1][3] = 0;
			uni_res[2][0] = matrix[2][0];
			uni_res[2][1] = matrix[2][1];
			uni_res[2][2] = matrix[2][2];
			uni_res[2][3] = 0;
			uni_res[3][0] = matrix[3][0];
			uni_res[3][1] = matrix[3][1];
			uni_res[3][2] = matrix[3][2];
			uni_res[3][3] = 1;

			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					uni_res_xnorm[i][j] = 0.0;
					for (int k = 0; k < 4; k++) {
						uni_res_xnorm[i][j] += Xnorm[matlevel][i][k] * uni_res[k][j];
					}
				}

			}

			float temp_k = 1.0 / sqrt((uni_res_xnorm[0][0] * uni_res_xnorm[0][0]) + (uni_res_xnorm[0][1] * uni_res_xnorm[0][1]) + (uni_res_xnorm[0][2] * uni_res_xnorm[0][2]));
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					uni_res_xnorm[i][j] *= temp_k;
				}
			}

			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					Xnorm[matlevel + 1][i][j] = (float)uni_res_xnorm[i][j];
				}
			}
		}

		matlevel++;

	}
	else
		return GZ_FAILURE;

	return GZ_SUCCESS;
}

int GzRender::GzPopMatrix()
{
	/* HW 3.10
	- pop a matrix off the Ximage stack
	- check for stack underflow
	*/

	if (matlevel < 0)
		return GZ_FAILURE;
	matlevel--;
	return GZ_SUCCESS;
}

int GzRender::GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
	/* HW1.4 write pixel values into the buffer */

	int indexI, indexJ;
	GzIntensity red, green, blue, alpha;
	GzDepth depth;

	indexI = max(min(i, xres - 1), 0);
	indexJ = max(min(j, yres - 1), 0);

	red = max(min(r, 4095), 0);
	green = max(min(g, 4095), 0);
	blue = max(min(b, 4095), 0);

	alpha = a;
	depth = z;

	if (pixelbuffer == NULL) return GZ_FAILURE;

	if (z < pixelbuffer[ARRAY(indexI, indexJ)].z) {
		pixelbuffer[ARRAY(indexI, indexJ)].red = red;
		pixelbuffer[ARRAY(indexI, indexJ)].green = green;
		pixelbuffer[ARRAY(indexI, indexJ)].blue = blue;
		pixelbuffer[ARRAY(indexI, indexJ)].alpha = alpha;
		pixelbuffer[ARRAY(indexI, indexJ)].z = depth;
	}

	return GZ_SUCCESS;
}


int GzRender::GzGet(int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z)
{
	/* HW1.5 retrieve a pixel information from the pixel buffer */

	int indexI = max(min(i, xres - 1), 0);
	int indexJ = max(min(j, yres - 1), 0);

	if (pixelbuffer == NULL) return GZ_FAILURE;

	*r = pixelbuffer[ARRAY(indexI, indexJ)].red;
	*g = pixelbuffer[ARRAY(indexI, indexJ)].green;
	*b = pixelbuffer[ARRAY(indexI, indexJ)].blue;
	*a = pixelbuffer[ARRAY(indexI, indexJ)].alpha;
	*z = pixelbuffer[ARRAY(indexI, indexJ)].z;

	return GZ_SUCCESS;
}


int GzRender::GzFlushDisplay2File(FILE* outfile)
{
	/* HW1.6 write image to ppm file -- "P6 %d %d 255\r" */

	fprintf(outfile, "P6 %d %d 255\n", xres, yres);

	if (pixelbuffer == NULL) return GZ_FAILURE;

	for (int i = 0; i <= (xres * yres); i++) {

		GzPixel p = pixelbuffer[i];
		char pred, pgreen, pblue;
		size_t size = 1, count = 1;

		pred = p.red >> 4;
		pgreen = p.green >> 4;
		pblue = p.blue >> 4;

		fwrite(&pred, size, count, outfile);
		fwrite(&pgreen, size, count, outfile);
		fwrite(&pblue, size, count, outfile);
	}

	return GZ_SUCCESS;
}
int GzRender::GzFlushDisplay2DepthFile(FILE* outfile)
{
	/* HW1.6 write image to ppm file -- "P6 %d %d 255\r" */

	fprintf(outfile, "P6 %d %d 255\n", xres, yres);

	CString msg;
	msg.Format(_T("z_map, %f\n"), z_map[0]);
	AfxMessageBox(msg);

	if (z_map == NULL) return GZ_FAILURE;
	float min_z = MAXINT;
	float span_z = MAXINT;
	GzDepth p;
	for (int i = 0; i <= (xres * yres); i++) {
		p = z_map[i];
		if (p < min_z)
			min_z = (float)p;
	}

	msg.Format(_T("min_z, %f\n"), min_z);
	AfxMessageBox(msg);

	span_z -= min_z;

	for (int i = 0; i <= (xres * yres); i++) {
		p = z_map[i];
		char pred, pgreen, pblue;
		size_t size = 1, count = 1;

		pred = (p - min_z) / span_z * 255.0;;
		pgreen = (p - min_z) / span_z * 255.0;;
		pblue = (p - min_z) / span_z * 255.0;

		fwrite(&pred, size, count, outfile);
		fwrite(&pgreen, size, count, outfile);
		fwrite(&pblue, size, count, outfile);
	}

	return GZ_SUCCESS;
}

int GzRender::GzFlushDisplay2FrameBuffer()
{
	/* HW1.7 write pixels to framebuffer:
		- put the pixels into the frame buffer
		- CAUTION: when storing the pixels into the frame buffer, the order is blue, green, and red
		- NOT red, green, and blue !!!
	*/

	int index = 0;

	if (pixelbuffer == NULL || framebuffer == NULL) return GZ_FAILURE;

	for (int i = 0; i <= (xres * yres); i++) {

		GzPixel p = pixelbuffer[i];
		char fred, fgreen, fblue;

		fred = p.red >> 4;
		fgreen = p.green >> 4;
		fblue = p.blue >> 4;

		framebuffer[index++] = fblue;
		framebuffer[index++] = fgreen;
		framebuffer[index++] = fred;
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

	/*
	- GzPutAttribute() must accept the following tokens/values:

	- GZ_RGB_COLOR					GzColor		default flat-shade color
	- GZ_INTERPOLATE				int			shader interpolation mode
	- GZ_DIRECTIONAL_LIGHT			GzLight
	- GZ_AMBIENT_LIGHT            	GzLight		(ignore direction)
	- GZ_AMBIENT_COEFFICIENT		GzColor		Ka reflectance
	- GZ_DIFFUSE_COEFFICIENT		GzColor		Kd reflectance
	- GZ_SPECULAR_COEFFICIENT       GzColor		Ks coef's
	- GZ_DISTRIBUTION_COEFFICIENT   float		spec power
	*/

	if (pixelbuffer == NULL || framebuffer == NULL || numAttributes == NULL || *nameList == NULL || *valueList == NULL) return GZ_FAILURE;

	int i = 0;
	//GzIntensity red, green, blue;

	while (i < numAttributes)
	{
		GzToken token = nameList[i];
		if (token == GZ_RGB_COLOR)
		{
			GzColor* c = (GzColor*)valueList[i];
			//float *c = (float*)*(valueList+i);
			flatcolor[0] = max(min(c[i][RED], 4095), 0);
			flatcolor[1] = max(min(c[i][GREEN], 4095), 0);
			flatcolor[2] = max(min(c[i][BLUE], 4095), 0);
		}
		else if (token == GZ_INTERPOLATE)
		{
			int *m = (int*)valueList[i];
			interp_mode = *m;
		}
		else if (token == GZ_DIRECTIONAL_LIGHT)
		{
			GzLight *dir = (GzLight*)valueList[i];
			lights[numlights].direction[0] = dir->direction[0];
			lights[numlights].direction[1] = dir->direction[1];
			lights[numlights].direction[2] = dir->direction[2];
			lights[numlights].color[0] = dir->color[0]; //le
			lights[numlights].color[1] = dir->color[1];
			lights[numlights].color[2] = dir->color[2];
			numlights++;

		}
		else if (token == GZ_AMBIENT_LIGHT)
		{
			GzLight *amb = (GzLight*)valueList[i];
			ambientlight.direction[0] = amb->direction[0];
			ambientlight.direction[1] = amb->direction[1];
			ambientlight.direction[2] = amb->direction[2];
			ambientlight.color[0] = amb->color[0]; //la
			ambientlight.color[1] = amb->color[1];
			ambientlight.color[2] = amb->color[2];
		}
		else if (token == GZ_AMBIENT_COEFFICIENT)
		{
			float *ac = (float*)valueList[i];
			Ka[0] = ac[0];
			Ka[1] = ac[1];
			Ka[2] = ac[2];
		}
		else if (token == GZ_DIFFUSE_COEFFICIENT)
		{
			float *dc = (float*)valueList[i];
			Kd[0] = dc[0];
			Kd[1] = dc[1];
			Kd[2] = dc[2];
		}
		else if (token == GZ_SPECULAR_COEFFICIENT)
		{
			float *sc = (float*)valueList[i];
			Ks[0] = sc[0];
			Ks[1] = sc[1];
			Ks[2] = sc[2];
		}
		else if (token == GZ_DISTRIBUTION_COEFFICIENT)
		{
			float *distrc = (float*)valueList[i];
			spec = *distrc;
		}
		else if (token == GZ_TEXTURE_MAP)
		{
			GzTexture uv = (GzTexture)valueList[i];
			tex_fun = uv;
		}
		i++;
	}

	return GZ_SUCCESS;
}

int GzRender::GzPutTriangle(int numParts, GzToken *nameList, GzPointer *valueList)
{
	for (int i = 0; i < numParts; i++)
	{
		if (pixelbuffer == NULL || framebuffer == NULL) return GZ_FAILURE;
		GzCoord* vP = (GzCoord*)valueList[0];	//position
		GzCoord* nP = (GzCoord*)valueList[1];	//norm
		GzTextureIndex* uvP = (GzTextureIndex*)valueList[2];

		float v4d[3][4];
		float n4d[3][4];
		//GzCoord v4d[4];
		v4d[0][0] = vP[0][X];
		v4d[0][1] = vP[0][Y];
		v4d[0][2] = vP[0][Z];
		v4d[1][0] = vP[1][X];
		v4d[1][1] = vP[1][Y];
		v4d[1][2] = vP[1][Z];
		v4d[2][0] = vP[2][X];
		v4d[2][1] = vP[2][Y];
		v4d[2][2] = vP[2][Z];
		v4d[0][3] = 1.0f;
		v4d[1][3] = 1.0f;
		v4d[2][3] = 1.0f;

		n4d[0][0] = nP[0][X];
		n4d[0][1] = nP[0][Y];
		n4d[0][2] = nP[0][Z];
		n4d[1][0] = nP[1][X];
		n4d[1][1] = nP[1][Y];
		n4d[1][2] = nP[1][Z];
		n4d[2][0] = nP[2][X];
		n4d[2][1] = nP[2][Y];
		n4d[2][2] = nP[2][Z];
		n4d[0][3] = 1.0f;
		n4d[1][3] = 1.0f;
		n4d[2][3] = 1.0f;

		float t4d[3][4];
		float nt4d[3][4];

		float s = 0;
		float a = 0;

		//matrix multiply
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 4; k++) {
					s = s + Ximage[matlevel][i][k] * v4d[j][k];
					a = a + Xnorm[matlevel][i][k] * n4d[j][k];
				}
				t4d[j][i] = s;
				nt4d[j][i] = a;
				s = 0;
				a = 0;
			}
		}

		GzCoord v1 = { t4d[0][0] / t4d[0][3], t4d[0][1] / t4d[0][3], t4d[0][2] / t4d[0][3] }; //x, y, z
		GzCoord v2 = { t4d[1][0] / t4d[1][3], t4d[1][1] / t4d[1][3], t4d[1][2] / t4d[1][3] };
		GzCoord v3 = { t4d[2][0] / t4d[2][3], t4d[2][1] / t4d[2][3], t4d[2][2] / t4d[2][3] };

		GzCoord n1 = { nt4d[0][0] / nt4d[0][3], nt4d[0][1] / nt4d[0][3], nt4d[0][2] / nt4d[0][3] }; //x, y, z
		GzCoord n2 = { nt4d[1][0] / nt4d[1][3], nt4d[1][1] / nt4d[1][3], nt4d[1][2] / nt4d[1][3] };
		GzCoord n3 = { nt4d[2][0] / nt4d[2][3], nt4d[2][1] / nt4d[2][3], nt4d[2][2] / nt4d[2][3] };

		GzCoord n1_bak = { n1[0] / 1 , n1[1] / 1 , n1[2] / 1 };
		GzCoord n2_bak = { n2[0] / 1 , n2[1] / 1 , n2[2] / 1 };
		GzCoord n3_bak = { n3[0] / 1 , n3[1] / 1 , n3[2] / 1 };

		GzColor color_v1 = { 0, 0, 0 };
		GzColor color_v2 = { 0, 0, 0 };
		GzColor color_v3 = { 0, 0, 0 };

		GzTextureIndex uv1 = { uvP[0][0], uvP[0][1] };
		GzTextureIndex uv2 = { uvP[1][0], uvP[1][1] };
		GzTextureIndex uv3 = { uvP[2][0], uvP[2][1] };

		if (tex_fun != NULL && interp_mode == GZ_COLOR)
		{
			for (int k = 0; k < 3; k++) {
				Ka[k] = 1.0;
				Ks[k] = 1.0;
				Kd[k] = 1.0;
			}
		}

		getcolor(numlights, spec, lights, ambientlight, Ka, Kd, Ks, n1, color_v1, 0, NULL);
		getcolor(numlights, spec, lights, ambientlight, Ka, Kd, Ks, n2, color_v2, 0, NULL);
		getcolor(numlights, spec, lights, ambientlight, Ka, Kd, Ks, n3, color_v3, 0, NULL);

		//sort_base_y(v1, v2, v3);
		if (v1[1] > v2[1]) {
			for (int i = 0; i < 3; i++) {
				swap(v1[i], v2[i]); //xyz
				swap(n1[i], n2[i]); //xyz vector
				if (interp_mode != GZ_FLAT)
					swap(color_v1[i], color_v2[i]); //rgb
			}
			swap(uv1[0], uv2[0]);
			swap(uv1[1], uv2[1]);
		}

		if (v2[1] > v3[1]) {
			for (int i = 0; i < 3; i++) {
				swap(v2[i], v3[i]);
				swap(n2[i], n3[i]);
				if (interp_mode != GZ_FLAT)
					swap(color_v2[i], color_v3[i]);
			}
			swap(uv3[0], uv2[0]);
			swap(uv3[1], uv2[1]);
		}

		if (v1[1] > v2[1]) {
			for (int i = 0; i < 3; i++) {
				swap(v1[i], v2[i]);
				swap(n1[i], n2[i]);
				if (interp_mode != GZ_FLAT)
					swap(color_v1[i], color_v2[i]);
			}
			swap(uv1[0], uv2[0]);
			swap(uv1[1], uv2[1]);
		}
		//sort_base_y end

		int flag = -1;

		//divide special case
		if (v1[Y] == v2[Y]) //flag invert tri
		{
			if (v2[X] > v1[X]) {
				for (int i = 0; i < 3; i++) {
					swap(v1[i], v2[i]);
					swap(n1[i], n2[i]);
					if (interp_mode != GZ_FLAT)
						swap(color_v1[i], color_v2[i]);
				}
				swap(uv1[0], uv2[0]);
				swap(uv1[1], uv2[1]);
			}

			flag = 0; //invert tri
		}
		else if (v2[Y] == v3[Y]) //flag invert tri
		{
			if (v2[X] > v3[X]) {
				//v1 always on the upper left
				for (int i = 0; i < 3; i++) {
					swap(v3[i], v2[i]);
					swap(n3[i], n2[i]);
					if (interp_mode != GZ_FLAT)
						swap(color_v2[i], color_v3[i]);
				}
				swap(uv3[0], uv2[0]);
				swap(uv3[1], uv2[1]);
			}

			flag = 1; //tri 
		}
		else {

			flag = 2;
		}

		float uvVz;
		GzTextureIndex Puv[3];

		Puv[0][0] = uv1[0] / (v1[2] / (MAXINT - v1[2]) + 1.0);
		Puv[0][1] = uv1[1] / (v1[2] / (MAXINT - v1[2]) + 1.0);

		Puv[1][0] = uv2[0] / (v2[2] / (MAXINT - v2[2]) + 1.0);
		Puv[1][1] = uv2[1] / (v2[2] / (MAXINT - v2[2]) + 1.0);

		Puv[2][0] = uv3[0] / (v3[2] / (MAXINT - v3[2]) + 1.0);
		Puv[2][1] = uv3[1] / (v3[2] / (MAXINT - v3[2]) + 1.0);

		// Get the plane to interpolate u
		float uX_1 = v2[0] - v1[0], uY_1 = v2[1] - v1[1], uZ_1 = Puv[1][0] - Puv[0][0];
		float uX_2 = v3[0] - v1[0], uY_2 = v3[1] - v1[1], uZ_2 = Puv[2][0] - Puv[0][0];
		float uPanelA = (uY_1 * uZ_2 - uZ_1 * uY_2);
		float uPanelB = (uZ_1 * uX_2 - uX_1 * uZ_2);
		float uPanelC = (uX_1 * uY_2 - uY_1 * uX_2);
		float uPanelD = -1.0f * (uPanelA * v1[0] + uPanelB * v1[1] + uPanelC * Puv[0][0]);
		// Get the plane to interpolate v
		float vX_1 = v2[0] - v1[0], vY_1 = v2[1] - v1[1], vZ_1 = Puv[1][1] - Puv[0][1];
		float vX_2 = v3[0] - v1[0], vY_2 = v3[1] - v1[1], vZ_2 = Puv[2][1] - Puv[0][1];
		float vPanelA = (vY_1 * vZ_2 - vZ_1 * vY_2);
		float vPanelB = (vZ_1 * vX_2 - vX_1 * vZ_2);
		float vPanelC = (vX_1 * vY_2 - vY_1 * vX_2);
		float vPanelD = -1.0f * (vPanelA * v1[0] + vPanelB * v1[1] + vPanelC * Puv[0][1]);


		//render based on case
		if (flag == 2) {
			float s12x = (v2[X] - v1[X]) / (v2[Y] - v1[Y]);
			float s12z = (v2[Z] - v1[Z]) / (v2[Y] - v1[Y]);

			float s13x = (v3[X] - v1[X]) / (v3[Y] - v1[Y]);
			float s13z = (v3[Z] - v1[Z]) / (v3[Y] - v1[Y]);

			float s23x = (v3[X] - v2[X]) / (v3[Y] - v2[Y]);
			float s23z = (v3[Z] - v2[Z]) / (v3[Y] - v2[Y]);

			float deltaY = 0;

			float s12c_r, s13c_r, s23c_r, s12c_g, s13c_g, s23c_g, s12c_b, s13c_b, s23c_b;
			float s12n_x, s13n_x, s23n_x, s12n_y, s13n_y, s23n_y, s12n_z, s13n_z, s23n_z;

			if (interp_mode == GZ_COLOR) {
				s12c_r = (color_v2[RED] - color_v1[RED]) / (v2[Y] - v1[Y]);
				s13c_r = (color_v3[RED] - color_v1[RED]) / (v3[Y] - v1[Y]);
				s23c_r = (color_v3[RED] - color_v2[RED]) / (v3[Y] - v2[Y]);

				s12c_g = (color_v2[GREEN] - color_v1[GREEN]) / (v2[Y] - v1[Y]);
				s13c_g = (color_v3[GREEN] - color_v1[GREEN]) / (v3[Y] - v1[Y]);
				s23c_g = (color_v3[GREEN] - color_v2[GREEN]) / (v3[Y] - v2[Y]);

				s12c_b = (color_v2[BLUE] - color_v1[BLUE]) / (v2[Y] - v1[Y]);
				s13c_b = (color_v3[BLUE] - color_v1[BLUE]) / (v3[Y] - v1[Y]);
				s23c_b = (color_v3[BLUE] - color_v2[BLUE]) / (v3[Y] - v2[Y]);
			}
			else if (interp_mode == GZ_NORMALS) {
				s12n_x = (n2[X] - n1[X]) / (v2[Y] - v1[Y]);
				s13n_x = (n3[X] - n1[X]) / (v3[Y] - v1[Y]);
				s23n_x = (n3[X] - n2[X]) / (v3[Y] - v2[Y]);

				s12n_y = (n2[Y] - n1[Y]) / (v2[Y] - v1[Y]);
				s13n_y = (n3[Y] - n1[Y]) / (v3[Y] - v1[Y]);
				s23n_y = (n3[Y] - n2[Y]) / (v3[Y] - v2[Y]);

				s12n_z = (n2[Z] - n1[Z]) / (v2[Y] - v1[Y]);
				s13n_z = (n3[Z] - n1[Z]) / (v3[Y] - v1[Y]);
				s23n_z = (n3[Z] - n2[Z]) / (v3[Y] - v2[Y]);
			}


			for (int j = ceil(v1[Y]); j <= v3[Y]; j++) { //y

				deltaY = j - v1[Y];

				if (deltaY <= v2[Y] - v1[Y]) {
					float x1 = v1[X] + s12x * deltaY;
					float x2 = v1[X] + s13x * deltaY;
					float zValue1 = v1[Z] + s12z * deltaY;
					float zValue2 = v1[Z] + s13z * deltaY;

					float c_rValue1, c_rValue2, c_gValue1, c_gValue2, c_bValue1, c_bValue2;
					float n_xValue1, n_xValue2, n_yValue1, n_yValue2, n_zValue1, n_zValue2;

					if (interp_mode == GZ_COLOR) {
						c_rValue1 = color_v1[RED] + s12c_r * deltaY;
						c_rValue2 = color_v1[RED] + s13c_r * deltaY;

						c_gValue1 = color_v1[GREEN] + s12c_g * deltaY;
						c_gValue2 = color_v1[GREEN] + s13c_g * deltaY;

						c_bValue1 = color_v1[BLUE] + s12c_b * deltaY;
						c_bValue2 = color_v1[BLUE] + s13c_b * deltaY;
					}
					else if (interp_mode == GZ_NORMALS) {
						n_xValue1 = n1[X] + s12n_x * deltaY;
						n_xValue2 = n1[X] + s13n_x * deltaY;

						n_yValue1 = n1[Y] + s12n_y * deltaY;
						n_yValue2 = n1[Y] + s13n_y * deltaY;

						n_zValue1 = n1[Z] + s12n_z * deltaY;
						n_zValue2 = n1[Z] + s13n_z * deltaY;
					}

					if (x1 > x2) {
						swap(x1, x2);
						swap(zValue1, zValue2);
						if (interp_mode == GZ_COLOR) {
							swap(c_rValue1, c_rValue2);
							swap(c_gValue1, c_gValue2);
							swap(c_bValue1, c_bValue2);
						}
						else if (interp_mode == GZ_NORMALS) {
							swap(n_xValue1, n_xValue2);
							swap(n_yValue1, n_yValue2);
							swap(n_zValue1, n_zValue2);
						}
					}

					//cal z slope in one y
					float slopeZX = (zValue2 - zValue1) / (x2 - x1);
					float slopec_rX, slopec_gX, slopec_bX;
					float slopen_xX, slopen_yX, slopen_zX;
					if (interp_mode == GZ_COLOR) {
						slopec_rX = (c_rValue2 - c_rValue1) / (x2 - x1);
						slopec_gX = (c_gValue2 - c_gValue1) / (x2 - x1);
						slopec_bX = (c_bValue2 - c_bValue1) / (x2 - x1);
					}
					else if (interp_mode == GZ_NORMALS) {
						slopen_xX = (n_xValue2 - n_xValue1) / (x2 - x1);
						slopen_yX = (n_yValue2 - n_yValue1) / (x2 - x1);
						slopen_zX = (n_zValue2 - n_zValue1) / (x2 - x1);
					}

					for (int i = ceil(x1); i <= x2; i++) { //x		

						GzIntensity r = 0, g = 0, b = 0, a = 0;
						GzDepth z = MAXINT;



						if (!(i >= xres || i < 0 || j >= yres || j < 0)) {
							GzGet(i, j, &r, &g, &b, &a, &z);
							float newZ = zValue1 + slopeZX * (i - x1);

							GzTextureIndex cUV;
							float texture[3];

							if (tex_fun != NULL) {
								float para_uv = (float)newZ / ((float)MAXINT - (float)newZ);
								cUV[0] = -(uPanelA*i + uPanelB * j + uPanelD) / uPanelC;
								cUV[1] = -(vPanelA*i + vPanelB * j + vPanelD) / vPanelC;

								cUV[0] *= (para_uv + 1.0);
								cUV[1] *= (para_uv + 1.0);
								tex_fun(cUV[0], cUV[1], texture);
							}



							if (newZ < z && newZ >= 0) {
								//if (newZ < z) {
								if (interp_mode == GZ_FLAT) {
									r = max(min((GzIntensity)ctoi(color_v1[0]), 4095), 0);
									g = max(min((GzIntensity)ctoi(color_v1[1]), 4095), 0);
									b = max(min((GzIntensity)ctoi(color_v1[2]), 4095), 0);
								}
								else if (interp_mode == GZ_COLOR) {
									float newR = c_rValue1 + slopec_rX * (i - x1);
									float newG = c_gValue1 + slopec_gX * (i - x1);
									float newB = c_bValue1 + slopec_bX * (i - x1);

									if (tex_fun != NULL) {
										newR *= texture[0];
										newG *= texture[1];
										newB *= texture[2];
									}

									r = max(min((GzIntensity)ctoi(newR), 4095), 0);
									g = max(min((GzIntensity)ctoi(newG), 4095), 0);
									b = max(min((GzIntensity)ctoi(newB), 4095), 0);
								}
								else if (interp_mode == GZ_NORMALS) {
									float newX = n_xValue1 + slopen_xX * (i - x1);
									float newY = n_yValue1 + slopen_yX * (i - x1);
									float newZ = n_zValue1 + slopen_zX * (i - x1);
									GzCoord n_temp = { newX, newY, newZ };
									GzColor color_temp = { 0, 0, 0 };

									getcolor(numlights, spec, lights, ambientlight, Ka, Kd, Ks, n_temp, color_temp, (tex_fun != NULL) ? 1 : 0, texture);

									r = max(min((GzIntensity)ctoi(color_temp[RED]), 4095), 0);
									g = max(min((GzIntensity)ctoi(color_temp[GREEN]), 4095), 0);
									b = max(min((GzIntensity)ctoi(color_temp[BLUE]), 4095), 0);
								}

								GzPut(i, j, r, g, b, 1, newZ);
							}

						}
					}
				}
				else {
					float x1 = v2[X] + s23x * (deltaY - (v2[Y] - v1[Y]));
					float x2 = v1[X] + s13x * deltaY;
					float zValue1 = v2[Z] + s23z * (deltaY - (v2[Y] - v1[Y]));
					float zValue2 = v1[Z] + s13z * deltaY;

					float c_rValue1, c_rValue2, c_gValue1, c_gValue2, c_bValue1, c_bValue2;
					float n_xValue1, n_xValue2, n_yValue1, n_yValue2, n_zValue1, n_zValue2;

					if (interp_mode == GZ_COLOR) {
						c_rValue1 = color_v2[RED] + s23c_r * (deltaY - (v2[Y] - v1[Y]));
						c_rValue2 = color_v1[RED] + s13c_r * deltaY;

						c_gValue1 = color_v2[GREEN] + s23c_g * (deltaY - (v2[Y] - v1[Y]));
						c_gValue2 = color_v1[GREEN] + s13c_g * deltaY;

						c_bValue1 = color_v2[BLUE] + s23c_b * (deltaY - (v2[Y] - v1[Y]));
						c_bValue2 = color_v1[BLUE] + s13c_b * deltaY;
					}
					else if (interp_mode == GZ_NORMALS) {
						n_xValue1 = n2[X] + s23n_x * (deltaY - (v2[Y] - v1[Y]));
						n_xValue2 = n1[X] + s13n_x * deltaY;

						n_yValue1 = n2[Y] + s23n_y * (deltaY - (v2[Y] - v1[Y]));
						n_yValue2 = n1[Y] + s13n_y * deltaY;

						n_zValue1 = n2[Z] + s23n_z * (deltaY - (v2[Y] - v1[Y]));
						n_zValue2 = n1[Z] + s13n_z * deltaY;
					}

					if (x1 > x2) {
						swap(x1, x2);
						swap(zValue1, zValue2);
						if (interp_mode == GZ_COLOR) {
							swap(c_rValue1, c_rValue2);
							swap(c_gValue1, c_gValue2);
							swap(c_bValue1, c_bValue2);
						}
						else if (interp_mode == GZ_NORMALS) {
							swap(n_xValue1, n_xValue2);
							swap(n_yValue1, n_yValue2);
							swap(n_zValue1, n_zValue2);
						}
					}

					//cal z slope in one y
					float slopeZX = (zValue2 - zValue1) / (x2 - x1);
					float slopec_rX, slopec_gX, slopec_bX;
					float slopen_xX, slopen_yX, slopen_zX;

					if (interp_mode == GZ_COLOR) {
						slopec_rX = (c_rValue2 - c_rValue1) / (x2 - x1);
						slopec_gX = (c_gValue2 - c_gValue1) / (x2 - x1);
						slopec_bX = (c_bValue2 - c_bValue1) / (x2 - x1);
					}
					else if (interp_mode == GZ_NORMALS) {
						slopen_xX = (n_xValue2 - n_xValue1) / (x2 - x1);
						slopen_yX = (n_yValue2 - n_yValue1) / (x2 - x1);
						slopen_zX = (n_zValue2 - n_zValue1) / (x2 - x1);
					}

					for (int i = ceil(x1); i <= x2; i++) { //x		

						GzIntensity r = 0, g = 0, b = 0, a = 0;
						GzDepth z = MAXINT;

						if (!(i >= xres || i < 0 || j >= yres || j < 0)) {
							GzGet(i, j, &r, &g, &b, &a, &z);
							float newZ = zValue1 + slopeZX * (i - x1);

							GzTextureIndex cUV;
							float texture[3];

							if (tex_fun != NULL) {
								float para_uv = (float)newZ / ((float)MAXINT - (float)newZ);
								cUV[0] = -(uPanelA*i + uPanelB * j + uPanelD) / uPanelC;
								cUV[1] = -(vPanelA*i + vPanelB * j + vPanelD) / vPanelC;

								cUV[0] *= (para_uv + 1.0);
								cUV[1] *= (para_uv + 1.0);
								tex_fun(cUV[0], cUV[1], texture);
							}


							if (newZ < z && newZ >= 0) {
								//if (newZ < z) {
								if (interp_mode == GZ_FLAT) {
									r = max(min((GzIntensity)ctoi(color_v1[0]), 4095), 0);
									g = max(min((GzIntensity)ctoi(color_v1[1]), 4095), 0);
									b = max(min((GzIntensity)ctoi(color_v1[2]), 4095), 0);
								}
								else if (interp_mode == GZ_COLOR) {
									float newR = c_rValue1 + slopec_rX * (i - x1);
									float newG = c_gValue1 + slopec_gX * (i - x1);
									float newB = c_bValue1 + slopec_bX * (i - x1);

									if (tex_fun != NULL) {
										newR *= texture[0];
										newG *= texture[1];
										newB *= texture[2];
									}

									r = max(min((GzIntensity)ctoi(newR), 4095), 0);
									g = max(min((GzIntensity)ctoi(newG), 4095), 0);
									b = max(min((GzIntensity)ctoi(newB), 4095), 0);
								}
								else if (interp_mode == GZ_NORMALS) {
									float newX = n_xValue1 + slopen_xX * (i - x1);
									float newY = n_yValue1 + slopen_yX * (i - x1);
									float newZ = n_zValue1 + slopen_zX * (i - x1);
									GzCoord n_temp = { newX, newY, newZ };
									GzColor color_temp = { 0, 0, 0 };

									getcolor(numlights, spec, lights, ambientlight, Ka, Kd, Ks, n_temp, color_temp, (tex_fun != NULL) ? 1 : 0, texture);

									r = max(min((GzIntensity)ctoi(color_temp[RED]), 4095), 0);
									g = max(min((GzIntensity)ctoi(color_temp[GREEN]), 4095), 0);
									b = max(min((GzIntensity)ctoi(color_temp[BLUE]), 4095), 0);
								}
								GzPut(i, j, r, g, b, 1, newZ);
							}

						}
					}
				}
				//
			}
		}

		if (flag == 1) {	//line V1; line v2, v3
			float s12x = (v2[X] - v1[X]) / (v2[Y] - v1[Y]);
			float s12z = (v2[Z] - v1[Z]) / (v2[Y] - v1[Y]);
			float s13x = (v3[X] - v1[X]) / (v3[Y] - v1[Y]);
			float s13z = (v3[Z] - v1[Z]) / (v3[Y] - v1[Y]);
			float deltaY = 0;

			float s12c_r, s13c_r, s12c_g, s13c_g, s12c_b, s13c_b;
			float s12n_x, s13n_x, s12n_y, s13n_y, s12n_z, s13n_z;


			if (interp_mode == GZ_COLOR) {
				s12c_r = (color_v2[RED] - color_v1[RED]) / (v2[Y] - v1[Y]);
				s13c_r = (color_v3[RED] - color_v1[RED]) / (v3[Y] - v1[Y]);

				s12c_g = (color_v2[GREEN] - color_v1[GREEN]) / (v2[Y] - v1[Y]);
				s13c_g = (color_v3[GREEN] - color_v1[GREEN]) / (v3[Y] - v1[Y]);

				s12c_b = (color_v2[BLUE] - color_v1[BLUE]) / (v2[Y] - v1[Y]);
				s13c_b = (color_v3[BLUE] - color_v1[BLUE]) / (v3[Y] - v1[Y]);
			}
			else if (interp_mode == GZ_NORMALS) {
				s12n_x = (n2[X] - n1[X]) / (v2[Y] - v1[Y]);
				s13n_x = (n3[X] - n1[X]) / (v3[Y] - v1[Y]);

				s12n_y = (n2[Y] - n1[Y]) / (v2[Y] - v1[Y]);
				s13n_y = (n3[Y] - n1[Y]) / (v3[Y] - v1[Y]);

				s12n_z = (n2[Z] - n1[Z]) / (v2[Y] - v1[Y]);
				s13n_z = (n3[Z] - n1[Z]) / (v3[Y] - v1[Y]);
			}

			for (int j = ceil(v1[Y]); j <= v3[Y]; j++) { //y

				deltaY = j - v1[Y];

				float x1 = v1[X] + s12x * deltaY;
				float x2 = v1[X] + s13x * deltaY;
				float zValue1 = v1[Z] + s12z * deltaY;
				float zValue2 = v1[Z] + s13z * deltaY;

				float c_rValue1, c_rValue2, c_gValue1, c_gValue2, c_bValue1, c_bValue2;
				float n_xValue1, n_xValue2, n_yValue1, n_yValue2, n_zValue1, n_zValue2;

				if (interp_mode == GZ_COLOR) {
					c_rValue1 = color_v1[RED] + s12c_r * deltaY;
					c_rValue2 = color_v1[RED] + s13c_r * deltaY;

					c_gValue1 = color_v1[GREEN] + s12c_g * deltaY;
					c_gValue2 = color_v1[GREEN] + s13c_g * deltaY;

					c_bValue1 = color_v1[BLUE] + s12c_b * deltaY;
					c_bValue2 = color_v1[BLUE] + s13c_b * deltaY;
				}
				else if (interp_mode == GZ_NORMALS) {
					n_xValue1 = n1[X] + s12n_x * deltaY;
					n_xValue2 = n1[X] + s13n_x * deltaY;

					n_yValue1 = n1[Y] + s12n_y * deltaY;
					n_yValue2 = n1[Y] + s13n_y * deltaY;

					n_zValue1 = n1[Z] + s12n_z * deltaY;
					n_zValue2 = n1[Z] + s13n_z * deltaY;
				}

				if (x1 > x2) {
					swap(x1, x2);
					swap(zValue1, zValue2);

					if (interp_mode == GZ_COLOR) {
						swap(c_rValue1, c_rValue2);
						swap(c_gValue1, c_gValue2);
						swap(c_bValue1, c_bValue2);
					}
					else if (interp_mode == GZ_NORMALS) {
						swap(n_xValue1, n_xValue2);
						swap(n_yValue1, n_yValue2);
						swap(n_zValue1, n_zValue2);
					}
				}

				//cal z slope in one y
				float slopeZX = (zValue2 - zValue1) / (x2 - x1);
				float slopec_rX, slopec_gX, slopec_bX;
				float slopen_xX, slopen_yX, slopen_zX;

				if (interp_mode == GZ_COLOR) {
					slopec_rX = (c_rValue2 - c_rValue1) / (x2 - x1);
					slopec_gX = (c_gValue2 - c_gValue1) / (x2 - x1);
					slopec_bX = (c_bValue2 - c_bValue1) / (x2 - x1);
				}
				else if (interp_mode == GZ_NORMALS) {
					slopen_xX = (n_xValue2 - n_xValue1) / (x2 - x1);
					slopen_yX = (n_yValue2 - n_yValue1) / (x2 - x1);
					slopen_zX = (n_zValue2 - n_zValue1) / (x2 - x1);
				}

				for (int i = ceil(x1); i <= x2; i++) { //x		

					GzIntensity r = 0, g = 0, b = 0, a = 0;
					GzDepth z = MAXINT;

					if (!(i >= xres || i < 0 || j >= yres || j < 0)) {
						GzGet(i, j, &r, &g, &b, &a, &z);
						float newZ = zValue1 + slopeZX * (i - x1);

						GzTextureIndex cUV;
						float texture[3];

						if (tex_fun != NULL) {
							float para_uv = (float)newZ / ((float)MAXINT - (float)newZ);
							cUV[0] = -(uPanelA*i + uPanelB * j + uPanelD) / uPanelC;
							cUV[1] = -(vPanelA*i + vPanelB * j + vPanelD) / vPanelC;

							cUV[0] *= (para_uv + 1.0);
							cUV[1] *= (para_uv + 1.0);
							tex_fun(cUV[0], cUV[1], texture);
						}


						if (newZ < z && newZ >= 0) {
							//if (newZ < z) {
							if (interp_mode == GZ_FLAT) {
								r = max(min((GzIntensity)ctoi(color_v1[0]), 4095), 0);
								g = max(min((GzIntensity)ctoi(color_v1[1]), 4095), 0);
								b = max(min((GzIntensity)ctoi(color_v1[2]), 4095), 0);
							}
							else if (interp_mode == GZ_COLOR) {
								float newR = c_rValue1 + slopec_rX * (i - x1);
								float newG = c_gValue1 + slopec_gX * (i - x1);
								float newB = c_bValue1 + slopec_bX * (i - x1);

								if (tex_fun != NULL) {
									newR *= texture[0];
									newG *= texture[1];
									newB *= texture[2];
								}

								r = max(min((GzIntensity)ctoi(newR), 4095), 0);
								g = max(min((GzIntensity)ctoi(newG), 4095), 0);
								b = max(min((GzIntensity)ctoi(newB), 4095), 0);
							}
							else if (interp_mode == GZ_NORMALS) {
								float newX = n_xValue1 + slopen_xX * (i - x1);
								float newY = n_yValue1 + slopen_yX * (i - x1);
								float newZ = n_zValue1 + slopen_zX * (i - x1);
								GzCoord n_temp = { newX, newY, newZ };
								GzColor color_temp = { 0, 0, 0 };

								getcolor(numlights, spec, lights, ambientlight, Ka, Kd, Ks, n_temp, color_temp, (tex_fun != NULL) ? 1 : 0, texture);

								r = max(min((GzIntensity)ctoi(color_temp[RED]), 4095), 0);
								g = max(min((GzIntensity)ctoi(color_temp[GREEN]), 4095), 0);
								b = max(min((GzIntensity)ctoi(color_temp[BLUE]), 4095), 0);
							}
							GzPut(i, j, r, g, b, 1, newZ);
						}

					}
				}
			}
		}

		if (flag == 0) { //line: V2, V1; line V3
			float s13x = (v3[X] - v1[X]) / (v3[Y] - v1[Y]);
			float s13z = (v3[Z] - v1[Z]) / (v3[Y] - v1[Y]);
			float s23x = (v3[X] - v2[X]) / (v3[Y] - v2[Y]);
			float s23z = (v3[Z] - v2[Z]) / (v3[Y] - v2[Y]);
			float deltaY = 0;

			float s13c_r, s23c_r, s13c_g, s23c_g, s13c_b, s23c_b;
			float s13n_x, s23n_x, s13n_y, s23n_y, s13n_z, s23n_z;


			if (interp_mode == GZ_COLOR) {
				s13c_r = (color_v3[RED] - color_v1[RED]) / (v3[Y] - v1[Y]);
				s23c_r = (color_v3[RED] - color_v2[RED]) / (v3[Y] - v2[Y]);

				s13c_g = (color_v3[GREEN] - color_v1[GREEN]) / (v3[Y] - v1[Y]);
				s23c_g = (color_v3[GREEN] - color_v2[GREEN]) / (v3[Y] - v2[Y]);

				s13c_b = (color_v3[BLUE] - color_v1[BLUE]) / (v3[Y] - v1[Y]);
				s23c_b = (color_v3[BLUE] - color_v2[BLUE]) / (v3[Y] - v2[Y]);
			}
			else if (interp_mode == GZ_NORMALS) {
				s13n_x = (n3[X] - n1[X]) / (v3[Y] - v1[Y]);
				s23n_x = (n3[X] - n2[X]) / (v3[Y] - v2[Y]);

				s13n_y = (n3[Y] - n1[Y]) / (v3[Y] - v1[Y]);
				s23n_y = (n3[Y] - n2[Y]) / (v3[Y] - v2[Y]);

				s13n_z = (n3[Z] - n1[Z]) / (v3[Y] - v1[Y]);
				s23n_z = (n3[Z] - n2[Z]) / (v3[Y] - v2[Y]);
			}

			for (int j = ceil(v1[Y]); j <= v3[Y]; j++) { //y

				deltaY = j - v1[Y];

				float x1 = v2[X] + s23x * deltaY;
				float x2 = v1[X] + s13x * deltaY;
				float zValue1 = v2[Z] + s23z * deltaY;
				float zValue2 = v1[Z] + s13z * deltaY;

				float c_rValue1, c_rValue2, c_gValue1, c_gValue2, c_bValue1, c_bValue2;
				float n_xValue1, n_xValue2, n_yValue1, n_yValue2, n_zValue1, n_zValue2;

				if (interp_mode == GZ_COLOR) {
					c_rValue1 = color_v2[RED] + s23c_r * deltaY;
					c_rValue2 = color_v1[RED] + s13c_r * deltaY;

					c_gValue1 = color_v2[GREEN] + s23c_g * deltaY;
					c_gValue2 = color_v1[GREEN] + s13c_g * deltaY;

					c_bValue1 = color_v2[BLUE] + s23c_b * deltaY;
					c_bValue2 = color_v1[BLUE] + s13c_b * deltaY;
				}
				else if (interp_mode == GZ_NORMALS) {
					n_xValue1 = n2[X] + s23n_x * deltaY;
					n_xValue2 = n1[X] + s13n_x * deltaY;

					n_yValue1 = n2[Y] + s23n_y * deltaY;
					n_yValue2 = n1[Y] + s13n_y * deltaY;

					n_zValue1 = n2[Z] + s23n_z * deltaY;
					n_zValue2 = n1[Z] + s13n_z * deltaY;
				}

				if (x1 > x2) {
					swap(x1, x2);
					swap(zValue1, zValue2);
					if (interp_mode == GZ_COLOR) {
						swap(c_rValue1, c_rValue2);
						swap(c_gValue1, c_gValue2);
						swap(c_bValue1, c_bValue2);
					}
					else if (interp_mode == GZ_NORMALS) {
						swap(n_xValue1, n_xValue2);
						swap(n_yValue1, n_yValue2);
						swap(n_zValue1, n_zValue2);
					}
				}

				//cal z slope in one y
				float slopeZX = (zValue2 - zValue1) / (x2 - x1);
				float slopec_rX, slopec_gX, slopec_bX;
				float slopen_xX, slopen_yX, slopen_zX;

				if (interp_mode == GZ_COLOR) {
					slopec_rX = (c_rValue2 - c_rValue1) / (x2 - x1);
					slopec_gX = (c_gValue2 - c_gValue1) / (x2 - x1);
					slopec_bX = (c_bValue2 - c_bValue1) / (x2 - x1);
				}
				else if (interp_mode == GZ_NORMALS) {
					slopen_xX = (n_xValue2 - n_xValue1) / (x2 - x1);
					slopen_yX = (n_yValue2 - n_yValue1) / (x2 - x1);
					slopen_zX = (n_zValue2 - n_zValue1) / (x2 - x1);
				}

				for (int i = ceil(x1); i <= x2; i++) { //x		

					GzIntensity r = 0, g = 0, b = 0, a = 0;
					GzDepth z = MAXINT;

					if (!(i >= xres || i < 0 || j >= yres || j < 0)) {
						GzGet(i, j, &r, &g, &b, &a, &z);
						float newZ = zValue1 + slopeZX * (i - x1);

						GzTextureIndex cUV;
						float texture[3];

						if (tex_fun != NULL) {
							float para_uv = (float)newZ / ((float)MAXINT - (float)newZ);
							cUV[0] = -(uPanelA*i + uPanelB * j + uPanelD) / uPanelC;
							cUV[1] = -(vPanelA*i + vPanelB * j + vPanelD) / vPanelC;

							cUV[0] *= (para_uv + 1.0);
							cUV[1] *= (para_uv + 1.0);
							tex_fun(cUV[0], cUV[1], texture);
						}


						if (newZ < z && newZ >= 0) {
							//if (newZ < z) {
							if (interp_mode == GZ_FLAT) {
								r = max(min((GzIntensity)ctoi(color_v1[0]), 4095), 0);
								g = max(min((GzIntensity)ctoi(color_v1[1]), 4095), 0);
								b = max(min((GzIntensity)ctoi(color_v1[2]), 4095), 0);
							}
							else if (interp_mode == GZ_COLOR) {
								float newR = c_rValue1 + slopec_rX * (i - x1);
								float newG = c_gValue1 + slopec_gX * (i - x1);
								float newB = c_bValue1 + slopec_bX * (i - x1);

								if (tex_fun != NULL) {
									newR *= texture[0];
									newG *= texture[1];
									newB *= texture[2];
								}

								r = max(min((GzIntensity)ctoi(newR), 4095), 0);
								g = max(min((GzIntensity)ctoi(newG), 4095), 0);
								b = max(min((GzIntensity)ctoi(newB), 4095), 0);
							}
							else if (interp_mode == GZ_NORMALS) {
								float newX = n_xValue1 + slopen_xX * (i - x1);
								float newY = n_yValue1 + slopen_yX * (i - x1);
								float newZ = n_zValue1 + slopen_zX * (i - x1);
								GzCoord n_temp = { newX, newY, newZ };
								GzColor color_temp = { 0, 0, 0 };
								getcolor(numlights, spec, lights, ambientlight, Ka, Kd, Ks, n_temp, color_temp, (tex_fun != NULL) ? 1 : 0, texture);
								r = max(min((GzIntensity)ctoi(color_temp[RED]), 4095), 0);
								g = max(min((GzIntensity)ctoi(color_temp[GREEN]), 4095), 0);
								b = max(min((GzIntensity)ctoi(color_temp[BLUE]), 4095), 0);
							}
							GzPut(i, j, r, g, b, 1, newZ);
						}

					}
				}
			}

		}
	}


	return GZ_SUCCESS;
}
int GzRender::GzPutZ(int i, int j, float z) {
	/* HW1.4 write pixel values into the buffer */

	int indexI, indexJ;

	GzDepth depth;

	indexI = max(min(i, xres - 1), 0);
	indexJ = max(min(j, yres - 1), 0);

	depth = z;

	if (z_map == NULL) return GZ_FAILURE;

	if (z < z_map[ARRAY(indexI, indexJ)]) {

		z_map[ARRAY(indexI, indexJ)] = depth;
	}

	CString msg;
	msg.Format(_T("GzPutZ, %f\n"), z_map[ARRAY(indexI, indexJ)]);
	//AfxMessageBox(msg);

	return GZ_SUCCESS;
}
int GzRender::GzGetZ(int i, int j, float *z) {

	int indexI = max(min(i, xres - 1), 0);
	int indexJ = max(min(j, yres - 1), 0);

	if (z_map == NULL) return GZ_FAILURE;

	*z = z_map[ARRAY(indexI, indexJ)];

	CString msg;
	msg.Format(_T("GzGetZ, %f\n"), z_map[ARRAY(indexI, indexJ)]);
	//AfxMessageBox(msg);

	return GZ_SUCCESS;

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

int GzRender::GzCalShadowDepth(int numParts, GzToken *nameList, GzPointer *valueList)
/* numParts - how many names and values */
{
	//for (int i = 0; i < numParts; i++)
	//{
	//	if (pixelbuffer == NULL || framebuffer == NULL) return GZ_FAILURE;

	//	GzCoord* vP = (GzCoord*)valueList[0];
	//	float v4d[3][4];

	//	//GzCoord v4d[4];
	//	v4d[0][0] = vP[0][X];
	//	v4d[0][1] = vP[0][Y];
	//	v4d[0][2] = vP[0][Z];
	//	v4d[1][0] = vP[1][X];
	//	v4d[1][1] = vP[1][Y];
	//	v4d[1][2] = vP[1][Z];
	//	v4d[2][0] = vP[2][X];
	//	v4d[2][1] = vP[2][Y];
	//	v4d[2][2] = vP[2][Z];
	//	v4d[0][3] = 1.0;
	//	v4d[1][3] = 1.0;
	//	v4d[2][3] = 1.0;

	//	float t4d[3][4];
	//	float s = 0;

	//	//matrix multiply
	//	for (int i = 0; i < 4; i++)
	//	{
	//		for (int j = 0; j < 3; j++)
	//		{
	//			for (int k = 0; k < 4; k++)
	//				s = s + XLight[lightmatlevel][i][k] * v4d[j][k];

	//			t4d[j][i] = s;
	//			s = 0;
	//		}
	//	}

	//	if (t4d[0][2] < m_light.position[2] || t4d[1][2] < m_light.position[2] || t4d[2][2] < m_light.position[2]) return GZ_FAILURE;

	//	GzCoord v1 = { t4d[0][0] / t4d[0][3], t4d[0][1] / t4d[0][3], t4d[0][2] / t4d[0][3] }; //x, y, z
	//	GzCoord v2 = { t4d[1][0] / t4d[1][3], t4d[1][1] / t4d[1][3], t4d[1][2] / t4d[1][3] };
	//	GzCoord v3 = { t4d[2][0] / t4d[2][3], t4d[2][1] / t4d[2][3], t4d[2][2] / t4d[2][3] };

	//	sort_base_y(v1, v2, v3);

	//	int flag = -1;

	//	if (v1[Y] == v2[Y]) //flag invert tri
	//	{
	//		if (v2[X] > v1[X]) {
	//			for (int i = 0; i < 3; i++) {
	//				swap(v1[i], v2[i]);

	//			}
	//		}

	//		flag = 0; //invert tri
	//	}
	//	else if (v2[Y] == v3[Y]) //flag invert tri
	//	{
	//		if (v2[X] > v3[X]) {
	//			//v1 always on the upper left
	//			for (int i = 0; i < 3; i++) {
	//				swap(v3[i], v2[i]);
	//			}
	//		}

	//		flag = 1; //tri 
	//	}
	//	else {

	//		flag = 2;
	//	}

	//	if (flag == 2) {
	//		float s12x = (v2[X] - v1[X]) / (v2[Y] - v1[Y]);
	//		float s12z = (v2[Z] - v1[Z]) / (v2[Y] - v1[Y]);
	//		float s13x = (v3[X] - v1[X]) / (v3[Y] - v1[Y]);
	//		float s13z = (v3[Z] - v1[Z]) / (v3[Y] - v1[Y]);
	//		float s23x = (v3[X] - v2[X]) / (v3[Y] - v2[Y]);
	//		float s23z = (v3[Z] - v2[Z]) / (v3[Y] - v2[Y]);
	//		float deltaY = 0;


	//		for (int j = ceil(v1[Y]); j <= v3[Y]; j++) { //y

	//			deltaY = j - v1[Y];

	//			if (deltaY <= v2[Y] - v1[Y]) {
	//				float x1 = v1[X] + s12x * deltaY;
	//				float x2 = v1[X] + s13x * deltaY;
	//				float zValue1 = v1[Z] + s12z * deltaY;
	//				float zValue2 = v1[Z] + s13z * deltaY;

	//				if (x1 > x2) {
	//					swap(x1, x2);
	//					swap(zValue1, zValue2);
	//				}

	//				//cal z slope in one y
	//				float slopeZX = (zValue2 - zValue1) / (x2 - x1);

	//				for (int i = ceil(x1); i <= x2; i++) { //x		

	//					GzIntensity r = 0, g = 0, b = 0, a = 0;
	//					GzDepth z = MAXINT;

	//					if (!(i >= xres || i < 0 || j >= yres || j < 0)) {
	//						GzGetZ(i, j, &z);
	//						float newZ = zValue1 + slopeZX * (i - x1);
	//						if (newZ < z && newZ >= 0) {
	//							GzPutZ(i, j, newZ);
	//						}

	//					}
	//				}
	//			}
	//			else {
	//				float x1 = v2[X] + s23x * (deltaY - (v2[Y] - v1[Y]));
	//				float x2 = v1[X] + s13x * deltaY;
	//				float zValue1 = v2[Z] + s23z * (deltaY - (v2[Y] - v1[Y]));
	//				float zValue2 = v1[Z] + s13z * deltaY;

	//				if (x1 > x2) {
	//					swap(x1, x2);
	//					swap(zValue1, zValue2);
	//				}

	//				//cal z slope in one y
	//				float slopeZX = (zValue2 - zValue1) / (x2 - x1);

	//				for (int i = ceil(x1); i <= x2; i++) { //x		

	//					GzIntensity r = 0, g = 0, b = 0, a = 0;
	//					GzDepth z = MAXINT;

	//					if (!(i >= xres || i < 0 || j >= yres || j < 0)) {
	//						GzGetZ(i, j, &z);
	//						float newZ = zValue1 + slopeZX * (i - x1);
	//						if (newZ < z && newZ >= 0) {
	//							GzPutZ(i, j, newZ);
	//						}

	//					}
	//				}
	//			}
	//			//
	//		}
	//	}

	//	if (flag == 1) {	//line V1; line v2, v3
	//		float s12x = (v2[X] - v1[X]) / (v2[Y] - v1[Y]);
	//		float s12z = (v2[Z] - v1[Z]) / (v2[Y] - v1[Y]);
	//		float s13x = (v3[X] - v1[X]) / (v3[Y] - v1[Y]);
	//		float s13z = (v3[Z] - v1[Z]) / (v3[Y] - v1[Y]);
	//		float deltaY = 0;

	//		for (int j = ceil(v1[Y]); j <= v3[Y]; j++) { //y

	//			deltaY = j - v1[Y];

	//			float x1 = v1[X] + s12x * deltaY;
	//			float x2 = v1[X] + s13x * deltaY;
	//			float zValue1 = v1[Z] + s12z * deltaY;
	//			float zValue2 = v1[Z] + s13z * deltaY;

	//			if (x1 > x2) {
	//				swap(x1, x2);
	//				swap(zValue1, zValue2);
	//			}

	//			//cal z slope in one y
	//			float slopeZX = (zValue2 - zValue1) / (x2 - x1);

	//			for (int i = ceil(x1); i <= x2; i++) { //x		

	//				GzIntensity r = 0, g = 0, b = 0, a = 0;
	//				GzDepth z = MAXINT;

	//				if (!(i >= xres || i < 0 || j >= yres || j < 0)) {
	//					GzGetZ(i, j, &z);
	//					float newZ = zValue1 + slopeZX * (i - x1);
	//					if (newZ < z && newZ >= 0) {
	//						GzPutZ(i, j, newZ);
	//					}

	//				}
	//			}
	//		}
	//	}

	//	if (flag == 0) { //line: V2, V1; line V3
	//		float s13x = (v3[X] - v1[X]) / (v3[Y] - v1[Y]);
	//		float s13z = (v3[Z] - v1[Z]) / (v3[Y] - v1[Y]);
	//		float s23x = (v3[X] - v2[X]) / (v3[Y] - v2[Y]);
	//		float s23z = (v3[Z] - v2[Z]) / (v3[Y] - v2[Y]);
	//		float deltaY = 0;

	//		for (int j = ceil(v1[Y]); j <= v3[Y]; j++) { //y

	//			deltaY = j - v1[Y];

	//			float x1 = v2[X] + s23x * deltaY;
	//			float x2 = v1[X] + s13x * deltaY;
	//			float zValue1 = v2[Z] + s23z * deltaY;
	//			float zValue2 = v1[Z] + s13z * deltaY;

	//			if (x1 > x2) {
	//				swap(x1, x2);
	//				swap(zValue1, zValue2);
	//			}

	//			//cal z slope in one y
	//			float slopeZX = (zValue2 - zValue1) / (x2 - x1);

	//			for (int i = ceil(x1); i <= x2; i++) { //x		

	//				GzIntensity r = 0, g = 0, b = 0, a = 0;
	//				GzDepth z = MAXINT;

	//				if (!(i >= xres || i < 0 || j >= yres || j < 0)) {
	//					GzGetZ(i, j, &z);
	//					float newZ = zValue1 + slopeZX * (i - x1);
	//					if (newZ < z && newZ >= 0) {
	//						GzPutZ(i, j, newZ);
	//					}

	//				}
	//			}
	//		}

	//	}
	//}

	typedef struct {
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

	GzDepth z_pre, z;
	float z_n;
	double Az, Bz, Cz, Dz,flag;

	GzCoord p1, p2, p3, N;
	Coord3* loc;



	loc = (Coord3*)*valueList;

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


		// calculate the light 
		for (int i = 0; i < 3; i++) {
			p1[i] = (*v1)[i];
			p2[i] = (*v2)[i];
			p3[i] = (*v3)[i];
		}
		calPlaneParams(p1, p2, p3, &Az, &Bz, &Cz, &Dz,&flag);

		// skip if any z value >0 => calculate
		  ////  // step 3 computing slope
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

				z = calZ(x, y, Az, Bz, Cz, Dz,flag);
				 // check the state
				GzGetZ(x, y, &z_n);

				if (z_n > z && z >= 0) {
					//  change u,v to the projective domain
					GzPutZ(x, y, z);
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
				z = calZ(x, y, Az, Bz, Cz, Dz, flag);
				// check the state
				GzGetZ(x, y, &z_n);

				if (z_n > z && z >= 0) {
					//  change u,v to the projective domain
					GzPutZ(x, y, z);

				}
			}

		}

	}
	return GZ_SUCCESS;
}

void rotation(GzCoord a, GzCoord b, GzCoord N, GzCoord*base1, GzCoord* base2) {
	float M[3][3];
	float alpha = rand_angle();
	for (int i = 0; i < 3; i++) {
		M[i][i] = pow(N[i], 2) + pow(1 - N[i], 2)*cos(alpha);
		(*base1)[i] = 0;
	}

	M[0][1] = N[0] * N[1] * (1 - cos(alpha)) - N[2] * sin(alpha);
	M[1][0] = N[0] * N[1] * (1 - cos(alpha)) + N[2] * sin(alpha);

	M[0][2] = N[0] * N[2] * (1 - cos(alpha)) + N[1] * sin(alpha);
	M[2][0] = N[0] * N[2] * (1 - cos(alpha)) - N[1] * sin(alpha);

	M[1][2] = N[1] * N[2] * (1 - cos(alpha)) - N[0] * sin(alpha);
	M[2][1] = N[1] * N[2] * (1 - cos(alpha)) + N[0] * sin(alpha);

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			(*base1)[i] += M[i][j] * a[j];
		}
		(*base2)[i] = -(*base1)[i];
		(*base1)[i] += b[i];
		(*base2)[i] += b[i];
	}

}

float** GzRender::GzAddGrassWithModelSpace(int numParts, GzToken *nameList, GzPointer *valueList) {


	typedef GzCoord   Coord3[3];
	int sorted_y[3] = { 0,1,2 };

	GzCoord *v1, *v2, *v3;
	GzCoord *n1, *n2, *n3;
	float offsetx, offsety;
	double A, B, C, D, flag;
	int temp, first_vertix;
	Coord3* loc;

	//vertex
	loc = (Coord3*)*valueList;

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

	//normal

	loc = (Coord3*)valueList[1];

	n1 = (*loc + sorted_y[0]);
	n2 = (*loc + sorted_y[1]);
	n3 = (*loc + sorted_y[2]);

	setNormal(*v1, *v2, *v3, *n1, *n2, *n3, normalTablePram);
	setZ(*v1, *v2, *v3, *n1, *n2, *n3, zTablePram);

	float midX = getMidX(*v1, *v2, *v3);
	float midY = getMidY(*v1, *v2, *v3);
	float midZ = getMidZ(*v1, *v2, *v3);

	float midNx = getNormalX(midX, midY, midZ, normalTablePram);
	float midNy = getNormalY(midX, midY, midZ, normalTablePram);
	float midNz = getNormalZ(midX, midY, midZ, normalTablePram);

	//need normalize
	float var;
	float height_scale = rand_height();
	float base_scale = rand_base();
	var = pow(pow(midNx, 2) + pow(midNy, 2) + pow(midNz, 2), 0.5);
	midNx /= var;
	midNy /= var;
	midNz /= var;


	GzCoord top = { midX + midNx * height_scale, midY + midNy * height_scale, midZ + midNz * height_scale };


	float dX, dY, dZ;
	dX = (*v1)[0] - midX;
	dY = (*v1)[1] - midY;
	dZ = (*v1)[2] - midZ;
	var = pow(pow(dX, 2) + pow(dY, 2) + pow(dZ, 2), 0.5);

	dX /= var;
	dY /= var;
	dZ /= var;
	GzCoord base = { dX*base_scale, dY*base_scale, dZ*base_scale };
	GzCoord N = { midNx, midNy, midNz };
	GzCoord mid = { midX, midY, midZ };

	GzCoord base1, base2;

	rotation(base, mid, N, &base1, &base2);

	calPlaneParams(base1, base2, top, &A, &B, &C, &D, &flag);

	midNx = A;
	midNy = B;
	midNz = C;
	var = pow(pow(midNx, 2) + pow(midNy, 2) + pow(midNz, 2), 0.5);
	midNx /= var;
	midNy /= var;
	midNz /= var;

	GzTextureIndex uvTop = { 0.5,0 };
	GzTextureIndex uvBase1 = { -0.5, 0 };
	GzTextureIndex uvBase2 = { 0, 1 };

	// get the 1/3 of the triangle
	GzCoord L32, L31, R32, R31;
	dX = top[0] - base1[0];
	dY = top[1] - base1[1];
	dZ = top[2] - base1[2];
	float c = rand_curve();
	float curve[5] = { 0,c / 200,c / 20,c / 10,c };
	float sign = rand_sign();

	L32[0] = dX * 1 / 3 + base1[0] + midNx * curve[0] * sign;
	L32[1] = dY * 1 / 3 + base1[1] + midNy * curve[0] * sign;
	L32[2] = dZ * 1 / 3 + base1[2] + midNz * curve[0] * sign;
	L31[0] = dX * 2 / 3 + base1[0] + midNx * curve[2] * sign;
	L31[1] = dY * 2 / 3 + base1[1] + midNy * curve[2] * sign;
	L31[2] = dZ * 2 / 3 + base1[2] + midNz * curve[2] * sign;
	dX = top[0] - base2[0];
	dY = top[1] - base2[1];
	dZ = top[2] - base2[2];
	R32[0] = dX * 1 / 3 + base2[0] + midNx * curve[1] * sign;
	R32[1] = dY * 1 / 3 + base2[1] + midNy * curve[1] * sign;
	R32[2] = dZ * 1 / 3 + base2[2] + midNz * curve[1] * sign;
	R31[0] = dX * 2 / 3 + base2[0] + midNx * curve[3] * sign;
	R31[1] = dY * 2 / 3 + base2[1] + midNy * curve[3] * sign;
	R31[2] = dZ * 2 / 3 + base2[2] + midNz * curve[3] * sign;
	top[0] += midNx * curve[4];
	top[1] += midNy * curve[4];
	top[2] += midNz * curve[4];

	// tri 1
	result[0][0] = L32[0];
	result[0][1] = L32[1];
	result[0][2] = L32[2];
	result[0][3] = midNx;
	result[0][4] = midNy;
	result[0][5] = midNz;
	result[0][6] = uvTop[0];
	result[0][7] = uvTop[1];

	result[1][0] = base1[0];
	result[1][1] = base1[1];
	result[1][2] = base1[2];
	result[1][3] = midNx;
	result[1][4] = midNy;
	result[1][5] = midNz;
	result[1][6] = uvBase1[0];
	result[1][7] = uvBase1[1];

	result[2][0] = base2[0];
	result[2][1] = base2[1];
	result[2][2] = base2[2];
	result[2][3] = midNx;
	result[2][4] = midNy;
	result[2][5] = midNz;
	result[2][6] = uvBase2[0];
	result[2][7] = uvBase2[1];
	//tri 2
	calPlaneParams(L32, base2, R32, &A, &B, &C, &D, &flag);

	midNx = A;
	midNy = B;
	midNz = C;
	var = pow(pow(midNx, 2) + pow(midNy, 2) + pow(midNz, 2), 0.5);
	midNx /= var;
	midNy /= var;
	midNz /= var;
	result[3][0] = L32[0];
	result[3][1] = L32[1];
	result[3][2] = L32[2];
	result[3][3] = midNx;
	result[3][4] = midNy;
	result[3][5] = midNz;
	result[3][6] = uvTop[0];
	result[3][7] = uvTop[1];

	result[4][0] = R32[0];
	result[4][1] = R32[1];
	result[4][2] = R32[2];
	result[4][3] = midNx;
	result[4][4] = midNy;
	result[4][5] = midNz;
	result[4][6] = uvBase1[0];
	result[4][7] = uvBase1[1];

	result[5][0] = base2[0];
	result[5][1] = base2[1];
	result[5][2] = base2[2];
	result[5][3] = midNx;
	result[5][4] = midNy;
	result[5][5] = midNz;
	result[5][6] = uvBase2[0];
	result[5][7] = uvBase2[1];
	//tri 3
	calPlaneParams(L32, L31, R32, &A, &B, &C, &D, &flag);

	midNx = A;
	midNy = B;
	midNz = C;
	var = pow(pow(midNx, 2) + pow(midNy, 2) + pow(midNz, 2), 0.5);
	midNx /= var;
	midNy /= var;
	midNz /= var;
	result[6][0] = L32[0];
	result[6][1] = L32[1];
	result[6][2] = L32[2];
	result[6][3] = midNx;
	result[6][4] = midNy;
	result[6][5] = midNz;
	result[6][6] = uvTop[0];
	result[6][7] = uvTop[1];

	result[7][0] = R32[0];
	result[7][1] = R32[1];
	result[7][2] = R32[2];
	result[7][3] = midNx;
	result[7][4] = midNy;
	result[7][5] = midNz;
	result[7][6] = uvBase1[0];
	result[7][7] = uvBase1[1];

	result[8][0] = L31[0];
	result[8][1] = L31[1];
	result[8][2] = L31[2];
	result[8][3] = midNx;
	result[8][4] = midNy;
	result[8][5] = midNz;
	result[8][6] = uvBase2[0];
	result[8][7] = uvBase2[1];

	//tri 4
	calPlaneParams(R31, L31, R32, &A, &B, &C, &D, &flag);

	midNx = A;
	midNy = B;
	midNz = C;
	var = pow(pow(midNx, 2) + pow(midNy, 2) + pow(midNz, 2), 0.5);
	midNx /= var;
	midNy /= var;
	midNz /= var;
	result[9][0] = R31[0];
	result[9][1] = R31[1];
	result[9][2] = R31[2];
	result[9][3] = midNx;
	result[9][4] = midNy;
	result[9][5] = midNz;
	result[9][6] = uvTop[0];
	result[9][7] = uvTop[1];

	result[10][0] = R32[0];
	result[10][1] = R32[1];
	result[10][2] = R32[2];
	result[10][3] = midNx;
	result[10][4] = midNy;
	result[10][5] = midNz;
	result[10][6] = uvBase1[0];
	result[10][7] = uvBase1[1];

	result[11][0] = L31[0];
	result[11][1] = L31[1];
	result[11][2] = L31[2];
	result[11][3] = midNx;
	result[11][4] = midNy;
	result[11][5] = midNz;
	result[11][6] = uvBase2[0];
	result[11][7] = uvBase2[1];
	//tri 5
	calPlaneParams(R31, L31, top, &A, &B, &C, &D, &flag);

	midNx = A;
	midNy = B;
	midNz = C;
	var = pow(pow(midNx, 2) + pow(midNy, 2) + pow(midNz, 2), 0.5);
	midNx /= var;
	midNy /= var;
	midNz /= var;
	result[12][0] = R31[0];
	result[12][1] = R31[1];
	result[12][2] = R31[2];
	result[12][3] = midNx;
	result[12][4] = midNy;
	result[12][5] = midNz;
	result[12][6] = uvTop[0];
	result[12][7] = uvTop[1];

	result[13][0] = top[0];
	result[13][1] = top[1];
	result[13][2] = top[2];
	result[13][3] = midNx;
	result[13][4] = midNy;
	result[13][5] = midNz;
	result[13][6] = uvBase1[0];
	result[13][7] = uvBase1[1];

	result[14][0] = L31[0];
	result[14][1] = L31[1];
	result[14][2] = L31[2];
	result[14][3] = midNx;
	result[14][4] = midNy;
	result[14][5] = midNz;
	result[14][6] = uvBase2[0];
	result[14][7] = uvBase2[1];
	return result;
}

void setNormal(GzCoord V1, GzCoord V2, GzCoord V3, GzCoord N1, GzCoord N2, GzCoord N3, double** normalTablePram) {

	GzCoord p1, p2, p3, N;
	for (int i = 0; i < 3; i++) {
		p1[i] = (V1)[i];
		p2[i] = (V2)[i];
		p3[i] = (V3)[i];
	}

	for (size_t i = 0; i < 3; i++)
	{
		p1[2] = (N1)[i];
		p2[2] = (N2)[i];
		p3[2] = (N3)[i];
		calPlaneParams(p1, p2, p3, &normalTablePram[i][0], &normalTablePram[i][1], &normalTablePram[i][2], &normalTablePram[i][3], &normalTablePram[i][4]);

	}
	return;
}

void setZ(GzCoord V1, GzCoord V2, GzCoord V3, GzCoord N1, GzCoord N2, GzCoord N3, double* zTablePram) {


	GzCoord p1, p2, p3, N;

	for (int i = 0; i < 3; i++) {
		p1[i] = (V1)[i];
		p2[i] = (V2)[i];
		p3[i] = (V3)[i];
	}

	calPlaneParams(p1, p2, p3, &zTablePram[0], &zTablePram[1], &zTablePram[2], &zTablePram[3], &zTablePram[4]);

	return;
}

float getNormalX(float x, float y, float z, double** normalTablePram) {
	return calZ(x, y, normalTablePram[0][0], normalTablePram[0][1], normalTablePram[0][2], normalTablePram[0][3], normalTablePram[0][4]);
}

float getNormalY(float x, float y, float z, double** normalTablePram) {
	return calZ(x, y, normalTablePram[1][0], normalTablePram[1][1], normalTablePram[1][2], normalTablePram[1][3], normalTablePram[1][4]);;
}

float getNormalZ(float x, float y, float z, double** normalTablePram) {
	return calZ(x, y, normalTablePram[2][0], normalTablePram[2][1], normalTablePram[2][2], normalTablePram[2][3], normalTablePram[2][4]);;
}

float getZ(float x, float y, float z, double* zTablePram) {
	return calZ(x, y, zTablePram[0], zTablePram[1], zTablePram[2], zTablePram[3], zTablePram[4]);
}

float getMidX(GzCoord V1, GzCoord V2, GzCoord V3) {
	return (V1[0] + V2[0] + V3[0]) / 3;
}

float getMidY(GzCoord V1, GzCoord V2, GzCoord V3) {
	return (V1[1] + V2[1] + V3[1]) / 3;
}

float getMidZ(GzCoord V1, GzCoord V2, GzCoord V3) {
	return (V1[2] + V2[2] + V3[2]) / 3;
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

void calPlaneParams(GzCoord v1, GzCoord v2, GzCoord v3, double* A, double* B, double* C, double* D, double* flag) {
	// calPlaneParams(p1,p2,p3,&A1,&B1,&C1,&D1);
	GzCoord norm;
	GzCoord e1, e2;
	if (v1[2] == v2[2] && v2[2] == v3[2]) {
		*flag = 0;
		*C = v2[2];

	}
	else {
		*flag = 1;
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

float calZ(float x, float y, double A, double B, double C, double D, double flag) {
	if (flag == 0) {
		return C;
	}
	if (C == 0) {
		C = 1;
	}
	return (-D - A * x - B * y) / C;
}


void swap(float& a, float& b)
{
	float c(std::move(a));
	a = std::move(b);
	b = std::move(c);
}

void sort_base_y(GzCoord v1, GzCoord v2, GzCoord v3) {

	if (v1[1] > v2[1]) {
		for (int i = 0; i < 3; i++) {
			swap(v1[i], v2[i]);
		}
	}

	if (v2[1] > v3[1]) {
		for (int i = 0; i < 3; i++) {
			swap(v2[i], v3[i]);
		}
	}

	if (v1[1] > v2[1]) {
		for (int i = 0; i < 3; i++) {
			swap(v1[i], v2[i]);
		}
	}
}

int getcolor(int numlights, float spec, GzLight* lights, GzLight ambientlight, GzColor Ka, GzColor Kd, GzColor Ks, GzCoord normal, GzColor color, int flag, GzColor texture) {

	GzCoord E = { 0, 0, -1 }, R;
	float nde, ndl, rde;
	GzColor diffuse_intens, specular_intens;
	for (int i = 0; i < 3; i++) {
		diffuse_intens[i] = 0;
		specular_intens[i] = 0;
	}

	nde = normal[0] * E[0] + normal[1] * E[1] + normal[2] * E[2];

	GzCoord newN = { 0, 0, 0 };

	for (int i = 0; i < numlights; i++) {
		ndl = normal[0] * lights[i].direction[0] + normal[1] * lights[i].direction[1] + normal[2] * lights[i].direction[2];
		if (ndl > 0 && nde > 0) {

			for (int j = 0; j < 3; j++)
				R[j] = 2 * ndl * normal[j] - lights[i].direction[j];
		}

		else if (ndl < 0 && nde < 0) {

			for (int j = 0; j < 3; j++)
				newN[j] = -1 * normal[j];

			ndl = newN[0] * lights[i].direction[0] + newN[1] * (lights[i].direction[1]) + newN[2] * (lights[i].direction[2]);
			for (int j = 0; j < 3; j++)
				R[j] = 2 * (ndl)* newN[j] - lights[i].direction[j];
		}
		else
		{
			//skip
			continue;
		}

		float modR = sqrt((R[0] * R[0]) + (R[1] * R[1]) + (R[2] * R[2]));

		for (int j = 0; j < 3; j++)
			R[j] = R[j] / modR;

		rde = R[0] * E[0] + R[1] * E[1] + R[2] * E[2];

		if (rde < 0) {
			rde = 0;
		}
		else if (rde > 1) {
			rde = 1;
		}

		for (int j = 0; j < 3; j++) {
			diffuse_intens[j] += ndl * lights[i].color[j];
			specular_intens[j] += pow(rde, spec) * lights[i].color[j];
		}
	}



	//color = specular + diffuse + ambient;
	if (flag == 0) {
		for (int i = 0; i < 3; i++)
			color[i] = Ks[i] * specular_intens[i] + Kd[i] * diffuse_intens[i] + Ka[i] * ambientlight.color[i];
	}
	else {
		for (int i = 0; i < 3; i++)
			color[i] = Ks[i] * specular_intens[i] + texture[i] * diffuse_intens[i] + texture[i] * ambientlight.color[i];
	}
	return GZ_SUCCESS;
}

