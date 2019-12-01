// Application5.cpp: implementation of the Application5 class.
//
//////////////////////////////////////////////////////////////////////

/*
 * application test code for homework assignment #5
*/

#include "stdafx.h"
#include "CS580HW.h"
#include "Application5.h"
#include "Gz.h"
#include "rend.h"
#include <windows.h>

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#define new DEBUG_NEW
#endif

#define INFILE  "ppot.asc"
#define OUTFILE "output.ppm"
#define OUTFILE2 "output_otherround.ppm"


extern int tex_fun(float u, float v, float w, GzColor color); /* image texture function */
extern int ptex_fun(float u, float v, float w, GzColor color); /* procedural texture function */
extern int GzFreeTexture();
void recursive(GzCoord* vertexList, GzCoord* normalList, GzTextureIndex* uvList, float* midPoint, GzRender* m_pRender, GzToken* nameListTriangle, int count);
void shade(GzCoord norm, GzCoord color);

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Application5::Application5()
{

}

Application5::~Application5()
{
	Clean();
}

int Application5::Initialize()
{
	GzCamera	camera;
	int		    xRes, yRes;	/* display parameters */

	GzToken		nameListShader[9]; 	    /* shader attribute names */
	GzPointer   valueListShader[9];		/* shader attribute pointers */
	GzToken     nameListLights[10];		/* light info */
	GzPointer   valueListLights[10];
	int			shaderType, interpStyle;
	float		specpower;
	int		status;

	status = 0;

	/*
	 * Allocate memory for user input
	 */
	m_pUserInput = new GzInput;

	/*
	 * initialize the display and the renderer
	 */
	m_nWidth = 900;		// frame buffer and display width
	m_nHeight = 900;    // frame buffer and display height

	m_pRender = new GzRender(m_nWidth, m_nHeight);
	m_pRender->GzDefault();

	m_pFrameBuffer = m_pRender->framebuffer;

	/* Translation matrix */
	GzMatrix	scale =
	{
		3.25,	0.0,	0.0,	0.0,
		0.0,	3.25,	0.0,	-3.25,
		0.0,	0.0,	3.25,	3.5,
		0.0,	0.0,	0.0,	1.0
	};

	GzMatrix	rotateX =
	{
		1.0,	0.0,	0.0,	0.0,
		0.0,	.7071,	.7071,	0.0,
		0.0,	-.7071,	.7071,	0.0,
		0.0,	0.0,	0.0,	1.0
	};

	GzMatrix	rotateY =
	{
		.866,	0.0,	-0.5,	0.0,
		0.0,	1.0,	0.0,	0.0,
		0.5,	0.0,	.866,	0.0,
		0.0,	0.0,	0.0,	1.0
	};

#if 1 	/* set up app-defined camera if desired, else use camera defaults */
	camera.position[X] = -3;
	camera.position[Y] = -25;
	camera.position[Z] = -4;

	camera.lookat[X] = 7.8;
	camera.lookat[Y] = 0.7;
	camera.lookat[Z] = 6.5;

	camera.worldup[X] = -0.2;
	camera.worldup[Y] = 1.0;
	camera.worldup[Z] = 0.0;

	camera.FOV = 63.7;              /* degrees *              /* degrees */

	status |= m_pRender->GzPutCamera(camera);
#endif 

	/* Start Renderer */
	status |= m_pRender->GzBeginRender();

	/* Light */
	GzLight	light1 = { {-0.7071, 0.7071, 0}, {0.5, 0.5, 0.9} };
	GzLight	light2 = { {0, -0.7071, -0.7071}, {0.9, 0.2, 0.3} };
	GzLight	light3 = { {0.7071, 0.0, -0.7071}, {0.2, 0.7, 0.3} };
	GzLight	ambientlight = { {0, 0, 0}, {0.5, 0.5, 0.5} };

	/* Material property */
	GzColor specularCoefficient = { 0.3, 0.3, 0.3 };
	GzColor ambientCoefficient = { 0.1, 0.5, 0.1 };
	GzColor diffuseCoefficient = { 0.1, 0.5, 0.1 };

	/*
	  renderer is ready for frame --- define lights and shader at start of frame
	*/

	/*
	 * Tokens associated with light parameters
	 */
	nameListLights[0] = GZ_DIRECTIONAL_LIGHT;
	valueListLights[0] = (GzPointer)&light1;
	nameListLights[1] = GZ_DIRECTIONAL_LIGHT;
	valueListLights[1] = (GzPointer)&light2;
	nameListLights[2] = GZ_DIRECTIONAL_LIGHT;
	valueListLights[2] = (GzPointer)&light3;
	status |= m_pRender->GzPutAttribute(3, nameListLights, valueListLights);

	nameListLights[0] = GZ_AMBIENT_LIGHT;
	valueListLights[0] = (GzPointer)&ambientlight;
	status |= m_pRender->GzPutAttribute(1, nameListLights, valueListLights);

	/*
	 * Tokens associated with shading
	 */
	nameListShader[0] = GZ_DIFFUSE_COEFFICIENT;
	valueListShader[0] = (GzPointer)diffuseCoefficient;

	/*
	* Select either GZ_COLOR or GZ_NORMALS as interpolation mode
	*/
	nameListShader[1] = GZ_INTERPOLATE;
	//interpStyle = GZ_COLOR;         /* Gouraud shading */
	interpStyle = GZ_NORMALS;         /* Phong shading */
	valueListShader[1] = (GzPointer)&interpStyle;

	nameListShader[2] = GZ_AMBIENT_COEFFICIENT;
	valueListShader[2] = (GzPointer)ambientCoefficient;
	nameListShader[3] = GZ_SPECULAR_COEFFICIENT;
	valueListShader[3] = (GzPointer)specularCoefficient;
	nameListShader[4] = GZ_DISTRIBUTION_COEFFICIENT;
	specpower = 32;
	valueListShader[4] = (GzPointer)&specpower;

	nameListShader[5] = GZ_TEXTURE_MAP;
#if 0  /* set up null texture function or valid pointer */
	valueListShader[5] = (GzPointer)0;
#else
	valueListShader[5] = (GzPointer)(ptex_fun);	/* or use ptex_fun */
#endif
	status |= m_pRender->GzPutAttribute(6, nameListShader, valueListShader);


	status |= m_pRender->GzPushMatrix(scale);
	status |= m_pRender->GzPushMatrix(rotateY);
	status |= m_pRender->GzPushMatrix(rotateX);

	if (status) exit(GZ_FAILURE);

	if (status)
		return(GZ_FAILURE);
	else
		return(GZ_SUCCESS);
}

int Application5::Render()
{
	GzToken		nameListTriangle[3]; 	/* vertex attribute names */

	GzPointer	valueListTriangle[3]; 	/* vertex attribute pointers */

	GzPointer*	valueListGrassTriangle[5];
	GzToken		nameListShader[9]; 	    /* shader attribute names */
	GzPointer   valueListShader[9];		/* shader attribute pointers */

	GzCoord		vertexList[3];	/* vertex position coordinates */
	GzCoord		normalList[3];	/* vertex normals */
	GzTextureIndex uvList[3];		/* vertex texture map indices */

	GzCoord*		vertexGrassList[5];	/* vertex position coordinates */
	GzCoord*		normalGrassList[5];	/* vertex normals */
	GzTextureIndex* uvGrassList[5];		/* vertex texture map indices */

	char		dummy[256];
	int			status;


	/* Initialize Display */
	status |= m_pRender->GzDefault();  /* init for new frame */

	/*
	* Tokens associated with triangle vertex values
	*/
	nameListTriangle[0] = GZ_POSITION;
	nameListTriangle[1] = GZ_NORMAL;
	nameListTriangle[2] = GZ_TEXTURE_INDEX;

	// I/O File open
	FILE *infile;
	if ((infile = fopen(INFILE, "r")) == NULL)
	{
		AfxMessageBox("The input file was not opened\n");
		return GZ_FAILURE;
	}

	FILE *outfile;
	FILE *outfile2;
	if ((outfile = fopen(OUTFILE, "wb")) == NULL)
	{
		AfxMessageBox("The output file was not opened\n");
		return GZ_FAILURE;
	}
	if ((outfile2 = fopen(OUTFILE2, "wb")) == NULL)
	{
		AfxMessageBox("The output file was not opened\n");
		return GZ_FAILURE;
	}

	/*
	* Walk through the list of triangles, set color
	* and render each triangle
	*/
	while (fscanf(infile, "%s", dummy) == 1) { 	/* read in tri word */
		fscanf(infile, "%f %f %f %f %f %f %f %f",
			&(vertexList[0][0]), &(vertexList[0][1]),
			&(vertexList[0][2]),
			&(normalList[0][0]), &(normalList[0][1]),
			&(normalList[0][2]),
			&(uvList[0][0]), &(uvList[0][1]));
		fscanf(infile, "%f %f %f %f %f %f %f %f",
			&(vertexList[1][0]), &(vertexList[1][1]),
			&(vertexList[1][2]),
			&(normalList[1][0]), &(normalList[1][1]),
			&(normalList[1][2]),
			&(uvList[1][0]), &(uvList[1][1]));
		fscanf(infile, "%f %f %f %f %f %f %f %f",
			&(vertexList[2][0]), &(vertexList[2][1]),
			&(vertexList[2][2]),
			&(normalList[2][0]), &(normalList[2][1]),
			&(normalList[2][2]),
			&(uvList[2][0]), &(uvList[2][1]));

		/*
		 * Set the value pointers to the first vertex of the
		 * triangle, then feed it to the renderer
		 * NOTE: this sequence matches the nameList token sequence
		 */

		uvList[0][2] = 0;
		uvList[1][2] = 0;
		uvList[2][2] = 0;

		//get first midpoint
		valueListTriangle[0] = (GzPointer)vertexList;
		valueListTriangle[1] = (GzPointer)normalList;
		valueListTriangle[2] = (GzPointer)uvList;

		float *midPoint = (float *)malloc(8 * sizeof(float));
		for (int i = 0; i < 5; i++) {
			vertexGrassList[i] = (GzCoord *)malloc(3 * sizeof(GzCoord));
			normalGrassList[i] = (GzCoord *)malloc(3 * sizeof(GzCoord));
			uvGrassList[i] = (GzTextureIndex*)malloc(3 * sizeof(GzTextureIndex));
			valueListGrassTriangle[i] = (GzPointer *)malloc(3 * sizeof(GzPointer));
		}

		m_pRender->GzAddGrassWithModelSpace(3, nameListTriangle, valueListTriangle, vertexGrassList, normalGrassList, uvGrassList, midPoint);

		for (int i = 0; i < 5; i++) {
			valueListGrassTriangle[i][0] = (GzPointer)vertexGrassList[i];
			valueListGrassTriangle[i][1] = (GzPointer)normalGrassList[i];
			valueListGrassTriangle[i][2] = (GzPointer)uvGrassList[i];
			m_pRender->GzPutTriangle(3, nameListTriangle, valueListGrassTriangle[i]);
		}


		//make each point interact with midpoint
		int count = 0;
		recursive(vertexList, normalList, uvList, midPoint, m_pRender, nameListTriangle, count);
		free(midPoint);

		m_pRender->GzPutTriangle(3, nameListTriangle, valueListTriangle);

	}

	m_pRender->GzFlushDisplay2File(outfile); 	/* write out or update display to file*/
	m_pRender->GzFlushDisplay2FrameBuffer();	// write out or update display to frame buffer


	if (fclose(infile))
		AfxMessageBox(_T("The input file was not closed\n"));

	if (fclose(outfile))
		AfxMessageBox(_T("The output file was not closed\n"));
	if (fclose(outfile2))
		AfxMessageBox(_T("The output file was not closed\n"));

	if (status)
		return(GZ_FAILURE);
	else
		return(GZ_SUCCESS);
}


float distance(float x1, float y1, float z1, float x2, float y2, float z2) {
	return pow(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2), 0.5);
}

void recursive(GzCoord* vertexList, GzCoord* normalList, GzTextureIndex* uvList, float* midPoint, GzRender* m_pRender, GzToken* nameListTriangle, int count) {
	if (count == 4 || distance(vertexList[0][0], vertexList[0][1], vertexList[0][2], midPoint[0], midPoint[1], midPoint[2]) < 0.1) return;
	GzPointer TempListTriangle[3];
	GzCoord		TempVertexList[3];
	GzCoord		TempNormalList[3];
	GzTextureIndex TempUVList[3];

	GzPointer	valueListGrassTriangle[5][3];

	GzCoord*		vertexGrassList[5];
	GzCoord*		normalGrassList[5];
	GzTextureIndex* uvGrassList[5];

	for (int i = 0; i < 5; i++) {
		vertexGrassList[i] = (GzCoord *)malloc(3 * sizeof(GzCoord));
		normalGrassList[i] = (GzCoord *)malloc(3 * sizeof(GzCoord));
		uvGrassList[i] = (GzTextureIndex*)malloc(3 * sizeof(GzTextureIndex));
	}

	for (int k = 0; k < 3; k++) {

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++) {
				TempVertexList[i][j] = vertexList[i][j];
				TempNormalList[i][j] = normalList[i][j];
				TempUVList[i][j] = uvList[i][j];
			}

		}

		for (int i = 0; i < 3; i++)
		{
			TempVertexList[k][i] = midPoint[i];
			TempNormalList[k][i] = midPoint[i + 3];
			TempUVList[k][i] = midPoint[i + 5];
		}

		TempListTriangle[0] = (GzPointer)TempVertexList;
		TempListTriangle[1] = (GzPointer)TempNormalList;
		TempListTriangle[2] = (GzPointer)TempUVList;

		m_pRender->GzAddGrassWithModelSpace(3, nameListTriangle, TempListTriangle, vertexGrassList, normalGrassList, uvGrassList, midPoint);

		for (int i = 0; i < 5; i++) {
			valueListGrassTriangle[i][0] = (GzPointer)vertexGrassList[i];
			valueListGrassTriangle[i][1] = (GzPointer)normalGrassList[i];
			valueListGrassTriangle[i][2] = (GzPointer)uvGrassList[i];
			m_pRender->GzPutTriangle(3, nameListTriangle, valueListGrassTriangle[i]);
		}

		recursive(TempVertexList, TempNormalList, TempUVList, midPoint, m_pRender, nameListTriangle, count + 1);
	}

}

int Application5::Clean()
{
	/*
	 * Clean up and exit
	 */
	int	status = 0;

	free(m_pRender);
	status |= GzFreeTexture();

	if (status)
		return(GZ_FAILURE);
	else
		return(GZ_SUCCESS);
}



