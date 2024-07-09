#pragma once
#include<iostream>
#include "GLFW/glfw3.h"
#include "FluidParticlesProperties&Init.h"
#include "DrawingHelpingFunc.h"
#include "SurfaceReconstructionComputations.h"
#include<execution>
#include "SurfaceReconstructionComputationsForAllParticles.h"
#include<algorithm>


inline int rez2D = 120;
inline float isoValue2D = 0;

static std::vector<std::vector<point>> gridPt2D(rez2D, std::vector<point>(rez2D));

//bottom left, bottom right, top right, top left
constexpr static int getState(float blf, float brf, float trf, float tlf) {
	int bl=0, br=0, tr=0,tl=0;
	if (blf > isoValue2D) bl = 1;
	if (brf > isoValue2D) br = 1;
	if (trf > isoValue2D) tr = 1;
	if (tlf > isoValue2D) tl = 1;
	return bl * 8 + br * 4 + tr * 2 + tl * 1;
}

static point InterpolationPoint2D(point a, point b) {
	point ip;
	ip.x = a.x + (isovalue - a.intensity) * (b.x - a.x) / (b.intensity - a.intensity);
	ip.y = a.y + (isovalue - a.intensity) * (b.y - a.y) / (b.intensity - a.intensity);
	
	return ip;
}

static void Grid2DIntensityUpdate( std::vector<particleIndexCellKeyPair>& spatialLookup,
	std::vector<int> startIndices, std::vector<particle>& particles, int& totalNumberOfParticles, float& smoothingRadius) {

	std::for_each(std::execution::par, gridPt2D.begin(), gridPt2D.end(), [&particles, &spatialLookup,&startIndices,
		&totalNumberOfParticles,&smoothingRadius](std::vector<point>& ptrow) {

		std::for_each(std::execution::par, ptrow.begin(), ptrow.end(), [&particles, &spatialLookup, &startIndices,
		&totalNumberOfParticles, &smoothingRadius](point &pt) {

			pt.intensity = ScalarFieldAtPosition2D(pt,spatialLookup,startIndices,particles,totalNumberOfParticles,smoothingRadius);
			});
		});
	/*for (int x = 0; x < rez2D; x++) {
		for (int y = 0; y < rez2D; y++) {
			gridPt2D[x][y].intensity = ScalarFieldAtPosition2D(gridPt2D[x][y],particles);

		}
	}*/
}

static void Grid2DInit() {
	for (int x = 0; x < rez2D; x++) {
		for (int y = 0; y < rez2D; y++) {
			gridPt2D[x][y].x = ((float)x / (float)rez2D - 0.5f) * 2;
			gridPt2D[x][y].y = ((float)y / (float)rez2D - 0.5f) * 2;
			//gridPt2D[x][y].intensity = ((((double)rand()) / RAND_MAX)-0.5);
			//gridPt2D[x][y].intensity = function(gridPt2D[x][y].x, gridPt2D[x][y].y, gridPt2D[x][y].z);

		}
	}
	

}


static void MarchingSquares()
{
	

	//printing the grid according to the intensity
	//cant draw in parallel, because parallel drawing does not work -> emperical observation
	for (int i = 0; i < rez2D; i++) {
		for (int ii = 0; ii < rez2D; ii++) {
			glPointSize(3.0f);
			float ratio = 20;
			glColor3f(gridPt2D[i][ii].intensity/ratio, gridPt2D[i][ii].intensity/ratio, gridPt2D[i][ii].intensity/ratio);
			glBegin(GL_POINTS);
			glVertex2f(gridPt2D[i][ii].x, gridPt2D[i][ii].y);
			glEnd();
		}
	}


	//main marching squares algorithm
	for (int x = 0; x < rez2D-1 ; x++) {
		for (int y = 0; y < rez2D-1 ; y++) {
			float xCoord = gridPt2D[x][y].x;
			float yCoord = gridPt2D[x][y].y;
			
			//interpolation does not work properly
			/*point a = InterpolationPoint2D(gridPt2D[x][y],     gridPt2D[x+1][y]), 
				  b = InterpolationPoint2D(gridPt2D[x+1][y],   gridPt2D[x + 1][y+1]), 
				  c = InterpolationPoint2D(gridPt2D[x+1][y+1], gridPt2D[x][y+1]), 
				  d = InterpolationPoint2D(gridPt2D[x][y],     gridPt2D[x][y+1]);*/
			point a, b, c, d;
			a.x = xCoord + 1 / ((float)rez2D);
			a.y = yCoord;

			b.x = xCoord + 2 / ((float)rez2D);
			b.y = yCoord + 1 / ((float)rez2D);

			c.x = xCoord + 1 / ((float)rez2D);
			c.y = yCoord + 2 / ((float)rez2D);

			d.x = xCoord;
			d.y = yCoord + 1 / ((float)rez2D);

			/*if (gridPt2D[x][y].intensity > 0) {
				DrawPoint(a,1,0,0,1);
				DrawPoint(b,0,1,0,1);
				DrawPoint(c,0,0,1,1);
				DrawPoint(d,0.51,0.51,0.51,1);
			}*/


			int state = getState(gridPt2D[x][y].intensity, gridPt2D[x+1][y].intensity, gridPt2D[x + 1][y + 1].intensity, gridPt2D[x][y+1].intensity);
			
			switch (state)
			{
			case 1:
				DrawLine(d, c);
				break;
			case 2:
				DrawLine(c, b);
				break;
			case 3:
				DrawLine(d, b);
				break;
			case 4:
				DrawLine(a, b);
				break;
			case 5:
				DrawLine(a, b);
				DrawLine(c, d);
				break;
			case 6:
				DrawLine(a, c);
				break;
			case 7:
				DrawLine(a, d);
				break;
			case 8:
				DrawLine(a, d);
				break;
			case 9:
				DrawLine(a, c);
				break;
			case 10:
				DrawLine(a, d);
				DrawLine(b, c);
				break;
			case 11:
				DrawLine(a, b);
				break;
			case 12:
				DrawLine(d, b);
				break;
			case 13:
				DrawLine(c, b);
				break;
			case 14:
				DrawLine(c, d);
				break;
			case 15:
				break;
			default:
				break;
			}

		}
	}

}


