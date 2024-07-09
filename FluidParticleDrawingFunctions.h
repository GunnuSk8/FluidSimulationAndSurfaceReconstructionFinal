#pragma once
#include "FluidParticlesComputations.h"
#include "FluidParticlesComputations2D.h"
#include "FluidParticlesComputations3D.h"
#include "happly.h"
#include<cmath>
static void InitialiseParticles(std::vector<particle> &particles, int &totalNumberOfParticles) {

	float outerBounds = 6.9f;
	float offSet = 0;
	float k = -0.5f, j = -0.5f, particleSpacing=0.001f;
	for (int i = 0; i < totalNumberOfParticles; i++) {
		 
		particles[i].p.x = (float)((double)std::rand() / (double)RAND_MAX - 0.5) * outerBounds / 5 + offSet;
		particles[i].p.y = (float)((double)std::rand() / (double)RAND_MAX - 0.5) * outerBounds / 5 + offSet;
	
		particles[i].particleProperty = ExampleFunction(particles[i].p);
		particles[i].particleIndex = i;
		
	}
	
}
happly::PLYData plyIn("combined_mesh.ply");

// Get mesh-style data from the object
std::vector<std::array<double, 3>> vPos = plyIn.getVertexPositions();
inline const float scale = 0.004f;
inline int numberOfParticles = (int)vPos.size();
static void InitialisePolytopeParticles(std::vector<particle> &particles) {

	for (int i = 0; i<numberOfParticles;i++) {
		particles[i].p.x = (float)vPos[i][0]* scale;
		particles[i].p.y = (float)vPos[i][1]* scale;
		particles[i].p.z = (float)vPos[i][2]* scale;
		particles[i].particleIndex = i;
		centroidPolytope.x += particles[i].p.x;
		centroidPolytope.y += particles[i].p.y;
		centroidPolytope.z += particles[i].p.z;
	}
	centroidPolytope.x /= numberOfParticles;
	centroidPolytope.y /= numberOfParticles;
	centroidPolytope.z /= numberOfParticles;
	for (int i = 0; i < numberOfParticles; i++) {
		particles[i].p.x -= centroidPolytope.x;
		particles[i].p.y -= centroidPolytope.y;
		particles[i].p.z -= centroidPolytope.z;
	}
}


static void InitialiseParticles3D(std::vector<particle>& particles, int& totalNumberOfParticles) {


	float outerBounds = 3.9f;
	float offSet = 0;
	for (int i = 0; i < totalNumberOfParticles; i++) {
		

	
		particles[i].p.x = (float)((double)std::rand() / (double)RAND_MAX - 0.5) * outerBounds / 5 + offSet;
		particles[i].p.y = (float)((double)std::rand() / (double)RAND_MAX - 0.5) * outerBounds / 5 + offSet;
		particles[i].p.z = (float)((double)std::rand() / (double)RAND_MAX - 0.5) * outerBounds / 5 + offSet;
		
		particles[i].particleProperty = ExampleFunction(particles[i].p);
		particles[i].particleIndex = i;

	}

}
static float radiusSphere = 0.5;
static void InitialiseParticles3DGradientTest(std::vector<particle> &particles) {
	
	for (int i = 0; i < particles.size(); i++) {
		float phi = (float)((double)std::rand() / (double)RAND_MAX - 0.5) * 3.14159f;
		float theta = (float)((double)std::rand() / (double)RAND_MAX) * 2 * 3.14159f;

		particles[i].p.x =  radiusSphere * cos(phi) * cos(theta);
		particles[i].p.y =  radiusSphere * cos(phi) * sin(theta);
		particles[i].p.z =  radiusSphere * sin(phi);
	
		particles[i].particleIndex = i;

	}
}

static void DrawParticles(std::vector<particle> &particles,  int &totalNumberOfParticles) {
	/*std::for_each(std::execution::par, particles.begin(), particles.end(),
		[](particle& par) {
			DrawParticle(par);
		});*/
	//here sequential rendering is showing better looking results at least when
	//the particles fall under gravity without inter particle interaction
	for (int i = 0; i < particles.size(); i++) {
		//DrawPoint(particles[i].p, 0, 0, 1,0.5);
		//DrawCircleBoundary(particles[i].p,smoothingRadius,1,1,1,1);
		//DrawCircleBoundary(particles[i].p,smoothingRadius*2,1,0,0,1);
		
		//updated kernel visualisation
		DrawPoint(particles[i].updatedKernelCenter, 1, 0, 1, 0.5);
		glLineWidth(1.0);
		glColor4f(1,1,1,0.3f);
		glBegin(GL_LINES);
		glVertex3f(particles[i].p.x, particles[i].p.y, particles[i].p.z);
		glVertex3f(particles[i].updatedKernelCenter.x, particles[i].updatedKernelCenter.y, particles[i].updatedKernelCenter.z);
		glEnd();
		
		DrawParticle(particles[i]);

		//DrawCircle(particles[i].p, smoothingRadius);
		
			
		
	}

}

static void DrawNbhdParticlesForParticle(particle &par, std::vector<particleIndexCellKeyPair>& spatialLookup, std::vector<int> startIndices,
	std::vector<particle> &particles,int &totalNumberOfParticles, float &smoothingRadius) {
	std::vector<int> nbhdParticleIndexes = GetNbhdParticlesForParticle(par,spatialLookup,startIndices,particles,totalNumberOfParticles,smoothingRadius);
	for (int i = 0; i < nbhdParticleIndexes.size();i++) {
		DrawPoint(particles[nbhdParticleIndexes[i]].p,0,1,0,1);
	}
	DrawPoint(par.p,1,0,0,1);
	//DrawCircleBoundary(par.p,smoothingRadius,1,1,1,1);
}

static void DrawNbhdParticlesForParticle3D(particle& par, std::vector<particleIndexCellKeyPair>& spatialLookup, std::vector<int> startIndices,
	std::vector<particle>& particles, int& totalNumberOfParticles, float& smoothingRadius) {
	std::vector<int> nbhdParticleIndexes = GetNbhdParticlesForParticle3D(par, spatialLookup, startIndices, particles, totalNumberOfParticles, smoothingRadius);
	std::cout << nbhdParticleIndexes.size()<<std::endl;
	for (int i = 0; i < nbhdParticleIndexes.size(); i++) {
		DrawPoint(particles[nbhdParticleIndexes[i]].p, 0, 1, 0, 1);
	}
	DrawPoint(par.p, 1, 0, 0, 1,2);
	DrawSphere(par.p);
	std::cout << par.Gi << std::endl;
	//DrawCircleBoundary(par.p, smoothingRadius, 1, 1, 1, 1);
}