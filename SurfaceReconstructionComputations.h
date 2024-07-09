#pragma once
#include "FluidParticlesProperties&Init.h"
#include "MathHelpingFunc.h"
#include<algorithm>
#include<execution>
#include<vector>
#include"FluidParticlesComputations.h"
#include "FluidParticlesComputationsForAllParticles.h"

inline float lambda = 0.9f;
inline float riScale = 2.05;
inline float ri = riScale * smoothingRadius; 
//i think for all these values we have to use updatekernel centres no?????
//Wij
static float IsotropicWeightingFunction(Vector3D &xi, Vector3D &xj, float &radius) {
	float distance = DistanceBetweenPoints(xi, xj);
	if (distance > radius) return 0;
	float ratio = distance / radius;
	return 1 - Power(ratio, 3);
}
//summationWij
static float GetSummationWij(particle &par, std::vector<particle> &particles) {
	float summationWij = 0;
	for (int j = 0; j < particles.size(); j++) {
		summationWij += IsotropicWeightingFunction(par.p, particles[j].p, ri);
	}
	return summationWij;
}
//summationWijTimesXj
static Vector3D GetSummationWijTimesXj(particle& par, std::vector<particle>& particles) {
	Vector3D summationWijTimesXj;
	for (int j = 0; j < particles.size(); j++) {
		summationWijTimesXj.x += IsotropicWeightingFunction(par.p, particles[j].p, ri) * particles[j].p.x;
		summationWijTimesXj.y += IsotropicWeightingFunction(par.p, particles[j].p, ri) * particles[j].p.y;
		summationWijTimesXj.z += IsotropicWeightingFunction(par.p, particles[j].p, ri) * particles[j].p.z;
	}
	return summationWijTimesXj;
}
//XiW
static Vector3D GetXiW(particle &par, std::vector<particle> &particles) {
	Vector3D summationWijtimesXij = GetSummationWijTimesXj(par,particles);
	float summationWij = GetSummationWij(par,particles);

	Vector3D XiW;
	XiW.x = summationWijtimesXij.x / summationWij;
	XiW.y = summationWijtimesXij.y / summationWij;
	XiW.z = summationWijtimesXij.z / summationWij;

	return XiW;
}












