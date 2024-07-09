#pragma once
#include "SurfaceReconstructionComputations.h"

static void CalculateCi2D(particle& par, std::vector<particle>& particles) {
	Vector3D XiW = GetXiW(par, particles);
	float summationWij = 0;
	for (int j = 0; j < particles.size(); j++) {
		Vector3D XjMinusXiW;
		XjMinusXiW.x = particles[j].p.x - XiW.x;
		XjMinusXiW.y = particles[j].p.y - XiW.y;

		float Wij = IsotropicWeightingFunction(par.p, particles[j].p, ri);
		summationWij += Wij;
		par.Ci2D(0, 0) += Wij * XjMinusXiW.x * XjMinusXiW.x; par.Ci2D(0, 1) += Wij * XjMinusXiW.x * XjMinusXiW.y;
		par.Ci2D(1, 0) += Wij * XjMinusXiW.y * XjMinusXiW.x; par.Ci2D(1, 1) += Wij * XjMinusXiW.y * XjMinusXiW.y;

	}

	par.Ci2D(0, 0) /= summationWij; par.Ci2D(0, 1) /= summationWij;
	par.Ci2D(1, 0) /= summationWij; par.Ci2D(1, 1) /= summationWij;

}
inline float Kr2D = 12;
inline float Ks2D = 5400;
inline float Kn2D = 5;
inline int Ne2D = 3;

static Eigen::Vector<float, 2> EditEigenValues2D(Eigen::Vector<float, 2>& singularValues, int& numNbhdParticles) {
	Eigen::Vector<float, 2> editedSingularValues;

	if (numNbhdParticles < Ne2D) {
		editedSingularValues.setConstant(Kn);
		return editedSingularValues;
	}
	float sigma1 = singularValues(0);
	editedSingularValues(0) = Ks2D * singularValues(0);
	editedSingularValues(1) = Ks2D * Max(singularValues(1), sigma1 / Kr2D);
	

	return editedSingularValues;
}

static void CalculateCi2DTilda(particle& par, std::vector<particleIndexCellKeyPair>& spatialLookup,
	std::vector<int> startIndices, std::vector<particle>& particles, int& totalNumberOfParticles, float& smoothingRadius) {
	Eigen::JacobiSVD<Eigen::Matrix<float, 2, 2>> svd(par.Ci2D, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix<float, 2, 2> U = svd.matrixU();
	Eigen::Matrix<float, 2, 2> V = svd.matrixV();
	Eigen::Vector<float, 2> singularValues = svd.singularValues();

	std::vector<int> nbhdParticles = GetNbhdParticlesForParticle(par, spatialLookup, startIndices, particles, totalNumberOfParticles, smoothingRadius);
	int numNbhdParticles = (int)nbhdParticles.size();
	Eigen::Vector<float, 2> singularValuesTilda = EditEigenValues2D(singularValues, numNbhdParticles);

	par.Ci2DTilda = U * singularValuesTilda.asDiagonal() * V.transpose();

}


static void CalculateGi2D(particle& par) {
	par.Gi2D = par.Ci2DTilda.inverse();
	par.Gi2D(0, 0) /= smoothingRadius; par.Gi2D(0, 1) /= smoothingRadius;
	par.Gi2D(1, 0) /= smoothingRadius; par.Gi2D(1, 1) /= smoothingRadius;
}

static int index = 69;
static float ScalarFieldAtPosition2D(point position, std::vector<particleIndexCellKeyPair>& spatialLookup,
	std::vector<int> startIndices, std::vector<particle>& particles, int& totalNumberOfParticles, float& smoothingRadius) {
	//std::cout << position.x << " , " << position.y;
	particle temp; temp.predictedPosition.x = position.x; temp.predictedPosition.y = position.y; temp.p.x = position.x; temp.p.y = position.y;
	//std::vector<int> nbhdParticles = GetNbhdParticlesForParticle(temp, spatialLookup, startIndices, particles, totalNumberOfParticles, smoothingRadius);
	float scalarField = 0;

	

	std::vector<int> nbhdParticleIndexes = GetNbhdParticlesForParticle(temp, spatialLookup, startIndices, particles, totalNumberOfParticles, smoothingRadius);

	//isotropic approach of density
	/*for (int i = 0; i < nbhdParticleIndexes.size(); i++) {
		float distance = DistanceBetweenPoints(temp.predictedPosition, particles[nbhdParticleIndexes[i]].predictedPosition);
		float influence = SmoothingKernel(smoothingRadius, distance);
		scalarField += mass * influence;
	}*/


	//getting scalar field due to only one particle
	float distance = DistanceBetweenPoints(particles[index].predictedPosition, temp.predictedPosition);
	//float influence = SmoothingKernel(smoothingRadius,distance);
	float influence = SmoothingKernel(distance, particles[index].Gi2D);// std::cout << particles[400].Gi2D << std::endl;
	scalarField= mass * influence;

	//ellipse try broski successfull babyboi, you bozo, this is just drawing the Gi2D for that particular particle not the freaking scalar field!!!
	/*Eigen::Matrix<float, 1, 2> pos;
	pos(0, 0) = position.x - particles[index].predictedPosition.x;
	pos(0, 1) = position.y - particles[index].predictedPosition.y;

	Eigen::MatrixXf answer = pos * particles[index].Gi2D * pos.transpose();
	scalarField = answer(0,0);*/
	//scalarField = 0;
	//for (int i = 0; i < nbhdParticleIndexes.size(); i++) {
	//	int particleIndex = nbhdParticleIndexes[i];
	//	float distance = DistanceBetweenPoints(particles[particleIndex].predictedPosition, temp.predictedPosition);
	//	//float influence = SmoothingKernel(smoothingRadius,distance);
	//	float influence = SmoothingKernel(distance, particles[particleIndex].Gi2D);// std::cout << particles[400].Gi2D << std::endl;
	//	scalarField = mass * influence;

	//	
	//}

	
	//scalarField = OptimisedCalculateDensityForParticle(temp, spatialLookup, startIndices, particles, totalNumberOfParticles, smoothingRadius);

	return scalarField;


}
