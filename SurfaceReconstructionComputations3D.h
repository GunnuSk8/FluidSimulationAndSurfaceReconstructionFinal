#pragma once
#include "SurfaceReconstructionComputations.h"

//Ci, the covariance matrix
static void CalculateCi(particle& par, std::vector<particle>& particles) {
	Vector3D XiW = GetXiW(par, particles);
	float summationWij = 0;
	for (int j = 0; j < particles.size(); j++) {
		Vector3D XjMinusXiW;
		XjMinusXiW.x = particles[j].p.x - XiW.x;
		XjMinusXiW.y = particles[j].p.y - XiW.y;
		XjMinusXiW.z = particles[j].p.z - XiW.z;

		float Wij = IsotropicWeightingFunction(par.p, particles[j].p, ri);
		summationWij += Wij;
		par.Ci(0, 0) += Wij * XjMinusXiW.x * XjMinusXiW.x; par.Ci(0, 1) += Wij * XjMinusXiW.x * XjMinusXiW.y; par.Ci(0, 2) += Wij * XjMinusXiW.x * XjMinusXiW.z;
		par.Ci(1, 0) += Wij * XjMinusXiW.y * XjMinusXiW.x; par.Ci(1, 1) += Wij * XjMinusXiW.y * XjMinusXiW.y; par.Ci(1, 2) += Wij * XjMinusXiW.y * XjMinusXiW.z;
		par.Ci(2, 0) += Wij * XjMinusXiW.z * XjMinusXiW.x; par.Ci(2, 1) += Wij * XjMinusXiW.z * XjMinusXiW.y; par.Ci(2, 2) += Wij * XjMinusXiW.z * XjMinusXiW.z;
	}
	par.Ci(0, 0) /= summationWij; par.Ci(0, 1) /= summationWij; par.Ci(0, 2) /= summationWij;
	par.Ci(1, 0) /= summationWij; par.Ci(1, 1) /= summationWij; par.Ci(1, 2) /= summationWij;
	par.Ci(2, 0) /= summationWij; par.Ci(2, 1) /= summationWij; par.Ci(2, 2) /= summationWij;

}

inline float Kr = 8;
inline float Ks = 8000;
inline float Kn = 10;
inline int Ne = 5;

static Eigen::Vector<float, 3> EditEigenValues(Eigen::Vector<float, 3>& singularValues, int& numNbhdParticles) {
	Eigen::Vector<float, 3> editedSingularValues;

	if (numNbhdParticles < Ne) {
		editedSingularValues.setConstant(Kn);
		return editedSingularValues;
	}
	float sigma1 = singularValues(0);
	editedSingularValues(0) = Ks * singularValues(0);
	editedSingularValues(1) = Ks * Max(singularValues(1), sigma1 / Kr);
	editedSingularValues(2) = Ks * Max(singularValues(2), sigma1 / Kr);

	return editedSingularValues;
}

static void CalculateCiTilda(particle& par, std::vector<particleIndexCellKeyPair>& spatialLookup,
	std::vector<int> startIndices, std::vector<particle>& particles, int& totalNumberOfParticles, float& smoothingRadius) {
	Eigen::JacobiSVD<Eigen::Matrix<float, 3, 3>> svd(par.Ci, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix<float, 3, 3> U = svd.matrixU();
	Eigen::Matrix<float, 3, 3> V = svd.matrixV();
	Eigen::Vector<float, 3> singularValues = svd.singularValues();

	//std::vector<int> nbhdParticles = GetNbhdParticlesForParticle3D(par, spatialLookup, startIndices, particles, totalNumberOfParticles, smoothingRadius);
	std::vector<int> nbhdParticles = GetNbhdPartilcesBruteForce(par, particles,smoothingRadius);
	int numNbhdParticles = (int)nbhdParticles.size();
	Eigen::Vector<float, 3> singularValuesTilda = EditEigenValues(singularValues, numNbhdParticles);

	par.CiTilda = U * singularValuesTilda.asDiagonal() * V.transpose();

}

//Gi
static void CalculateGi(particle& par) {
	//par.Gi = par.Ci.inverse();
	par.Gi = par.CiTilda.inverse();
	par.Gi(0, 0) /= smoothingRadius; par.Gi(0, 1) /= smoothingRadius; par.Gi(0, 2) /= smoothingRadius;
	par.Gi(1, 0) /= smoothingRadius; par.Gi(1, 1) /= smoothingRadius; par.Gi(1, 2) /= smoothingRadius;
	par.Gi(2, 0) /= smoothingRadius; par.Gi(2, 1) /= smoothingRadius; par.Gi(2, 2) /= smoothingRadius;
}
static int index3D = 0;
static float ScalarFieldAtPosition(point position, std::vector<particleIndexCellKeyPair>& spatialLookup,
	std::vector<int> startIndices, std::vector<particle>& particles3D, int& totalNumberOfParticles, float& smoothingRadius) {

	//std::cout << position.x << " , " << position.y << " , " << position.z ;
	particle temp; temp.predictedPosition.x = position.x; temp.predictedPosition.y = position.y; temp.predictedPosition.z = position.z;
	float scalarField = 0;
	//std::vector<int> nbhdParticles = GetNbhdParticlesForParticle3D(temp, spatialLookup, startIndices, particles3D, totalNumberOfParticles, smoothingRadius);
	//std::cout << nbhdParticles.size() << std::endl;
	//for (int j = 0; j < nbhdParticles.size(); j++) {
	////	std::cout << "scalar field inside loop:" << scalarField << std::endl;
	//	int indexOfnbhdParticle = nbhdParticles[j];
	//	float distance = DistanceBetweenPoints(temp.predictedPosition,particles3D[indexOfnbhdParticle].predictedPosition);
	//	scalarField += SmoothingKernel(distance, particles3D[indexOfnbhdParticle].Gi) / particles3D[indexOfnbhdParticle].density;
	//	//std::cout << SmoothingKernel(distance, particles3D[indexOfnbhdParticle].Gi) / particles3D[indexOfnbhdParticle].density << std::endl;
	//}


	//scalarField = ForeachPointWithinRadiusCalculateDensity3D(temp, spatialLookup, startIndices, particles3D, totalNumberOfParticles, smoothingRadius);
	////std::cout << "scalar field:" << scalarField << std::endl;
	////exit(0);
	//float distance = DistanceBetweenPoints(particles3D[0].predictedPosition, temp.predictedPosition);
	//float influence = SmoothingKernel(smoothingRadius,distance);
	//float influence = SmoothingKernel(distance, particles3D[0].Gi);// std::cout << particles[400].Gi2D << std::endl;
	//scalarField = mass * influence;
	//scalarField = 0;
	//for (int i = 0; i < nbhdParticles.size();i++) {
	//	int particleIndex = nbhdParticles[i];
	//	float distance = DistanceBetweenPoints(particles3D[particleIndex].predictedPosition, temp.predictedPosition);
	//	//float influence = SmoothingKernel(smoothingRadius,distance);
	//	float influence = SmoothingKernel(distance, particles3D[particleIndex].Gi);// std::cout << particles[400].Gi2D << std::endl;
	//	if (influence < 0) {
	//		std::cout << "scalar field values are negative ffs";
	//		exit(0);
	//	}
	//	scalarField += mass * influence;
	//}
	//
	
	//scalarField = 0;
	//particle tempIndex;
	//tempIndex.p.x = particles3D[index3D].p.x; tempIndex.p.y = particles3D[index3D].p.y; tempIndex.p.z = particles3D[index3D].p.z;
	//tempIndex.predictedPosition.x = particles3D[index3D].predictedPosition.x; 
	//tempIndex.predictedPosition.y = particles3D[index3D].predictedPosition.y; 
	//tempIndex.predictedPosition.z = particles3D[index3D].predictedPosition.z;
	////getting the cell where the samplePoint is (this will be the centre of 3X3 block)
	//Vector3D centre = PositionToCellCoord(temp, smoothingRadius);
	//float sqrRadius = smoothingRadius * smoothingRadius;
	////if (samplePoint.particleIndex == 6488) std::cout << "before the loop" << std::endl;
	////Loop over all the cells in the  3X3X3 block around the centre cell
	//for (int i = (int)centre.x - 1; i <= (int)centre.x + 1; i++) {
	//	//if (samplePoint.particleIndex == 6488) std::cout << "entering the first loop" << std::endl;
	//	for (int j = (int)centre.y - 1; j <= (int)centre.y + 1; j++) {
	//		for (int k = (int)centre.z - 1; k <= (int)centre.z + 1; k++) {
	//			//std::cout << "we are inside the loop" << std::endl;
	//			//if (samplePoint.particleIndex == 6488) std::cout << "loop enetered dawg" << std::endl;
	//			//get key from the current cell, then loop over all the points that share that key
	//			Vector3D cellCoord; cellCoord.x = (float)i; cellCoord.y = (float)j; cellCoord.z = (float)k;
	//			uint64_t hash = HashCell(cellCoord);
	//			uint64_t key = GetKeyFromHash(hash, totalNumberOfParticles);
	{

		//			int cellStartIndex = startIndices[key];

		//			for (int iter = cellStartIndex; iter < totalNumberOfParticles; iter++) {
		//				//exit the loop if we're no longer looking at the correct cell
		//				if (spatialLookup[iter].cellkey != key) break;

		//				int particleIndex = spatialLookup[iter].particleIndex;
		//				float sqrDistance = SquareMagnitude(temp.predictedPosition, particles3D[particleIndex].predictedPosition);

		//				//test if the point is inside the radius
		//				if (sqrDistance <= sqrRadius) {
		//					float distance = DistanceBetweenPoints(temp.predictedPosition, particles3D[particleIndex].predictedPosition);
		//					Eigen::Matrix<float, 3, 1> r;
		//					r << temp.predictedPosition.x - particles3D[particleIndex].predictedPosition.x,
		//						temp.predictedPosition.y - particles3D[particleIndex].predictedPosition.y,
		//						temp.predictedPosition.z - particles3D[particleIndex].predictedPosition.z;
		//					float influence = SmoothingKernel(r, particles3D[particleIndex].Gi);
		//					//float influence = SmoothingKernel(distance, particles[particleIndex].Gi);
		//					scalarField += mass * influence;
		//				}
		//			}
		//		}


		//	}
		//}
	}

	temp.p.x = position.x; temp.p.y = position.y; temp.p.z = position.z;
	std::vector<int> nbhdParticles = GetNbhdPartilcesBruteForce(temp,particles3D,smoothingRadius);
	//std::cout << nbhdParticles.size()<<"mnbhd size \n";
	for (int i = 0; i < nbhdParticles.size();i++) {
		int particleIndex = nbhdParticles[i];
		//float distance = DistanceBetweenPoints(temp.predictedPosition, particles3D[particleIndex].predictedPosition);
		Eigen::Matrix<float, 3, 1> r;
		r << (temp.predictedPosition.x - particles3D[particleIndex].predictedPosition.x),
			 (temp.predictedPosition.y - particles3D[particleIndex].predictedPosition.y),
			 (temp.predictedPosition.z - particles3D[particleIndex].predictedPosition.z);
		float influence = SmoothingKernel(r, particles3D[particleIndex].Gi);
		if (influence == 0) {
			std::cout << (particles3D[particleIndex].Gi.determinant()) << std::endl;
		}
		//float influence = SmoothingKernel(distance, particles[particleIndex].Gi);
		scalarField += mass * influence; //std::cout << "particle sf:"<<scalarField << std::endl;
		//if (scalarField < 0) exit(420);
	}
	//std::cout << "Full scalar field: " << scalarField << std::endl;


	return scalarField;


}


