#pragma once
#include "SurfaceReconstructionComputations.h"
#include "SurfaceReconstructionComputations3D.h"
#include "SurfaceReconstructionComputations2D.h"


static void UpdateKernelCenters(std::vector<particle>& particles) {
	std::for_each(std::execution::par, particles.begin(), particles.end(),
		[&particles](particle& par) {
			float summationWij = GetSummationWij(par, particles);

			Vector3D summationWijTimesXj = GetSummationWijTimesXj(par, particles);

			par.updatedKernelCenter.x = (1.0f - lambda) * par.p.x + lambda * summationWijTimesXj.x / summationWij;
			par.updatedKernelCenter.y = (1.0f - lambda) * par.p.y + lambda * summationWijTimesXj.y / summationWij;
			par.updatedKernelCenter.z = (1.0f - lambda) * par.p.z + lambda * summationWijTimesXj.z / summationWij;
		});

}

static void CalculateCiForAllParticles(std::vector<particle>& particles3D) {
	std::for_each(std::execution::par, particles3D.begin(), particles3D.end(),
		[&particles3D](particle& par) {
			CalculateCi(par, particles3D);
		});
}
static void CalculateCiTildaForAllParticles(std::vector<particleIndexCellKeyPair>& spatialLookup,
	std::vector<int> startIndices, std::vector<particle>& particles3D, int& totalNumberOfParticles, float& smoothingRadius) {

	std::for_each(std::execution::par, particles3D.begin(), particles3D.end(),
		[&particles3D, &spatialLookup, &startIndices, &totalNumberOfParticles, &smoothingRadius](particle& par) {
			CalculateCiTilda(par, spatialLookup, startIndices, particles3D, totalNumberOfParticles, smoothingRadius);
		});
}
static void CalculateGiForAllParticles(std::vector<particle>& particles3D) {
	std::for_each(std::execution::par, particles3D.begin(), particles3D.end(),
		[](particle& par) {
			CalculateGi(par);
		});
}

static void CalculateCi2DForAllParticles(std::vector<particle>& particles3D) {
	std::for_each(std::execution::par, particles3D.begin(), particles3D.end(),
		[&particles3D](particle& par) {
			CalculateCi2D(par, particles3D);
		});
}
static void CalculateCi2DTildaForAllParticles(std::vector<particleIndexCellKeyPair>& spatialLookup,
	std::vector<int> startIndices, std::vector<particle>& particles3D, int& totalNumberOfParticles, float& smoothingRadius) {

	std::for_each(std::execution::par, particles3D.begin(), particles3D.end(),
		[&particles3D, &spatialLookup, &startIndices, &totalNumberOfParticles, &smoothingRadius](particle& par) {
			CalculateCi2DTilda(par, spatialLookup, startIndices, particles3D, totalNumberOfParticles, smoothingRadius);
		});
}
static void CalculateGi2DForAllParticles(std::vector<particle>& particles3D) {
	std::for_each(std::execution::par, particles3D.begin(), particles3D.end(),
		[](particle& par) {
			CalculateGi2D(par);
		});
}
