#pragma once
#include "FluidParticlesComputations2D.h"
#include "FluidParticlesComputations3D.h"


static void UpdateDensities(std::vector<particleIndexCellKeyPair>& spacialLookup, std::vector<int>& startIndices,
    std::vector<particle>& particles, int& totalNumberOfParticles, float& smoothingRadius) {

    std::for_each(std::execution::par, particles.begin(), particles.end(),
        [&spacialLookup, &startIndices, &particles, &totalNumberOfParticles, &smoothingRadius](particle& par) {
            par.density = OptimisedCalculateDensityForParticle(par, spacialLookup, startIndices, particles, totalNumberOfParticles, smoothingRadius);
        });

}
static void UpdateDensities3D(std::vector<particleIndexCellKeyPair>& spacialLookup, std::vector<int>& startIndices,
    std::vector<particle>& particles, int& totalNumberOfParticles, float& smoothingRadius) {


    std::for_each(std::execution::par, particles.begin(), particles.end(),
        [&spacialLookup, &startIndices, &particles, &totalNumberOfParticles, &smoothingRadius](particle& par) {
           
            par.density = ForeachPointWithinRadiusCalculateDensity3D(par, spacialLookup, startIndices, particles, totalNumberOfParticles, smoothingRadius);
            //par.density = CalculateDensity(par.predictedPosition,particles,totalNumberOfParticles,smoothingRadius);
            //par.density = CalculateDensity(par.p, particles, totalNumberOfParticles, smoothingRadius);
            //std::cout << par.density<<std::endl;

        });

}





static void UpdateParticles(std::vector<particleIndexCellKeyPair> spacialLookup, std::vector<int> startIndices,
    std::vector<particle>& particles, int& totalNumberOfParticles, float& smoothingRadius) {

    std::for_each(std::execution::par, particles.begin(), particles.end(),
        [&spacialLookup, &startIndices, &particles, &totalNumberOfParticles, &smoothingRadius](particle& par) {
            UpdateParticle(par, spacialLookup, startIndices, particles, totalNumberOfParticles, smoothingRadius);
        });

}
static void UpdateParticles3D(std::vector<particleIndexCellKeyPair> spacialLookup, std::vector<int> startIndices,
    std::vector<particle>& particles, int& totalNumberOfParticles, float& smoothingRadius) {

    std::for_each(std::execution::par, particles.begin(), particles.end(),
        [&spacialLookup, &startIndices, &particles, &totalNumberOfParticles, &smoothingRadius](particle& par) {
            UpdateParticle3D(par, spacialLookup, startIndices, particles, totalNumberOfParticles, smoothingRadius);
        });

}


static void UpdateParticlesPredictedPosition(std::vector<particle>& particles) {
    std::for_each(std::execution::par, particles.begin(), particles.end(),
        [](particle& par) {
            //par.v.y += deltaTime * g;

            par.predictedPosition.x = par.p.x + par.v.x * deltaTime;
            par.predictedPosition.y = par.p.y + par.v.y * deltaTime;
            par.predictedPosition.z = par.p.z + par.v.z * deltaTime;
        });
}
