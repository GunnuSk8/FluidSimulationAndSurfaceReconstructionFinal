#pragma once
#include "FluidParticlesComputations.h"

static std::vector<int> GetNbhdParticlesForParticle(particle& samplePoint, std::vector<particleIndexCellKeyPair>& spatialLookup,
    std::vector<int> startIndices, std::vector<particle>& particles, int& totalNumberOfParticles, float& smoothingRadius) {

    std::vector<int> nbhdParticleIndexes;

    //getting the cell where the samplePoint is (this will be the centre of 3X3 block)
    Vector3D centre = PositionToCellCoord(samplePoint, smoothingRadius);
    float sqrRadius = smoothingRadius * smoothingRadius;

    //Loop over all the cells in the  3X3 block around the centre cell
    for (int i = (int)centre.x - 1; i <= (int)centre.x + 1; i++) {
        for (int j = (int)centre.y - 1; j <= (int)centre.y + 1; j++) {

            //get key from the current cell, then loop over all the points that share that key
            Vector3D cellCoord; cellCoord.x = (float)i; cellCoord.y = (float)j;
            uint64_t hash = HashCell(cellCoord);
            uint64_t key = GetKeyFromHash(hash, totalNumberOfParticles);

            int cellStartIndex = startIndices[key];

            for (int iter = cellStartIndex; iter < totalNumberOfParticles; iter++) {
                //exit the loop if we're no longer looking at the correct cell
                if (spatialLookup[iter].cellkey != key) break;

                int particleIndex = spatialLookup[iter].particleIndex;
                float sqrDistance = SquareMagnitude(samplePoint.predictedPosition, particles[particleIndex].predictedPosition);

                //test if the point is inside the radius
                if (sqrDistance <= sqrRadius) {
                    //do something with the particleIndex
                    //writing the code directly here
                    //getting the nbhd particles here
                    nbhdParticleIndexes.push_back(particleIndex);

                }

            }



        }

    }

    return nbhdParticleIndexes;


}

static float OptimisedCalculateDensityForParticle(particle& samplePoint, std::vector<particleIndexCellKeyPair>& spatialLookup,
    std::vector<int> startIndices, std::vector<particle>& particles, int& totalNumberOfParticles, float& smoothingRadius) {

    float density = 0; //density of that sample point

    std::vector<int> nbhdParticleIndexes = GetNbhdParticlesForParticle(samplePoint, spatialLookup, startIndices, particles, totalNumberOfParticles, smoothingRadius);

    for (int i = 0; i < nbhdParticleIndexes.size(); i++) {
        float distance = DistanceBetweenPoints(samplePoint.predictedPosition, particles[nbhdParticleIndexes[i]].predictedPosition);
        float influence = SmoothingKernel(smoothingRadius, distance);
        density += mass * influence;
    }


    return density;
}

static Vector3D OptimisedCalculatePressureForceForParticle(particle& samplePoint, std::vector<particleIndexCellKeyPair>& spatialLookup,
    std::vector<int> startIndices, std::vector<particle>& particles, int& totalNumberOfParticles, float& smoothingRadius) {

    Vector3D pressureForce;

    //std::vector<int> nbhdParticleIndexes = GetNbhdParticlesForParticle(samplePoint,spatialLookup,startIndices,particles,totalNumberOfParticles,smoothingRadius);
    //for (int i = 0; i < nbhdParticleIndexes.size(); i++) {
    //    int particleIndex = nbhdParticleIndexes[i];
    //    float distance = DistanceBetweenPoints(samplePoint.predictedPosition, particles[particleIndex].predictedPosition);
    //    //if both the particles are at the same location, then assign a random direction
    //    Vector3D direction;
    //    if (distance < 0.0001) {
    //        //std::cout << distance << "crash alert" << std::endl;
    //        //assigning a random direction if the distance is zero
    //        //this random direction is not working very well apparently
    //        direction.x = ((double)std::rand() / (double)RAND_MAX - 0.5) * 2;
    //        direction.y = ((double)std::rand() / (double)RAND_MAX - 0.5) * 2;
    //        //normalization, making it into a unit vector
    //        direction.x = direction.x / Magnitude(direction);
    //        direction.y = direction.y / Magnitude(direction);
    //    }
    //    else {
    //        direction.x = (samplePoint.predictedPosition.x - particles[particleIndex].predictedPosition.x) / distance;
    //        direction.y = (samplePoint.predictedPosition.y - particles[particleIndex].predictedPosition.y) / distance;
    //    }
    //    float slope = SmoothingKernelDerivative(smoothingRadius, distance);
    //    float sharedPressure = CalculateSharedPressure(particles[particleIndex].density, samplePoint.density);
    //    pressureForce.x += sharedPressure * /*ConvertDensityToPressure(particles[i][j].density) **/ direction.x *
    //        (slope)*mass / (particles[particleIndex].density);
    //    pressureForce.y += sharedPressure */* ConvertDensityToPressure(particles[i][j].density) **/ direction.y *
    //        (slope)*mass / (particles[particleIndex].density);
    //}
    
    //getting the cell where the samplePoint is (this will be the centre of 3X3 block)
    Vector3D centre = PositionToCellCoord(samplePoint, smoothingRadius);
    float sqrRadius = smoothingRadius * smoothingRadius;
    //Loop over all the cells in the  3X3 block around the centre cell
    for (int i = (int)centre.x - 1; i <= (int)centre.x + 1; i++) {
        for (int j = (int)centre.y - 1; j <= (int)centre.y + 1; j++) {
            //get key from the current cell, then loop over all the points that share that key
            Vector3D cellCoord; cellCoord.x = (float)i; cellCoord.y = (float)j;
            uint64_t hash = HashCell(cellCoord);
            uint64_t key = GetKeyFromHash(hash, totalNumberOfParticles);
            int cellStartIndex = startIndices[key];
            for (int k = cellStartIndex; k < totalNumberOfParticles; k++) {
                //exit the loop if we're no longer looking at the correct cell
                if (spatialLookup[k].cellkey != key) break;
                if (samplePoint.particleIndex == spatialLookup[k].particleIndex) continue;
                int particleIndex = spatialLookup[k].particleIndex;
                float sqrDistance = SquareMagnitude(samplePoint.predictedPosition, particles[particleIndex].predictedPosition);
                //test if the point is inside the radius
                if (sqrDistance <= sqrRadius) {
                    float distance = DistanceBetweenPoints(samplePoint.predictedPosition, particles[particleIndex].predictedPosition);
                    //if both the particles are at the same location, then assign a random direction
                    Vector3D direction;
                    if (distance < 0.0001) {
                        //std::cout << distance << "crash alert" << std::endl;
                        //assigning a random direction if the distance is zero
                        //this random direction is not working very well apparently
                        direction.x = (float)((double)std::rand() / (double)RAND_MAX - 0.5) * 2;
                        direction.y = (float)((double)std::rand() / (double)RAND_MAX - 0.5) * 2;
                        //normalization, making it into a unit vector
                        direction.x = direction.x / Magnitude(direction);
                        direction.y = direction.y / Magnitude(direction);
                    }
                    else {
                        direction.x = (samplePoint.predictedPosition.x - particles[particleIndex].predictedPosition.x) / distance;
                        direction.y = (samplePoint.predictedPosition.y - particles[particleIndex].predictedPosition.y) / distance;
                    }
                    float slope = SmoothingKernelDerivative(smoothingRadius, distance);
                    float sharedPressure = CalculateSharedPressure(particles[particleIndex].density, samplePoint.density);
                    pressureForce.x += sharedPressure * /*ConvertDensityToPressure(particles[i][j].density) **/ direction.x *
                        (slope)*mass / (particles[particleIndex].density);
                    pressureForce.y += sharedPressure */* ConvertDensityToPressure(particles[i][j].density) **/ direction.y *
                        (slope)*mass / (particles[particleIndex].density);
                }
            }
        }
    }
    return pressureForce;

}

static void UpdateParticle(particle& par, std::vector<particleIndexCellKeyPair> spacialLookup, std::vector<int> startIndices,
    std::vector<particle>& particles, int& totalNumberOfParticles, float& smoothingRadius) {

    //getting the pressure force
    par.pressureForce = OptimisedCalculatePressureForceForParticle(par, spacialLookup, startIndices, particles, totalNumberOfParticles, smoothingRadius);
    //par.pressureForce = CalculatePressureForce(par.particleIndex, par.p, particles, totalNumberOfParticles, smoothingRadius);
    //getting the acceleration
    par.a.x = par.pressureForce.x / par.density;
    par.a.y = par.pressureForce.y / par.density;
    par.a.z = par.pressureForce.z / par.density;
    //updating the particle velocities
    par.v.x += deltaTime * par.a.x;
    par.v.y += deltaTime * par.a.y;
    par.v.z += deltaTime * par.a.z;
    //external forces
    par.v.y += deltaTime * g;
    //updating positions
    par.p.x += deltaTime * par.v.x;
    par.p.y += deltaTime * par.v.y;
    par.p.z += deltaTime * par.v.z;

    ResolveCollisions(par);
}