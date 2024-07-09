#pragma once
#include "FluidParticlesComputations.h"

static std::vector<int> GetNbhdPartilcesBruteForce(particle &samplePoint, std::vector<particle> &particles, float &h) {
    std::vector<int> nbhdParticles;
    for (int i = 0; i < particles.size();i++) {
        if (DistanceBetweenPoints(samplePoint.p, particles[i].p) <= h 
            /*&& i!=samplePoint.particleIndex*/) {
            nbhdParticles.push_back(i);
        }
    }

    return nbhdParticles;

}

static std::vector<int> GetNbhdParticlesForParticle3D(particle& samplePoint, std::vector<particleIndexCellKeyPair>& spatialLookup,
    std::vector<int> startIndices, std::vector<particle>& particles, int& totalNumberOfParticles, float& smoothingRadius) {

    std::vector<int> nbhdParticleIndexes;

    //getting the cell where the samplePoint is (this will be the centre of 3X3 block)
    Vector3D centre = PositionToCellCoord(samplePoint, smoothingRadius);
    float sqrRadius = smoothingRadius * smoothingRadius;

    //Loop over all the cells in the  3X3 block around the centre cell
    for (float i = centre.x - smoothingRadius; i <= centre.x + smoothingRadius; i+=smoothingRadius) {
        for (float j = centre.y - smoothingRadius; j <= centre.y + smoothingRadius; j+=smoothingRadius) {
            for (float k = centre.y - smoothingRadius; k <= centre.y + smoothingRadius; k+=smoothingRadius) {
                //std::cout << "we are in the cubes \n";
                //get key from the current cell, then loop over all the points that share that key
                Vector3D cellCoord; cellCoord.x = (float)i; cellCoord.y = (float)j; cellCoord.z = (float)k;
                uint64_t hash = HashCell(cellCoord);
                uint64_t key = GetKeyFromHash(hash, totalNumberOfParticles);
                
                         //   std::cout << key<<" this is the key\n";
                int cellStartIndex = startIndices[key];
                            //std::cout << cellStartIndex<<" i knew it\n";

                for (int iter = cellStartIndex; iter < totalNumberOfParticles; iter++) {
                            //std::cout << "entered key hehe\n";
                    //exit the loop if we're no longer looking at the correct cell
                    if (spatialLookup[iter].cellkey != key) break;

                    int particleIndex = spatialLookup[iter].particleIndex;
                    float sqrDistance = SquareMagnitude(samplePoint.predictedPosition, particles[particleIndex].predictedPosition);

                    //test if the point is inside the radius
                    if (sqrDistance <= sqrRadius) {
                        auto it = std::find(nbhdParticleIndexes.begin(), nbhdParticleIndexes.end(), particleIndex);

                        if (it != nbhdParticleIndexes.end()) {
                            //std::cout << "came to useless area\n";
                            // Element found
                            //std::cout << "particle already there" << std::endl;
                        }
                        else {
                            //std::cout << "letsgooo??\n";
                            // Element not found
                            //std::cout << "particle not there" << std::endl;
                            if (particleIndex != samplePoint.particleIndex) nbhdParticleIndexes.push_back(particleIndex);
                        }

                    }


                }



            }
        }

    }
    
    for (float i = centre.x - smoothingRadius; i <= centre.x + smoothingRadius; i += smoothingRadius) {
        for (float j = centre.y - smoothingRadius; j <= centre.y + smoothingRadius; j += smoothingRadius) {
            for (float k = centre.y - smoothingRadius; k <= centre.y + smoothingRadius; k += smoothingRadius) {
                std::cout << centre.x<<" , " << centre.y << " , " << centre.z << " , " << " entered the drawing grid loop\n";
                particle a; a.p.x = i, a.p.y = j, a.p.z = k;
                DrawPoint(a.p,1,0,1,1,5);
            }
        }
    }

    return nbhdParticleIndexes;


}

static int GetNumNbhdParticles(particle& samplePoint, std::vector<particleIndexCellKeyPair>& spatialLookup,
    std::vector<int> startIndices, std::vector<particle>& particles, int& totalNumberOfParticles, float& smoothingRadius) {
    std::vector<int> nbhdParticles = GetNbhdParticlesForParticle3D(samplePoint,spatialLookup,
        startIndices,particles,totalNumberOfParticles,smoothingRadius);

    return (int)nbhdParticles.size();
}

static Vector3D CalculatePropertyGradient3D(Vector3D& p, std::vector<particleIndexCellKeyPair>& spatialLookup,
    std::vector<int> startIndices, std::vector<particle>& particles, int& totalNumberOfParticles, float& smoothingRadius) {
    Vector3D propertyGradient;
    particle temp; temp.predictedPosition.x = p.x; temp.predictedPosition.y = p.y; temp.predictedPosition.z = p.z;
    temp.p.x = p.x; temp.p.y = p.y; temp.p.z = p.z;
    std::vector<int> nbhdParticles = GetNbhdParticlesForParticle3D(temp, spatialLookup, startIndices, particles, totalNumberOfParticles, smoothingRadius);

    for (int i = 0; i < nbhdParticles.size(); i++) {
        int particleIndex = nbhdParticles[i];

        float distance = DistanceBetweenPoints(p, particles[particleIndex].p);
        Vector3D direction;
        //TODO: we might be dividing by 0 here, so do SOMETHING! to solve it and prevent crashing
        direction.x = (p.x - particles[particleIndex].p.x) / distance;
        direction.y = (p.y - particles[particleIndex].p.y) / distance;
        direction.z = (p.z - particles[particleIndex].p.z) / distance;
        float slope = SmoothingKernelDerivative(smoothingRadius, distance);

        propertyGradient.x += particles[particleIndex].density * direction.x * slope * mass / particles[particleIndex].density;
        propertyGradient.y += particles[particleIndex].density * direction.y * slope * mass / particles[particleIndex].density;
        propertyGradient.z += particles[particleIndex].density * direction.z * slope * mass / particles[particleIndex].density;


    }
    return propertyGradient;
}



static float ForeachPointWithinRadiusCalculateDensity3D(particle& samplePoint, std::vector<particleIndexCellKeyPair>& spatialLookup,
    std::vector<int> startIndices, std::vector<particle>& particles, int& totalNumberOfParticles, float& smoothingRadius) {

    float density = 0; //density of that sample point

    //getting the cell where the samplePoint is (this will be the centre of 3X3 block)
    Vector3D centre = PositionToCellCoord(samplePoint, smoothingRadius);
    float sqrRadius = smoothingRadius * smoothingRadius;
    //if (samplePoint.particleIndex == 6488) std::cout << "before the loop" << std::endl;
    //Loop over all the cells in the  3X3X3 block around the centre cell
    for (int i = (int)centre.x - 1; i <= (int)centre.x + 1; i++) {
        //if (samplePoint.particleIndex == 6488) std::cout << "entering the first loop" << std::endl;
        for (int j = (int)centre.y - 1; j <= (int)centre.y + 1; j++) {
            for (int k = (int)centre.z - 1; k <= (int)centre.z + 1; k++) {
                //std::cout << "we are inside the loop" << std::endl;
                //if (samplePoint.particleIndex == 6488) std::cout << "loop enetered dawg" << std::endl;
                //get key from the current cell, then loop over all the points that share that key
                Vector3D cellCoord; cellCoord.x = (float)i; cellCoord.y = (float)j; cellCoord.z = (float)k;
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
                        float distance = DistanceBetweenPoints(samplePoint.predictedPosition, particles[particleIndex].predictedPosition);
                        float influence = SmoothingKernel(smoothingRadius, distance);
                        //float influence = SmoothingKernel(distance, particles[particleIndex].Gi);
                        density += mass * influence;
                    }
                }
            }


        }
    }

    return density;
}





static Vector3D ForeachPointWithinRadiusCalculatePressureForce3D(particle& samplePoint, std::vector<particleIndexCellKeyPair>& spatialLookup,
    std::vector<int> startIndices, std::vector<particle>& particles, int& totalNumberOfParticles, float& smoothingRadius) {

    Vector3D pressureForce;

    //getting the cell where the samplePoint is (this will be the centre of 3X3 block)
    Vector3D centre = PositionToCellCoord(samplePoint, smoothingRadius);
    float sqrRadius = smoothingRadius * smoothingRadius;

    //Loop over all the cells in the  3X3X3 block around the centre cell
    for (int i = (int)centre.x - 1; i <= (int)centre.x + 1; i++) {
        for (int j = (int)centre.y - 1; j <= (int)centre.y + 1; j++) {
            for (int k = (int)centre.z - 1; k <= (int)centre.z + 1; k++) {
                //get key from the current cell, then loop over all the points that share that key
                Vector3D cellCoord; cellCoord.x = (float)i; cellCoord.y = (float)j; cellCoord.z = (float)k;
                uint64_t hash = HashCell(cellCoord);
                uint64_t key = GetKeyFromHash(hash, totalNumberOfParticles);

                int cellStartIndex = startIndices[key];

                for (int iter = cellStartIndex; iter < totalNumberOfParticles; iter++) {
                    //exit the loop if we're no longer looking at the correct cell
                    if (spatialLookup[iter].cellkey != key) break;

                    if (samplePoint.particleIndex == spatialLookup[iter].particleIndex) continue;

                    int particleIndex = spatialLookup[iter].particleIndex;
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
                            direction.z = (float)((double)std::rand() / (double)RAND_MAX - 0.5) * 2;
                            //normalization, making it into a unit vector
                            direction.x = direction.x / Magnitude(direction);
                            direction.y = direction.y / Magnitude(direction);
                            direction.z = direction.z / Magnitude(direction);
                        }
                        else {
                            direction.x = (samplePoint.predictedPosition.x - particles[particleIndex].predictedPosition.x) / distance;
                            direction.y = (samplePoint.predictedPosition.y - particles[particleIndex].predictedPosition.y) / distance;
                            direction.z = (samplePoint.predictedPosition.z - particles[particleIndex].predictedPosition.z) / distance;
                        }

                        float slope = SmoothingKernelDerivative(smoothingRadius, distance);
                        float sharedPressure = CalculateSharedPressure(particles[particleIndex].density, samplePoint.density);
                        //std::cout << particles[particleIndex].density << std::endl;
                        pressureForce.x += sharedPressure * /*ConvertDensityToPressure(particles[i][j].density) **/ direction.x *
                            (slope)*mass / (particles[particleIndex].density);
                        pressureForce.y += sharedPressure */* ConvertDensityToPressure(particles[i][j].density) **/ direction.y *
                            (slope)*mass / (particles[particleIndex].density);
                        pressureForce.z += sharedPressure */* ConvertDensityToPressure(particles[i][j].density) **/ direction.z *
                            (slope)*mass / (particles[particleIndex].density);

                    }

                }

            }
        }

    }
    return pressureForce;

}

static void UpdateParticle3D(particle& par, std::vector<particleIndexCellKeyPair> spacialLookup, std::vector<int> startIndices,
    std::vector<particle>& particles, int& totalNumberOfParticles, float& smoothingRadius) {

    //getting the pressure force
    par.pressureForce = ForeachPointWithinRadiusCalculatePressureForce3D(par, spacialLookup, startIndices, particles, totalNumberOfParticles, smoothingRadius);
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

    ResolveCollisions3D(par);
}