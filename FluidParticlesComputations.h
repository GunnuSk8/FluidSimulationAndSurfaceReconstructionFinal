#pragma once
#include "DrawingHelpingFunc.h"
#include<iostream>
#include<vector>
#include<algorithm>
#include<execution>
#include "MathHelpingFunc.h"
#include<chrono>
#include "SmoothingFunctions.h"

static void ResolveCollisions(particle &pt) {
    if (abs(pt.p.x) > halfBoundSizeX - centerPtToBoundaryDist) {
        pt.p.x = (halfBoundSizeX - centerPtToBoundaryDist) * (pt.p.x > 0 ? 1 : -1);
        pt.v.x *= -1 * collisionDamping;
    }
    if (abs(pt.p.y) > halfBoundSizeY - centerPtToBoundaryDist) {
        pt.p.y = (halfBoundSizeY - centerPtToBoundaryDist) * (pt.p.y > 0 ? 1 : -1);
        pt.v.y *= -1 * collisionDamping;
    }
}
static void ResolveCollisions3D(particle &pt) {
    if (abs(pt.p.x) > halfBoundSizeX3D - centerPtToBoundaryDist) {
        pt.p.x = (halfBoundSizeX3D - centerPtToBoundaryDist) * (pt.p.x > 0 ? 1 : -1);
        pt.v.x *= -1 * collisionDamping;
    }
    if (abs(pt.p.y) > halfBoundSizeY3D - centerPtToBoundaryDist) {
        pt.p.y = (halfBoundSizeY3D - centerPtToBoundaryDist) * (pt.p.y > 0 ? 1 : -1);
        pt.v.y *= -1 * collisionDamping;
    } 
    if (abs(pt.p.z) > halfBoundSizeZ3D - centerPtToBoundaryDist) {
        pt.p.z = (halfBoundSizeZ3D - centerPtToBoundaryDist) * (pt.p.z > 0 ? 1 : -1);
        pt.v.z *= -1 * collisionDamping;
    }
}


//this function can run for both 3D and 2D
constexpr static float CalculateDensity(Vector3D &p, std::vector<particle> &particles, int &totalNumberOfParticles, float &smoothingRadius) {
    
    float density = 0;
    
    //we cant use parallel processing here because we will get wrong values
    for (int i = 0; i < totalNumberOfParticles; i++) {
            float distance = DistanceBetweenPoints(p, particles[i].predictedPosition);
            
            //float distance = DistanceBetweenPoints(p, particles[i].p); //this is for density visualisation under static conditions, 
            //because then the predicted positions are not updated every frame
            
            float influence = SmoothingKernel(smoothingRadius, distance);
            density += mass * influence;
    }
    
    return density;
}

//can be used for both 2D and 3D
static float CalculateProperty(Vector3D &p,std::vector<particle> &particles, int &totalNumberOfParticles, float &smoothingRadius) {
    float property = 0;
    
    for (int i = 0; i < totalNumberOfParticles; i++) {
       
            float distance = DistanceBetweenPoints(p, particles[i].p);
            float influence = SmoothingKernel(smoothingRadius, distance);
            property += particles[i].particleProperty * influence * mass / particles[i].density;

    }
    return property;
}

//can be used for 2d and 3d
static Vector3D CalculatePropertyGradient(Vector3D &p, std::vector<particle> &particles, int &totalNumberOfParticles, float &smoothingRadius) {
    Vector3D propertyGradient;
    

    for (int i = 0; i < totalNumberOfParticles; i++) {
        
            
            float distance = DistanceBetweenPoints(p, particles[i].p);
            Vector3D direction;
            //TODO: we might be dividing by 0 here, so do SOMETHING! to solve it and prevent crashing
            direction.x = (p.x - particles[i].p.x )/distance;
            direction.y = (p.y - particles[i].p.y)/distance;
            float slope = SmoothingKernelDerivative(smoothingRadius, distance);
            
            propertyGradient.x += particles[i].particleProperty * direction.x * slope * mass / particles[i].density  ; 
            propertyGradient.y += particles[i].particleProperty * direction.y * slope * mass / particles[i].density  ; 

        
    }
    return propertyGradient;
}

static float ConvertDensityToPressure(float &density) {
    float densityError = targetDensity - density;
    float pressure = densityError * pressureMultiplier;
    return (pressure); //because the particles were being attracted to each other, have to check later, what is the exact issue
    //use abs(pressure) if things mess up, that will make the particles go into empty spaces only
}

static float CalculateSharedPressure(float &densityA, float &densityB) {
    float pressureA = ConvertDensityToPressure(densityA);
    float pressureB = ConvertDensityToPressure(densityB);
    return ( pressureA + pressureB ) / 2;
}


//can be used for 2d and 3d
static Vector3D CalculatePressureForce(int &particleIndex, Vector3D &p, std::vector<particle> &particles,
    int &totalNumberOfParticles, float &smoothingRadius) {
    Vector3D pressureForce;
   
    for (int i = 0; i < totalNumberOfParticles; i++) {
        
            if (i == particleIndex) continue;

            float distance = DistanceBetweenPoints(p, particles[i].predictedPosition);
         
          
            Vector3D direction;
            if (distance < 0.0001) {
                
                //assigning a random direction if the distance is zero
                //this random direction is not working very well apparently
                direction.x = (float)((double)std::rand() / (double)RAND_MAX - 0.5f) * 2;
                direction.y = (float)((double)std::rand() / (double)RAND_MAX - 0.5f) * 2;
                direction.z = (float)((double)std::rand() / (double)RAND_MAX - 0.5f) * 2;
                //normalization, making it into a unit vector
                direction.x = direction.x / Magnitude(direction);
                direction.y = direction.y / Magnitude(direction);
                direction.z = direction.z / Magnitude(direction);
            }
            else {
                direction.x = (p.x - particles[i].predictedPosition.x) / distance;
                direction.y = (p.y - particles[i].predictedPosition.y) / distance;
                direction.z = (p.z - particles[i].predictedPosition.z) / distance;
            }
            
            float slope = SmoothingKernelDerivative(smoothingRadius, distance);
            float sharedPressure = CalculateSharedPressure(particles[particleIndex].density, particles[i].density); 

            pressureForce.x += sharedPressure * /*ConvertDensityToPressure(particles[i].density) **/ direction.x * 
                (slope) * mass / (particles[i].density);
            pressureForce.y += sharedPressure */* ConvertDensityToPressure(particles[i].density) **/ direction.y * 
                (slope) * mass / (particles[i].density);
            pressureForce.z += sharedPressure */* ConvertDensityToPressure(particles[i].density) **/ direction.z * 
                (slope) * mass / (particles[i].density);

        
    }
    return pressureForce;
}

static struct particleIndexCellKeyPair {
    int particleIndex=0;
    uint32_t cellkey=0;
};

static Vector3D PositionToCellCoord(particle &point, float &smoothingRadius) {
    Vector3D cellCoord;
    cellCoord.x = (point.predictedPosition.x / smoothingRadius);
    cellCoord.y = (point.predictedPosition.y / smoothingRadius);
    cellCoord.z = (point.predictedPosition.z / smoothingRadius);
    return cellCoord;
}

static uint64_t HashCell(Vector3D &Coord) {
    uint64_t a = (uint64_t)abs(Coord.x) * 15823;
    uint64_t b = (uint64_t)abs(Coord.y) * 9737333;
    uint64_t c = (uint64_t)abs(Coord.z) * 15823;
    return a + b + c;

}

static uint64_t GetKeyFromHash(uint64_t &hash, int totalNumberOfParticles) {
    return hash % (uint64_t)totalNumberOfParticles;
}

static void UpdateSpacialLookup(std::vector<particleIndexCellKeyPair> &spacialLookup, std::vector<int> &startIndices,
    std::vector<particle> &particles,int &totalNumberOfParticles, float &smoothingRadius) {

    //creating unordered spacial lookup
    std::for_each(std::execution::par,particles.begin(),particles.end(),
        [&spacialLookup,&startIndices,&totalNumberOfParticles,&smoothingRadius](particle &par) {
            Vector3D cellCoord = PositionToCellCoord(par, smoothingRadius);
            uint64_t cellHash = HashCell(cellCoord);
            uint64_t cellKey = GetKeyFromHash(cellHash, totalNumberOfParticles);

            spacialLookup[par.particleIndex].particleIndex = par.particleIndex;
            spacialLookup[par.particleIndex].cellkey = cellKey;

            startIndices[par.particleIndex] = INT16_MAX;

        });
    //creating unordered spatial lookup
    //for (int i = 0; i < totalNumberOfParticles; i++) {
    //    Vector2D cellCoord = PositionToCellCoord(particles[i],smoothingRadius);
    //    uint32_t cellHash = HashCell(cellCoord);
    //    uint32_t cellKey = GetKeyFromHash(cellHash,totalNumberOfParticles);
    //    //particleIndexCellKeyPair temp;
    //    spacialLookup[i].particleIndex = i;
    //    spacialLookup[i].cellkey = cellKey;
    //    //spacialLookup.push_back(temp);
    //    startIndices[i] = (INT16_MAX); //reset start index
    //}
  
    //sort by cell key
    std::sort(spacialLookup.begin(), spacialLookup.end(), [](const auto& a, const auto& b) {
        return a.cellkey < b.cellkey;
        });
  
    std::for_each(std::execution::par, particles.begin(),particles.end(),
        [&spacialLookup, &startIndices](particle &par) {
            int i = par.particleIndex;
            uint32_t key = spacialLookup[i].cellkey;
            //uint32_t keyPrev = i == 0 ? UINT32_MAX : spacialLookup[i - 1].cellkey;
            uint32_t keyPrev = i == 0 ? UINT32_MAX : spacialLookup[static_cast<std::vector<particleIndexCellKeyPair, std::allocator<particleIndexCellKeyPair>>::size_type>(i) - 1].cellkey;
            if (key != keyPrev) {
                startIndices[key] = i;
            }
        });
    //calculate start indices of each unique cell key in the spacial lookup
    //for (int i = 0; i < totalNumberOfParticles; i++) {
    //    uint32_t key = spacialLookup[i].cellkey;
    //    uint32_t keyPrev = i == 0 ? UINT32_MAX : spacialLookup[i - 1].cellkey;
    //    if (key != keyPrev) {
    //        startIndices[key] = i;
    //    }
    //}
  
}









static float ExampleFunction(Vector3D &p) {
    float b = 20, d = 3.6f, c = 6, a = 0.3f;
    return cos(b * p.y - a + c * sin(d * p.x));
}