#pragma once

inline int totalNumberOfParticles = 421;
//for 2D water simulation and surface
std::vector<particle> particles(totalNumberOfParticles);
std::vector<particleIndexCellKeyPair> spacialLookup(totalNumberOfParticles);
std::vector<int> startIndices(totalNumberOfParticles);
//for 3D water simulation and surface
std::vector<particle> particles3D(totalNumberOfParticles);
std::vector<particleIndexCellKeyPair> spacialLookup3D(totalNumberOfParticles);
std::vector<int> startIndices3D(totalNumberOfParticles);


static void DisplayFunc3D() {

    DrawBoundingBox3D();

    UpdateKernelCenters(particles3D);

    UpdateParticlesPredictedPosition(particles3D);
    UpdateSpacialLookup(spacialLookup3D, startIndices3D, particles3D, totalNumberOfParticles, smoothingRadius);

    CalculateCiForAllParticles(particles3D);
    CalculateCiTildaForAllParticles(spacialLookup3D, startIndices3D, particles3D, totalNumberOfParticles, smoothingRadius);
    CalculateGiForAllParticles(particles3D);

    UpdateDensities3D(spacialLookup3D, startIndices3D, particles3D, totalNumberOfParticles, smoothingRadius);

    //UpdateParticles3D(spacialLookup3D, startIndices3D, particles3D, totalNumberOfParticles, smoothingRadius);
    DrawParticles(particles3D, totalNumberOfParticles);


    //DrawNbhdParticlesForParticle3D(particles3D[index3D], spacialLookup3D, startIndices3D, particles3D, totalNumberOfParticles, smoothingRadius);

    GridIntensityUpdate(spacialLookup3D, startIndices3D, particles3D, totalNumberOfParticles, smoothingRadius);
    MarchingCubes();
}
static void FluidSimulationAndSurface3D() {
    static float cRot = 20;
    glLoadIdentity();
    glPushMatrix();
    glTranslatef(0.0, 0.0, 0.0);
    glScalef(1, 1, 1);
    glRotatef(cRot, 1.0, 1.0, 0.0);

    DisplayFunc3D();

    cRot += 01;
    glPopMatrix();
    glFlush();

}
static void FluidSimulationAndSurface2D() {

    DrawBoundingBox();

    UpdateKernelCenters(particles);

    UpdateParticlesPredictedPosition(particles);
    UpdateSpacialLookup(spacialLookup, startIndices, particles, totalNumberOfParticles, smoothingRadius);
    UpdateDensities(spacialLookup, startIndices, particles, totalNumberOfParticles, smoothingRadius);

    CalculateCi2DForAllParticles(particles);
    CalculateCi2DTildaForAllParticles(spacialLookup, startIndices, particles, totalNumberOfParticles, smoothingRadius);
    CalculateGi2DForAllParticles(particles);


    //UpdateParticles(spacialLookup, startIndices, particles, totalNumberOfParticles, smoothingRadius);
    //std::cout << "velocity of indexed particle" << Magnitude(particles[index].v) << std::endl;
    DrawParticles(particles, totalNumberOfParticles);
    DrawCircle(particles[index].p, smoothingRadius);
    DrawNbhdParticlesForParticle(particles[index], spacialLookup, startIndices, particles, totalNumberOfParticles, smoothingRadius);

    Grid2DIntensityUpdate(spacialLookup, startIndices, particles, totalNumberOfParticles, smoothingRadius);
    MarchingSquares();

}
