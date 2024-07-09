#pragma once
#include "FluidParticlesComputations.h"

static void DrawSmoothingKernelFunction() {

    glColor3f(1.0f, 1.0f, 1.0f);
    float smoothnessFactor = 0.01f;
    float radius = 1;
    glBegin(GL_LINE_STRIP);
    for (float i = 0; i <= 1; i += smoothnessFactor) {
        particle temp;
        temp.p.x = i;
        
        temp.p.y = SmoothingKernel(radius, i)/ 2;
        glVertex2f(temp.p.x, temp.p.y);
    }
    glEnd();
    glBegin(GL_LINE_STRIP);
    for (float i = 0; i <= 1; i += smoothnessFactor) {
        particle temp;
        temp.p.x = i;
        temp.p.y = SmoothingKernel(radius, i) / 2;
        glVertex2f(-temp.p.x, temp.p.y);
    }
    glEnd();
    glColor3f(0.0, 1.0, 1.0);
    glBegin(GL_LINE_STRIP);
    glVertex2f(-1.0, 0.0);
    glVertex2f(1.0, 0.0);
    glEnd();glBegin(GL_LINE_STRIP);
    glVertex2f(0.0, 1.0);
    glVertex2f(0.0, -1.0);
    glEnd();
   
}

static void visualiseDensityField(std::vector<particleIndexCellKeyPair>& spatialLookup,
    std::vector<int> startIndices, std::vector<particle>& particles, int &totalNumberOfParticles, float &smoothingRadius) {
    for (float i = -0.8f; i < 0.801f; i += 0.05f) {
        for (float j = -0.8f; j < 0.801f; j += 0.05f) {
            
            Vector3D p;
            p.x = i; p.y = j;
            particle par; par.predictedPosition.x = i; par.predictedPosition.y = j;
            //float density = CalculateDensity(p, particles, totalNumberOfParticles,  smoothingRadius);
            float density = OptimisedCalculateDensityForParticle(par,spatialLookup,startIndices,particles,totalNumberOfParticles,smoothingRadius);
            //std::cout << density << std::endl;
            
           
          /*  float c = Min(density/100, 1);
            glColor4f(c,c,c, 0.5);
            DrawCircle(p,0.05,1);*/
            
            glPointSize(5.0f);
            glEnable(GL_POINT_SMOOTH);
            glBegin(GL_POINTS);
            glVertex2f(p.x, p.y);
            glEnd();
            if (density < targetDensity - targetDensityThreshhold) {
                glColor4f(0.1266f, 0.5330f, 0.6886f, 0.95f);
                glPointSize(3.0f);
                glBegin(GL_POINTS); 
                glVertex2f(p.x, p.y);
                glEnd();
                //DrawCircle(p, 0.01, 1);
            }
            else if (density > targetDensity + targetDensityThreshhold) {
                glColor4f(0.8301f, 0.2914f, 0.2075f, 0.95f);
                glPointSize(3.0f);
                glBegin(GL_POINTS);
                glVertex2f(p.x, p.y);
                glEnd();
                //DrawCircle(p, 0.01, 1);
            }
            else {
                glColor4f(1.0f, 1.0f, 1.0f, 0.95f);
                glPointSize(3.0f);
                glBegin(GL_POINTS);
                glVertex2f(p.x, p.y);
                glEnd();
                //std::cout << "SomePoints have target density";
                //DrawCircle(p, 0.01, 1);
            }
            

        }
    }
}

static void visualisePropertyField(std::vector<particle> &particles, int &totalNumberOfParticles, float &smoothingRadius) {
    //to visualise the example function properly,
    //decrease the smoothing radius,=> decrease the area that influences the property at a particular point
    //and increase the total number of particles
    for (float i = -0.8f; i < 0.801; i += 0.01f) {
        for (float j = -0.8f; j < 0.801; j += 0.01f) {

            Vector3D p;
            p.x = i; p.y = j;
            float property = CalculateProperty(p, particles, totalNumberOfParticles, smoothingRadius);
            //std::cout << property << std::endl;
            float c = Min((property*2) , 1);
            
            glColor4f(c, c, c, 0.1f);
            DrawCircle(p, 0.01f, 1);
         

        }
    }
}

static void visualisePropertyGradientField(std::vector<particle> &particles, int &totalNumberOfParticles, float &smoothingRadius) {
    for (float i = -0.8f; i < 0.801; i += 0.1f) {
        for (float j = -0.8f; j < 0.801; j += 0.1f) {
            Vector3D p;
            p.x = i; p.y = j;
            Vector3D gradient = CalculatePropertyGradient(p,particles,totalNumberOfParticles,smoothingRadius);
            DrawArrow(p,gradient);

        }
    }
}

static void DrawGradient3DAtAPoint(point p, std::vector<particleIndexCellKeyPair>& spatialLookup,
    std::vector<int> startIndices, std::vector<particle>& particles, int& totalNumberOfParticles, float& smoothingRadius) {
    Vector3D pos; pos.x = p.x; pos.y = p.y; pos.z = p.z;
    Vector3D grad = CalculatePropertyGradient3D(pos,spatialLookup,startIndices,particles,totalNumberOfParticles,smoothingRadius);
    DrawArrow(pos,grad);
}

static void visualiseExampleFunction() {

    for (float i = -0.8f; i < 0.801; i += 0.01f) {
        for (float j = -0.8f; j < 0.801; j += 0.01f) {
            Vector3D p;
            p.x = i; p.y = j;
            float exampleFunction = ExampleFunction(p);
            glColor4f(exampleFunction,exampleFunction,exampleFunction,0.9f);
            glPointSize(1.0f);
            glBegin(GL_POINTS);
            glVertex2f(p.x,p.y);
            glEnd();
            //DrawCircle(p, 0.1, 1);

        }
    }
}