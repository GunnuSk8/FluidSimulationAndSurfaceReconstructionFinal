#pragma once
#include<cmath>
#include<vector>
#include<Eigen/Dense>


static struct Vector3D {
    float x = 0;
    float y = 0;
    float z = 0;
};

//this for marching cubes 
struct point {
    float x=0;
    float y=0;
    float z=0;
    float intensity=0;
};

static struct particle {
    int particleIndex = 0 ;
    Vector3D p; //position
    Vector3D v; //velocity
    float particleProperty = 0;
    float density = 1;
    Vector3D a; //acceleration
    Vector3D pressureForce;
    Vector3D predictedPosition;
    Vector3D updatedKernelCenter;
    Eigen::Matrix<float, 3, 3> Ci = Eigen::Matrix<float, 3, 3>::Zero();
    Eigen::Matrix<float, 3, 3> CiTilda = Eigen::Matrix<float, 3, 3>::Zero();
    Eigen::Matrix<float, 3, 3> Gi = Eigen::Matrix<float, 3, 3>::Zero();
    Eigen::Matrix<float, 2, 2> Ci2D = Eigen::Matrix<float, 2, 2>::Zero();
    Eigen::Matrix<float, 2, 2> Ci2DTilda = Eigen::Matrix<float, 2, 2>::Zero();
    Eigen::Matrix<float, 2, 2> Gi2D = Eigen::Matrix<float, 2, 2>::Zero();
};

inline float particleSize = 10.0f;
inline float centerPtToBoundaryDist = 0.0125f;
inline float g = -10.0f; //gravity acceleration
inline float deltaTime = 0.005f;//constant for time step for calculatin acceleration, as used in the video
inline float collisionDamping = 0.500f;
inline float targetDensity = 10.2f;
inline float targetDensityThreshhold = 0.1f * targetDensity;
inline float pressureMultiplier = 10.0f;
inline float mass = 1.0f;
inline float smoothingRadius = 0.055f;
inline Vector3D centroidPolytope;
