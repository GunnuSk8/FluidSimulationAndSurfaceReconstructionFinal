#pragma once
#include "MathHelpingFunc.h"
#include<Eigen/Dense>

constexpr static float SmoothingKernel(float& radius, float& distance) {
    if (distance >= radius) return 0;
    float volume = 3.14f * Power(radius, 4) / 6;

    return (radius - distance) * (radius - distance) / volume;
}

constexpr static float SmoothingKernelDerivative(float& radius, float& distance) {
    if (distance >= radius) return 0;

    float scale = -12 / (3.14f * Power(radius, 4));
    return scale * (radius - distance);
    //this is a negative value
    //means
    //as the distance increases, the property value decreases
}
inline float sigma = 1;
static float P(float Gr) {
    if (Gr > smoothingRadius) return 0;
    return 315 * pow((smoothingRadius * smoothingRadius - Gr * Gr), 3) / (64 * 3.14f * pow(smoothingRadius, 9));
}
static float G(float Gr, float r) {
    if (r > smoothingRadius) return 0; // this thing needs to be looked upon my friend
    float volume = 3.14f * Power(smoothingRadius, 4) / 6;

    return (smoothingRadius - Gr) * (smoothingRadius - Gr) / volume;
}
static float SmoothingKernel(Eigen::Matrix<float, 3, 1>& r, Eigen::Matrix<float, 3, 3>& Gi) {
    Eigen::Matrix<float, 3, 1> Gr = Gi * r;
  
    float smooth = sigma * Gi.determinant() * G(Gr.norm(), r.norm());


    return smooth;

}
static float SmoothingKernel(float& r, Eigen::Matrix<float, 3, 3>& Gi) {
    if (r > smoothingRadius) return 0;
    float smooth;
    Eigen::Matrix<float, 3, 3> Gr = Gi * r;
    /*   if (G.determinant()<0) exit(0);*/
    smooth = sigma * Gi.determinant() * G(Gr.norm(), r);
    //std::cout << "\nGr:"<< Gr;

    return smooth;

}
static float SmoothingKernel(float& r, Eigen::Matrix<float, 2, 2>& Gi) {
    if (r > smoothingRadius) return 0;
    float smooth;
    Eigen::Matrix<float, 2, 2> Gr = Gi * r;
    /*   if (G.determinant()<0) exit(0);*/
    smooth = sigma * Gi.determinant() * P(Gr.norm());
    //std::cout << "\nGr:"<< Gr;

    return smooth;

}
