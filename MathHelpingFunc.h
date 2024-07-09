#pragma once
#include "FluidParticlesProperties&Init.h"
#include<cmath>
#include<iostream>
constexpr static float Max(float a, float b) {
    if (a > b)
        return a;
    else
        return b;

}

constexpr static float Min(float a, float b) {
    if (a < b)
        return a;
    else
        return b;

}

static float Magnitude(Vector3D &a) {
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

constexpr static float SquareMagnitude(Vector3D &a, Vector3D &b) {
    return ((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}

static float DistanceBetweenPoints(Vector3D &a, Vector3D &b) {
    Vector3D c;
    c.x = a.x - b.x;
    c.y = a.y - b.y;
    c.z = a.z - b.z;
    return Magnitude(c);
}

constexpr static float Power(float a, int b) {
    float ans = 1;
    for (int i = 0; i < b; i++) {
        ans *= a;
    }
    return ans;
}
