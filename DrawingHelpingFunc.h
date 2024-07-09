#pragma once
#include<GLFW/glfw3.h>
#include "FluidParticlesProperties&Init.h"
#include "MathHelpingFunc.h"
#include<iostream>

static void DrawPoint(Vector3D position, float r=1, float g=1, float b=1, float a=1,float particleSize = 10 ) {
    glPointSize(particleSize);
    glEnable(GL_POINT_SMOOTH);
    glBegin(GL_POINTS);
    glColor4f(r, g, b, a);
    glVertex3f(position.x, position.y, position.z);
    glEnd();
}

static void DrawPoint(point position, float r=1, float g=1, float b=1, float a=1) {
    glPointSize(4);
    glEnable(GL_POINT_SMOOTH);
    glBegin(GL_POINTS);
    glColor4f(r, g, b, a);
    glVertex3f(position.x, position.y, position.z);
    glEnd();
}

static void DrawParticle(particle par) {
    float speed = Magnitude(par.v);
   
    glPointSize(particleSize);
    glEnable(GL_POINT_SMOOTH);
    glBegin(GL_POINTS);
    //glColor4f(46.0f  / 256.0f, 163.0f / 256.0f, 242.0f / 256.0f, 0.3);
    glColor4f(speed,speed,0.7f,0.7f);
    glVertex3f(par.p.x, par.p.y,par.p.z);
    glEnd();

}

static void DrawLine(Vector3D &p1, Vector3D &p2, float r=1, float g=1, float b=1, float a=1) {
    glColor4f(r, g, b,a);
    glBegin(GL_LINES);
    glVertex3f(p1.x, p1.y,p1.z);
    glVertex3f(p2.x, p2.y,p2.z);
    glEnd();
}

static void DrawLine(point &p1, point &p2, float r=1, float g=1, float b=1, float a=1) {
    glColor4f(r, g, b,a);
    glBegin(GL_LINES);
    glVertex3f(p1.x, p1.y,p1.z);
    glVertex3f(p2.x, p2.y,p2.z);
    glEnd();
}

static float halfBoundSizeX = 0.8f;
static float halfBoundSizeY = 0.8f;
static void DrawBoundingBox() {
    Vector3D a, b, c, d;
    a.x = halfBoundSizeX; a.y = halfBoundSizeY;
    b.x = -halfBoundSizeX; b.y = halfBoundSizeY;
    c.x = -halfBoundSizeX; c.y = -halfBoundSizeY;
    d.x = halfBoundSizeX; d.y = -halfBoundSizeY;
    DrawLine(a, b);
    DrawLine(b, c);
    DrawLine(c, d);
    DrawLine(d, a);
}

static float halfBoundSizeX3D = 0.5;
static float halfBoundSizeY3D = 0.5;
static float halfBoundSizeZ3D = 0.5;

static void DrawBoundingBox3D() {
    
        glLineWidth(1.0f);
        glColor4f(255.0 / 256.0, 192.0 / 256.0, 203.0 / 256.0, 0.8f);
        glBegin(GL_LINE_LOOP);
        glVertex3f( halfBoundSizeX3D,  halfBoundSizeY3D,  halfBoundSizeZ3D);
        glVertex3f(-halfBoundSizeX3D,  halfBoundSizeY3D,  halfBoundSizeZ3D);
        glVertex3f(-halfBoundSizeX3D, -halfBoundSizeY3D,  halfBoundSizeZ3D);
        glVertex3f( halfBoundSizeX3D, -halfBoundSizeY3D,  halfBoundSizeZ3D);
        glEnd(); glBegin(GL_LINE_LOOP);
        glVertex3f( halfBoundSizeX3D,  halfBoundSizeY3D, -halfBoundSizeZ3D);
        glVertex3f(-halfBoundSizeX3D,  halfBoundSizeY3D, -halfBoundSizeZ3D);
        glVertex3f(-halfBoundSizeX3D, -halfBoundSizeY3D, -halfBoundSizeZ3D);
        glVertex3f( halfBoundSizeX3D, -halfBoundSizeY3D, -halfBoundSizeZ3D);
        glEnd();
        glBegin(GL_LINES);
        glVertex3f( halfBoundSizeX3D,  halfBoundSizeY3D,  halfBoundSizeZ3D);
        glVertex3f( halfBoundSizeX3D,  halfBoundSizeY3D, -halfBoundSizeZ3D);
        glEnd();
        glBegin(GL_LINES);
        glVertex3f(-halfBoundSizeX3D,  halfBoundSizeY3D,  halfBoundSizeZ3D);
        glVertex3f(-halfBoundSizeX3D,  halfBoundSizeY3D, -halfBoundSizeZ3D);
        glEnd(); glBegin(GL_LINES);
        glVertex3f(-halfBoundSizeX3D, -halfBoundSizeY3D,  halfBoundSizeZ3D);
        glVertex3f(-halfBoundSizeX3D, -halfBoundSizeY3D, -halfBoundSizeZ3D);
        glEnd(); glBegin(GL_LINES);
        glVertex3f( halfBoundSizeX3D, -halfBoundSizeY3D, halfBoundSizeZ3D);
        glVertex3f( halfBoundSizeX3D, -halfBoundSizeY3D, -halfBoundSizeZ3D);
        glEnd();
}


static void DrawCircle(Vector3D p, float radius = 0.2f,int mode=0) {
    if (mode == 0)
    {
        glColor4f(12.0f / 256.0f, 113.0f / 256.0f, 195.0f / 256.0f, 0.50);
    }
    glBegin(GL_POLYGON);
    for (int i = 0; i < 50; i++) {
        float theta = 2.0f * 3.14f * (float)i / (float)50;
        float x = p.x + radius * cos(theta);
        float y = p.y + radius * sin(theta);
        glVertex2d(x, y);
    }
    glEnd();
}

static void DrawCircleBoundary(Vector3D p, float radius, float r, float g, float b, float a) {
    glColor4f(r,g,b, a);
    glBegin(GL_LINE_LOOP);
    for (int i = 0; i < 50; i++) {
        float theta = 2.0f * 3.14f * (float)i / (float)50;
        float x = p.x + radius * cos(theta);
        float y = p.y + radius * sin(theta);
        glVertex2d(x, y);
    }
    glEnd();
}

static void DrawArrow(Vector3D &p, Vector3D &gradient) {
    glColor4f(1.0f,0.0f,0.0f,0.6f);
    Vector3D otherPoint; otherPoint.x = p.x + gradient.x; otherPoint.y = p.y + gradient.y; otherPoint.z = p.z + gradient.z;
    float L = Magnitude(gradient);
    float d = 0.05f;
    float dx, dy, dz;
    if (d < L) {
        dx = p.x + (gradient.x * d / L);
        dy = p.y + (gradient.y * d / L);
        dz = p.z + (gradient.z * d / L);

    }
    else {
        dx = p.x + (gradient.x );
        dy = p.y + (gradient.y);
        dz = p.z + (gradient.z);
    }
        
    glLineWidth(5.0);
    glBegin(GL_LINES);
    glVertex3f(p.x,p.y,p.z);
    glVertex3f( dx , dy, dz);
    glEnd();
    glColor4f(1.0f,1.0f,0.0f,0.80f);
    glPointSize(10.0f);
    glEnable(GL_POINT_SMOOTH);
    glBegin(GL_POINTS);
    glVertex3f(p.x, p.y, p.z);
    glEnd();

}

static void DrawSphere(Vector3D p) {
    float phi, theta;
    DrawPoint(p,1,0,0,1);
    for (phi = -3.14f/2; phi < 3.14f/2 + 0.01; phi+=0.01f) {
        for (theta = 0; theta < 2*3.14f; theta+=0.01f) {
            Vector3D pointPosition;
            pointPosition.x = p.x + smoothingRadius * cos(phi) * cos(theta);
            pointPosition.y = p.y + smoothingRadius * cos(phi) * sin(theta);
            pointPosition.z = p.z + smoothingRadius * sin(phi);


            glPointSize(3);
            glEnable(GL_POINT_SMOOTH);
            glBegin(GL_POINTS);
            glColor4f(1.0f, 1.0f, 1.0f, 0.01f);
            glVertex3f(pointPosition.x, pointPosition.y, pointPosition.z);
            glEnd();
        }
    }
}

