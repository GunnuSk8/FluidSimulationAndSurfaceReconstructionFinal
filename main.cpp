#pragma once
#include<GL/glew.h>
#include<GLFW/glfw3.h>
#include<vector>
#include "FluidParticleDrawingFunctions.h"
#include "VisualiseDrawingFunctions.h"
#include "SurfaceReconstructionComputations.h"
#include "MarRChingCBueS.h"
#include"MaRchingSQuares.h"
#include<fstream>
#include "happly.h"
#include "imgui.h"
#include "imgui_impl_glfw_gl3.h"
#include "FluidParticlesComputationsForAllParticles.h"
#include "SurfaceReconstructionComputationsForAllParticles.h"
#include "config.h"
#include "FluidSimulationTests3DAnd2D.h"
#include "FunctionAndIdeaTests.h"
#include "PolytopeSurfaceReconstructionTests.h"





int main(void)
{
    SetParameters();
 
    InitialiseParticles(particles, totalNumberOfParticles);
    InitialiseParticles3D(particles3D,totalNumberOfParticles);
    InitialisePolytopeParticles(polytopeParticles);
    InitialiseParticles3DGradientTest(particles3DGradientTest);

    initnbhdCheck();
    
    GridInit();
    Grid2DInit();

   

   
    /*code to write points to a .csv file*/    
    //std::fstream myFile;
    //myFile.open("originalPoints2D.csv",std::ios::out); //write mode
    //if (myFile.is_open()) {
    //    for (int i = 0; i < totalNumberOfParticles; i++) {
    //        myFile << particles[i].p.x <<"," << particles[i].p.y<<"," << particles[i].p.z<<std::endl;
    //    }
    //    myFile.close();
    //}
    //myFile.open("updatedKernelCentersPoints2D.csv",std::ios::out);
    //if (myFile.is_open()) {
    //    for (int i = 0; i < totalNumberOfParticles; i++) {
    //        myFile << particles[i].updatedKernelCenter.x << "," << particles[i].updatedKernelCenter.y << "," << particles[i].updatedKernelCenter.z << std::endl;
    //    }
    //    myFile.close();
    //}
    //myFile.open("originalPoints3D.csv", std::ios::out); //write mode
    //if (myFile.is_open()) {
    //    for (int i = 0; i < totalNumberOfParticles; i++) {
    //        myFile << particles3D[i].p.x << "," << particles3D[i].p.y << "," << particles3D[i].p.z << std::endl;
    //    }
    //    myFile.close();
    //}
    //myFile.open("updatedKernelCentersPoints3D.csv", std::ios::out);
    //if (myFile.is_open()) {
    //    for (int i = 0; i < totalNumberOfParticles; i++) {
    //        myFile << particles3D[i].updatedKernelCenter.x << "," << particles3D[i].updatedKernelCenter.y << "," << particles3D[i].updatedKernelCenter.z << std::endl;
    //    }
    //    myFile.close();
    //}

   
    GLFWwindow* window;

    /* Initialize the library */
    if (!glfwInit())
        return -1;

    float windowSize = 1;
    /* Create a windowed mode window and its OpenGL context */
    window = glfwCreateWindow((int)(900 * windowSize), (int)(900 * windowSize), "Simulation", NULL, NULL);


    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    glewInit(); 
    ImGui::CreateContext();
    ImGui_ImplGlfwGL3_Init(window, true);
    ImGui::StyleColorsDark();

    std::vector<float> srList = {0.02,0.025,0.03,0.035,0.040,0.045,0.05,0.055,0.06,0.065,0.07,0.075};
    int ind = 0;
    /* Loop until the user closes the window */
    while (!glfwWindowShouldClose(window))
    {
        ImGui_ImplGlfwGL3_NewFrame();
        {
            ImGui::Text("Hello, world!");
        }
        /* Render here */
        glClear(GL_COLOR_BUFFER_BIT);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        //smoothingRadius = srList[ind];


        //my functions

        //displayNbhdCheck();


     /*   FluidSimulationAndSurface2D();
        {
            ImGui::InputFloat("Smoothing Radius : ", &smoothingRadius);
            ImGui::InputFloat("Kr2D: ", &Kr2D);
            ImGui::InputFloat("Ks2D: ", &Ks2D);
            ImGui::InputFloat("Kn2D: ", &Kn2D);
            ImGui::InputInt("Ne2D: ", &Ne2D);
            ImGui::InputInt("index: ", &index);
            ImGui::InputFloat("isovalue2D: ", &isoValue2D);
        }*/
       /* FluidSimulationAndSurface3D();
        {
            ImGui::InputFloat("isovalue: ", &isovalue);
            ImGui::InputFloat("Kr: ", &Kr);
            ImGui::InputFloat("Ks: ", &Ks);
            ImGui::InputFloat("Kn: ", &Kn);
            ImGui::InputInt("Ne: ", &Ne);
            ImGui::InputInt("index: ", &index3D);
        }*/

        DisplayPolytope();
        {
            ImGui::InputFloat("isovalue: ", &isovalue);
            ImGui::SliderFloat("h: ", &smoothingRadius, -1.0, 1.0);
            ImGui::InputFloat("Kr: ", &Kr);
            ImGui::InputFloat("Ks: ", &Ks);
            ImGui::InputFloat("Kn: ", &Kn);
            ImGui::InputInt("Ne: ", &Ne);
            ImGui::InputInt("index: ", &index3D);
        }
        //MC();
        
        /*DrawEllipse();
        {
            ImGui::SliderFloat("a: ", &a, -1.0, 1.0);
            ImGui::SliderFloat("b: ", &b, -1.0, 1.0);
            ImGui::SliderFloat("c: ", &c, -1.0, 1.0);
            ImGui::SliderFloat("d: ", &d, -1.0, 1.0);
            ImGui::SliderFloat("xpos: ", &xpos, -1.0, 1.0);
            ImGui::SliderFloat("ypos: ", &ypos, -1.0, 1.0);
        }*/
  
        //DrawEllipse3D();
        /*{
            ImGui::SliderFloat("a3D: ", &a3D, -1.0, 1.0);
            ImGui::SliderFloat("b3D: ", &b3D, -1.0, 1.0);
            ImGui::SliderFloat("c3D: ", &c3D, -1.0, 1.0);
            ImGui::SliderFloat("d3D: ", &d3D, -1.0, 1.0);
            ImGui::SliderFloat("e3D: ", &e3D, -1.0, 1.0);
            ImGui::SliderFloat("f3D: ", &f3D, -1.0, 1.0);
            ImGui::SliderFloat("g3D: ", &g3D, -1.0, 1.0);
            ImGui::SliderFloat("h3D: ", &h3D, -1.0, 1.0);
            ImGui::SliderFloat("i3D: ", &i3D, -1.0, 1.0);
            ImGui::SliderFloat("xpos3D: ", &xpos3D, -1.0, 1.0);
            ImGui::SliderFloat("ypos3D: ", &ypos3D, -1.0, 1.0);
            ImGui::SliderFloat("zpos3D: ", &zpos3D, -1.0, 1.0);
        }*/
        //DrawEllipse2D(c,b,a,d,20);
        {
            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
        }
        ImGui::Render();
        ImGui_ImplGlfwGL3_RenderDrawData(ImGui::GetDrawData());
        /* Swap front and back buffers */
        glfwSwapBuffers(window);
  
        /* Poll for and process events */
        glfwPollEvents();
        ind++;

        //if (ind >= srList.size()) exit(0);
    }

    ImGui_ImplGlfwGL3_Shutdown();
    ImGui::DestroyContext();
    glfwTerminate();
    return 0;
}