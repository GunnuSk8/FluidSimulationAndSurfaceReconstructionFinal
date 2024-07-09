#pragma once
#include<string>
//for writing the .ply file of the surface that I extracted



inline int numPT = numberOfParticles;
//for surface of organ point cloud
std::vector<particle> polytopeParticles(numPT);
std::vector<particleIndexCellKeyPair> spacialLookupPT(numPT);
std::vector<int> startIndicesPT(numPT);
//for updated kernel approach to test if we get the surface or not, this needs to be cleared up dawg
std::vector<particle> polytopeParticlesUpdatedKernel(numberOfParticles);
std::vector<particleIndexCellKeyPair> spacialLookupPTUpdatedKernel(numberOfParticles);
std::vector<int> startIndicesPTUpdatedKernel(numberOfParticles);

//declaring stuff for the isolated points with 0 nbhd particles
std::vector<particle> isolatedPoints(1482);
std::vector<particleIndexCellKeyPair> spatialLookupIsolatedPoints(1482);
std::vector<int> startIndicesIsolatedPoints(1482);

inline int count = 0;
//write function to write a .ply file for the mesh that I have extracted

static std::string GetFileName() {
    std::string fileName = "CorrectedAsliAnisotropic_";
    fileName += "smoothingRadius"; fileName += std::to_string(smoothingRadius);
    fileName += "rez";             fileName += std::to_string(rez);
    fileName += "isoValue";        fileName += std::to_string(isovalue);
    fileName += "lambda";          fileName += std::to_string(lambda);
    fileName += "riScale";         fileName += std::to_string(riScale);
    fileName += "Kr";              fileName += std::to_string(Kr);
    fileName += "Ks";              fileName += std::to_string(Ks);
    fileName += "Kn";              fileName += std::to_string(Kn);
    fileName += "Ne";              fileName += std::to_string(Ne);
   

    // Iterate over each character in the string
    for (char& c : fileName) {
        if (c == '.') { // If the character is a '.'
            c = '_';    // Replace it with '_'
        }
    }

    std::cout << fileName << std::endl; // Output: This_is_a_test_string
    
    return fileName;
}


static void DisplayPolytope() {
    static float cRotx = 20;
    static float cRoty = 20;
    static float cRotz = 20;
    glLoadIdentity();
    glPushMatrix();
    glTranslatef(0.0, 0.0, 0.0);

    static float sc = 1;
    glScalef(sc, sc, sc);
    glRotatef(cRotx, cRoty, cRotz, 0.0);

    //DrawBoundingBox3D();

    UpdateKernelCenters(polytopeParticles);
    UpdateParticlesPredictedPosition(polytopeParticles);
    //UpdateSpacialLookup(spacialLookupPT, startIndicesPT, polytopeParticles, numPT,smoothingRadius);

    CalculateCiForAllParticles(polytopeParticles);
    CalculateCiTildaForAllParticles(spacialLookupPT, startIndicesPT, polytopeParticles, numPT, smoothingRadius);
    CalculateGiForAllParticles(polytopeParticles);

    //Drawing the polytope particles
    //std::vector<int> nbhdPar = GetNbhdParticlesForParticle3D(polytopeParticles[index3D], spacialLookupPT, startIndicesPT,
      //  polytopeParticles, numPT, smoothingRadius);
    //drawing the neighbourhood particles for a particular particle
    
    
    float minSV = INT16_MAX;
    float maxSV = 0;
    //printing all the particles
    for (int i = 0; i < polytopeParticles.size(); i++) {
        DrawPoint(polytopeParticles[i].p, 1, 0, 0, 0.9f, 3);
        point a; 
        a.x = polytopeParticles[i].p.x;
        a.y = polytopeParticles[i].p.y;
        a.z = polytopeParticles[i].p.z;
        float scalarField = ScalarFieldAtPosition(a,spacialLookupPT,startIndicesPT,polytopeParticles,numPT,smoothingRadius);
        if (scalarField < minSV) minSV = scalarField;
        if (scalarField > maxSV) maxSV = scalarField;
    }
       std::cout << "min scalar field at original particles: " << minSV << std::endl;
       std::cout << "max scalar field at original particles: " << maxSV << std::endl;

    //    std::vector<int> nbhdPar = GetNbhdPartilcesBruteForce(polytopeParticles[index3D],
    //        polytopeParticles, smoothingRadius);
    //    std::cout << nbhdPar.size() << std::endl;
    ////printing the index particle
    //DrawPoint(polytopeParticles[index3D].p, 1, 0, 0, 1, 5);
    //DrawSphere(polytopeParticles[index3D].p);
    ////drawing the nbhd points
    //for (int i = 0; i < nbhdPar.size(); i++) {
    //    DrawPoint(polytopeParticles[nbhdPar[i]].p, 0, 1, 1, 1, 3);
    //}

   
    //copying the updated kernel centre position to polytopparticleUpdatedKernel
    for (int i = 0; i < numberOfParticles; i++) {
        polytopeParticlesUpdatedKernel[i].p.x = polytopeParticles[i].updatedKernelCenter.x;
        polytopeParticlesUpdatedKernel[i].p.y = polytopeParticles[i].updatedKernelCenter.y;
        polytopeParticlesUpdatedKernel[i].p.z = polytopeParticles[i].updatedKernelCenter.z;
        polytopeParticlesUpdatedKernel[i].particleIndex = i;
    }
    //UpdateKernelCenters(polytopeParticlesUpdatedKernel);
    auto start = std::chrono::high_resolution_clock::now();
    UpdateParticlesPredictedPosition(polytopeParticlesUpdatedKernel);
    //UpdateSpacialLookup(spacialLookupPTUpdatedKernel, startIndicesPTUpdatedKernel, polytopeParticlesUpdatedKernel, numPT, smoothingRadius);

    CalculateCiForAllParticles(polytopeParticlesUpdatedKernel);
    CalculateCiTildaForAllParticles(spacialLookupPTUpdatedKernel, startIndicesPTUpdatedKernel, polytopeParticlesUpdatedKernel, numPT, smoothingRadius);
    CalculateGiForAllParticles(polytopeParticlesUpdatedKernel);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    // Print the duration in seconds
    std::cout << "Elapsed time: " << duration.count() << " seconds" << std::endl;
    minSV = INT16_MAX;
    maxSV = 0;
    float avgSV = 0;
    //printing all the particles ke updated kernels
    for (int i = 0; i < polytopeParticlesUpdatedKernel.size(); i++) {
        point a;
        a.x = polytopeParticlesUpdatedKernel[i].p.x;
        a.y = polytopeParticlesUpdatedKernel[i].p.y;
        a.z = polytopeParticlesUpdatedKernel[i].p.z;
        float scalarField = ScalarFieldAtPosition(a, spacialLookupPT, startIndicesPT, polytopeParticlesUpdatedKernel, numPT, smoothingRadius);
        float color = scalarField / (140000);
        DrawPoint(polytopeParticlesUpdatedKernel[i].p, color, 1, 1, 0.9f, 2);
        if (scalarField < minSV) minSV = scalarField;
        if (scalarField > maxSV) maxSV = scalarField;
        avgSV += scalarField;
    }
    avgSV = avgSV / polytopeParticlesUpdatedKernel.size();
    std::cout << "min scalar field updated kernel: " << minSV << std::endl;
    std::cout << "max scalar field updated kernel: " << maxSV << std::endl;
    std::cout << "avg scalar field updated kernel: " << avgSV << std::endl;

    //isovalue = minSV;

    //UpdateDensities3D(spacialLookupPTUpdatedKernel, startIndicesPTUpdatedKernel, polytopeParticlesUpdatedKernel, numPT, smoothingRadius);

    
    
   
 
  
    std::cout << riScale << std::endl;


    //GridIntensityUpdate(spacialLookupPT, startIndicesPT, polytopeParticlesUpdatedKernel, numPT, smoothingRadius);
    //MarchingCubes();

   


    glPopMatrix();
    glFlush();
    
    ImGui::SliderFloat("scale: ", &sc, 0, 2.0f);
    ImGui::SliderFloat("cRotx: ", &cRotx, -50.0f, 50.0f);
    //ImGui::InputFloat("cRoty: ", &cRoty);
    ImGui::SliderFloat("cRoty: ", &cRoty, -50.0f, 50.0f);
    //ImGui::InputFloat("cRotz: ", &cRotz);
    ImGui::SliderFloat("cRotz: ", &cRotz, -50.0f, 50.0f);
    //index3D++;
    count++; //for writing into the .ply file
    //writing the polytope into a .ply
    happly::PLYData plyOut;
    std::vector<std::array<double, 3>> meshVertexPositions;
    std::vector<std::vector<size_t>> meshFaceIndices;

   // std::string fileNamecsv = GetFileName();
   // fileNamecsv += ".csv";
   // std::fstream myFilePT;
   // myFilePT.open(fileNamecsv, std::ios::out); //write mode
   //if (myFilePT.is_open()) {
   //    myFilePT << "x Coord" << "," << "y Coord" << "," << "z coord"<<","<<"num nbhd particles"<<","<<"scalar value"<<","<<"||CTilda||"<<","<<"||KsCsTilda||" << std::endl;
   //    for (int i = 0; i < numPT; i++) {
   //        std::vector<int> nbhdParticels = GetNbhdPartilcesBruteForce(polytopeParticles[i], polytopeParticles, smoothingRadius); int size = nbhdParticels.size();
   //        point a;
   //        a.x = polytopeParticles[i].p.x;
   //        a.y = polytopeParticles[i].p.y;
   //        a.z = polytopeParticles[i].p.z;
   //        float scalar = ScalarFieldAtPosition(a, spacialLookupPT, startIndicesPT, polytopeParticles, numPT, smoothingRadius);
   //        float normCs = polytopeParticles[i].CiTilda.norm();
   //        float normKsCs = (polytopeParticles[i].CiTilda * Ks).norm();
   //        myFilePT << polytopeParticles[i].p.x << "," << polytopeParticles[i].p.y << "," << polytopeParticles[i].p.z<<","<<size<<","<<scalar<<","<<normCs<<","<<normKsCs << std::endl;
   //    }
   //    myFilePT.close();
   //}
   //std::string fileNamecsvUK = GetFileName();
   //fileNamecsvUK += "UpdatedKernel";
   //fileNamecsvUK += ".csv";
   //std::fstream myFilePTUK;
   //myFilePTUK.open(fileNamecsvUK, std::ios::out); //write mode
   //if (myFilePTUK.is_open()) {
   //    myFilePTUK << "x Coord" << "," << "y Coord" << "," << "z coord" << "," << "num nbhd particles" << "," << "scalar value" << "," << "||CTilda||" << "," << "||KsCsTilda||" << std::endl;
   //    for (int i = 0; i < numPT; i++) {
   //        std::vector<int> nbhdParticels = GetNbhdPartilcesBruteForce(polytopeParticlesUpdatedKernel[i], polytopeParticlesUpdatedKernel, smoothingRadius); int size = nbhdParticels.size();
   //        point a;
   //        a.x = polytopeParticlesUpdatedKernel[i].p.x;
   //        a.y = polytopeParticlesUpdatedKernel[i].p.y;
   //        a.z = polytopeParticlesUpdatedKernel[i].p.z;
   //        float scalar = ScalarFieldAtPosition(a, spacialLookupPTUpdatedKernel, startIndicesPTUpdatedKernel, polytopeParticlesUpdatedKernel, numPT, smoothingRadius);
   //        float normCs = polytopeParticlesUpdatedKernel[i].CiTilda.norm();
   //        float normKsCs = (polytopeParticlesUpdatedKernel[i].CiTilda * Ks).norm();
   //        myFilePTUK << polytopeParticlesUpdatedKernel[i].p.x << "," << polytopeParticlesUpdatedKernel[i].p.y << "," << polytopeParticlesUpdatedKernel[i].p.z << "," << size << "," << scalar << "," << normCs << "," << normKsCs << std::endl;
   //    }
   //    myFilePTUK.close();
   //}
   //

   /*if(count == 1)
   {
       meshVertexPositions.clear();
       meshFaceIndices.clear();
       std::cout << "get surface mesh is started\n";
       GetSurfaceMesh(meshVertexPositions,meshFaceIndices);
       std::cout << "get surface mesh is done\n";
       plyOut.addVertexPositions(meshVertexPositions);
       std::cout << "add vertex positions done\n";
       plyOut.addFaceIndices(meshFaceIndices);
       std::cout << "add mesh positions done\n";
       std::string fileName = GetFileName();
       fileName += ".ply";
       plyOut.write(fileName, happly::DataFormat::ASCII);
       std::cout << ".ply file ready \n";
       count = 0;
       exit(90);
   }*/

}
