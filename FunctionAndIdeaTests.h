#pragma once


static int numNbhd = 5000;
std::vector<particle> checkNbhdParticlesFunc(numNbhd);
std::vector<particleIndexCellKeyPair> spatialLookupcheckNbhdParticlesFunc(numNbhd);
std::vector<int> startIndicescheckNbhdParticlesFunc(numNbhd);
void initnbhdCheck() {
    checkNbhdParticlesFunc[0].p.x = 0; checkNbhdParticlesFunc[0].p.y = 0; checkNbhdParticlesFunc[0].p.z = 0;
    for (int i = 1; i < numNbhd; i++) {
        checkNbhdParticlesFunc[i].p.x = (float)i / 50;

        checkNbhdParticlesFunc[i].p.y = (float)((double)std::rand() / (double)RAND_MAX - 0.5) * 5 / 5;
        checkNbhdParticlesFunc[i].p.z = (float)((double)std::rand() / (double)RAND_MAX - 0.5) * 5 / 5;
        checkNbhdParticlesFunc[i].p.x = (float)((double)std::rand() / (double)RAND_MAX - 0.5) * 5 / 5;
    }
    for (int i = 0; i < numNbhd; i++) {
        checkNbhdParticlesFunc[i].particleIndex = i;
    }
}
static float x1 = 0, y11 = 0, z1 = 0;
static float x2 = 0, y2 = 0, z2 = 0;
static float x3 = 0, y3 = 0, z3 = 0;
void displayNbhdCheck() {
    static float cRotx = 20;
    static float cRoty = 20;
    static float cRotz = 20;
    glLoadIdentity();
    glPushMatrix();
    glTranslatef(0.0, 0.0, 0.0);
    static float sc = 0.2;
    glScalef(sc, sc, sc);
    glRotatef(cRotx, cRoty, cRotz, 0.0);
    DrawBoundingBox3D();

    checkNbhdParticlesFunc[1].p.x = x1; checkNbhdParticlesFunc[1].p.y = y11; checkNbhdParticlesFunc[1].p.z = z1;
    checkNbhdParticlesFunc[2].p.x = x2; checkNbhdParticlesFunc[2].p.y = y2; checkNbhdParticlesFunc[2].p.z = z2;
    checkNbhdParticlesFunc[3].p.x = x3; checkNbhdParticlesFunc[3].p.y = y3; checkNbhdParticlesFunc[3].p.z = z3;

    //UpdateParticlesPredictedPosition(checkNbhdParticlesFunc);
    //UpdateSpacialLookup(spatialLookupcheckNbhdParticlesFunc, startIndicescheckNbhdParticlesFunc, checkNbhdParticlesFunc, numNbhd, smoothingRadius);

    //std::vector<int> nbhdPar = GetNbhdParticlesForParticle3D(checkNbhdParticlesFunc[index3D], spatialLookupcheckNbhdParticlesFunc, startIndicescheckNbhdParticlesFunc,
      //  checkNbhdParticlesFunc, numNbhd, smoothingRadius);
    std::vector<int> nbhdPar = GetNbhdPartilcesBruteForce(checkNbhdParticlesFunc[index3D],
        checkNbhdParticlesFunc,smoothingRadius);


    std::cout << nbhdPar.size() << std::endl;

    /*for (int i = 0; i < spatialLookupcheckNbhdParticlesFunc.size(); i++) {
        std::cout << spatialLookupcheckNbhdParticlesFunc[i].particleIndex << " , ";
    }
    std::cout << std::endl;
    std::cin.get();
    std::cout << "startind: \n";
    for (int i = 0; i < startIndicescheckNbhdParticlesFunc.size(); i++) {
        std::cout << startIndicescheckNbhdParticlesFunc[i] << " , ";
    }
    std::cout << std::endl;
    std::cin.get();*/


    for (int i = 0; i < checkNbhdParticlesFunc.size(); i++) {
        DrawPoint(checkNbhdParticlesFunc[i].p,1,1,1,0.4f,2);
        /*std::vector<int> nbhdParLocal = GetNbhdParticlesForParticle3D(checkNbhdParticlesFunc[i], spatialLookupcheckNbhdParticlesFunc, startIndicescheckNbhdParticlesFunc,
            checkNbhdParticlesFunc, numNbhd, smoothingRadius);
        if (nbhdParLocal.size() == 0) DrawPoint(checkNbhdParticlesFunc[i].p, 1, 0, 0, 0.4f, 2);*/
    }
   

    DrawPoint(checkNbhdParticlesFunc[index3D].p, 1, 0, 0, 1,5);
    DrawSphere(checkNbhdParticlesFunc[index3D].p);
    for (int i = 0; i < nbhdPar.size(); i++) {
        DrawPoint(checkNbhdParticlesFunc[nbhdPar[i]].p, 0, 1, 1, 1,3);
    }
   

    glPopMatrix();
    glFlush();

    //ImGui::InputFloat("cRotx: ", &cRotx);
    ImGui::SliderFloat("h: ", &smoothingRadius, -1.0f, 2.0f);
    ImGui::SliderFloat("scale: ", &sc, 0, 2.0f);
    ImGui::SliderFloat("cRotx: ", &cRotx, -50.0f, 50.0f);
    //ImGui::InputFloat("cRoty: ", &cRoty);
    ImGui::SliderFloat("cRoty: ", &cRoty, -50.0f, 50.0f);
    //ImGui::InputFloat("cRotz: ", &cRotz);
    ImGui::SliderFloat("cRotz: ", &cRotz, -50.0f, 50.0f);
    ImGui::SliderFloat("x1: ", &x1, -0.4f, 0.4f);
    ImGui::SliderFloat("y1: ", &y11, -0.4f, 0.4f);
    ImGui::SliderFloat("z1: ", &z1, -0.4f, 0.4f);
    ImGui::SliderFloat("x2: ", &x2, -0.4f, 0.4f);
    ImGui::SliderFloat("y2: ", &y2, -0.4f, 0.4f);
    ImGui::SliderFloat("z2: ", &z2, -0.4f, 0.4f);
    ImGui::SliderFloat("x3: ", &x3, -0.4f, 0.4f);
    ImGui::SliderFloat("y3: ", &y3, -0.4f, 0.4f);
    ImGui::SliderFloat("z3: ", &z3, -0.4f, 0.4f);
    ImGui::InputInt("index: ", &index3D);
    
}





//to test the gradient approach to try to remove the inner surface that we get for static models
std::vector<particle> particles3DGradientTest(totalNumberOfParticles);
std::vector<particleIndexCellKeyPair> spacialLookup3DGradientTest(totalNumberOfParticles);
std::vector<int> startIndices3DGradientTest(totalNumberOfParticles);

static void GradientTest() {
    static float cRot = 20;
    glLoadIdentity();
    glPushMatrix();
    glTranslatef(0.0, 0.0, 0.0);
    glScalef(1, 1, 1);
    glRotatef(cRot, 1.0, 1.0, 0.0);
    //if i run the below commented code then the speed of the particles keep on increasing and the surface changes for some weird reason, i dont know what is happening
    /*UpdateParticlesPredictedPosition(particles3DGradientTest);
    UpdateSpacialLookup(spacialLookup3DGradientTest, startIndices3DGradientTest, particles3DGradientTest, totalNumberOfParticles, smoothingRadius);
    UpdateDensities3D(spacialLookup3DGradientTest, startIndices3DGradientTest, particles3DGradientTest, totalNumberOfParticles, smoothingRadius);*/
    GridIntensityUpdate(spacialLookup3DGradientTest, startIndices3DGradientTest, particles3DGradientTest, totalNumberOfParticles, smoothingRadius);
    MarchingCubes();
    DrawParticles(particles3DGradientTest, totalNumberOfParticles);


    cRot += 01;
    glPopMatrix();
    glFlush();
}


//ellipse testing
static float a = 0.3f, b = 0, c = 0, d = 0.5f;
static float xpos = 0.1f, ypos = 0.1f;


static void DrawEllipse2D(float cx, float cy, float rx, float ry, int num_segments)
{
    float theta = 2 * 3.1415926f / float(num_segments);
    float c = cosf(theta);//precalculate the sine and cosine
    float s = sinf(theta);
    float t;

    float x = 1;//we start at angle = 0 
    float y = 0;

    glBegin(GL_LINE_LOOP);
    for (int ii = 0; ii < num_segments; ii++)
    {
        //apply radius and offset
        glVertex2f(x * rx + cx, y * ry + cy);//output vertex 

        //apply the rotation matrix
        t = x;
        x = c * x - s * y;
        y = s * t + c * y;
    }
    glEnd();
}

static void DrawEllipse() {
    Eigen::Matrix<float, 2, 2> ellipseMatrix;
    Eigen::Matrix<float, 1, 2> pos;
    ellipseMatrix << a, b,
        c, d;
    Vector3D origin, pt; pt.x = xpos, pt.y = ypos;
    for (int x = 0; x < rez2D; x++) {
        for (int y = 0; y < rez2D; y++) {
            pos(0, 0) = gridPt2D[x][y].x - pt.x;
            pos(0, 1) = gridPt2D[x][y].y - pt.y;
            Eigen::MatrixXf answer = pos * ellipseMatrix * pos.transpose();
            gridPt2D[x][y].intensity = answer(0, 0);
            //gridPt2D[x][y].intensity = function(gridPt2D[x][y].x, gridPt2D[x][y].y, gridPt2D[x][y].z);

        }
    }

    DrawPoint(origin);
    DrawPoint(pt, 1, 0, 0);
    MarchingSquares();
}
static float a3D = 0.3f, b3D = 0, c3D = 0,
d3D = 0, e3D = 0.3f, f3D = 0,
g3D = 0, h3D = 0, i3D = 0.3f;
static float xpos3D = 0.1f, ypos3D = 0.1f, zpos3D = 0.1f;
static void DrawEllipse3D() {
    static float cRot = 20;
    glLoadIdentity();
    glPushMatrix();
    glTranslatef(0.0, 0.0, 0.0);
    glScalef(1, 1, 1);
    glRotatef(cRot, 1.0, 1.0, 0.0);

    Eigen::Matrix<float, 3, 3> ellipseMatrix;
    Eigen::Matrix<float, 1, 3> pos;
    ellipseMatrix << a3D, b3D, c3D,
        d3D, e3D, f3D,
        g3D, h3D, i3D;
    Vector3D origin, pt; pt.x = xpos3D, pt.z = zpos3D;
    for (int x = 0; x < rez; x++) {
        for (int y = 0; y < rez; y++) {
            for (int z = 0; z < rez; z++) {
                pos(0, 0) = gridPt[x][y][z].x - pt.x;
                pos(0, 1) = gridPt[x][y][z].y - pt.y;
                pos(0, 2) = gridPt[x][y][z].z - pt.z;
                Eigen::MatrixXf answer = pos * ellipseMatrix * pos.transpose();
                gridPt[x][y][z].intensity = answer(0, 0);


            }
        }

    }


    DrawPoint(origin);
    DrawPoint(pt, 1, 0, 0);
    MarchingCubes();




    cRot += 0.01f;
    glPopMatrix();
    glFlush();
}


//marching cubes test
static void MC() {
    static float cRot = 20;
    glLoadIdentity();
    glPushMatrix();
    glTranslatef(0.0, 0.0, 0.0);
    glScalef(1, 1, 1);
    glRotatef(cRot, 1.0, 1.0, 0.0);

    MarchingCubes();

    cRot += 01;
    glPopMatrix();
    glFlush();
}
