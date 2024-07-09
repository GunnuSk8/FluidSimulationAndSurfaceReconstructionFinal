#pragma once
#include "FluidParticlesProperties&Init.h"
#include "SmoothingFunctions.h"
#include "MaRchingSQuares.h"
#include "MarRChingCBueS.h"
#include "SurfaceReconstructionComputations.h"
#include "SurfaceReconstructionComputations2D.h"
#include "SurfaceReconstructionComputations3D.h"

static void SetParameters() {

	smoothingRadius = 0.01;

	sigma = 1;

	rez2D = 120;
	isoValue2D = 0;

	rez = 100;
	isovalue = 32767/*20836.75*//*20000*//*21673.5*//*19951*/;

	lambda = 0.9;
	riScale = 6;
	ri = riScale * smoothingRadius;

	Kr2D = 12;
	Ks2D = 5400;
	Kn2D = 5;
	Ne2D = 3;

	Kr = 8;
	Ks = 10;
	Kn = 5;
	Ne = 1;
}