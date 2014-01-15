
#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

#include "camera.h"
#include "paramset.h"
#include "film.h"
using namespace std;

struct Len {
	bool isStop;
	// z: the location of the surface
	// aperture: it is the half of the length in the file
	float z, radius, n, aperture;
}; 

// RealisticCamera Declarations
class RealisticCamera : public Camera {
public:
	// RealisticCamera Public Methods
	RealisticCamera(const AnimatedTransform &cam2world,
						float hither, float yon, float sopen,
						float sclose, float filmdistance, float aperture_diameter, string specfile,
						float filmdiag, Film *film);
	float GenerateRay(const CameraSample &sample, Ray *) const;
  
private:
	float Hither, Yon, distance;
	vector<Len> lens;
	Transform RasterToCamera;
	bool ParseSpec(const string &file);
	bool Propagating(Point& O, Vector& D, Len surface, float n2) const;
};


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);


#endif	// PBRT_CAMERAS_REALISTIC_H