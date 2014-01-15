
#include "stdafx.h"
#include "cameras/realistic.h"
#include <fstream>
#include <iostream>

RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
				 float hither, float yon, 
				 float sopen, float sclose, 
				 float filmdistance, float aperture_diameter, string specfile, 
				 float filmdiag, Film *f)
	: Camera(cam2world, sopen, sclose, f) // pbrt-v2 doesnot specify hither and yon
{
	distance = 0;
	Hither = hither;
	Yon = yon;
	// Parse the Specfile
	if (!ParseSpec(specfile)) {
		cerr << "Error in parsing " << specfile << endl;
		exit(1);
	}
	// now, distance represents the true location of film 
	distance -= filmdistance;
	// build the transform
	float diag = sqrtf(f->xResolution * f->xResolution + f->yResolution * f->yResolution);
	float ratio = filmdiag / diag;
	// X and Y is in the Screen Space coordinate
	float X = ratio * 0.5 * f->xResolution;
	float Y = ratio * 0.5 * f->yResolution;
	RasterToCamera = Translate(Vector(0.f, 0.f, distance)) * 
					 Translate(Vector(X, -Y, 0.f)) *
					 Scale(ratio, ratio, 1) * Scale(-1.f, 1.f, 1.f);
}

float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
	// Generate raster and camera samples
    Point Pras(sample.imageX, sample.imageY, 0);
    Point Pcamera;
    RasterToCamera(Pras, &Pcamera);

	// Set the iterator
	int it = lens.size() - 1;
	
	float z, lensU, lensV, lensZ;
	// Sample point on lens - Method1
    //ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
    //lensU *= lens[it].aperture;
    //lensV *= lens[it].aperture;	

	// Sample point on lens - Method2
	float r = lens[it].aperture * sqrtf(sample.lensU), theta = 2 * M_PI * sample.lensV;
	lensU = r * cosf(theta);
	lensV = r * sinf(theta);

	z = sqrtf(lens[it].radius * lens[it].radius - lensV * lensV - lensU * lensU);
	lensZ = lens[it].z - lens[it].radius - ((lens[it].radius < 0) ? z : -z);
	Point hit = Point(lensU, lensV, lensZ);
	Vector T = hit - Pcamera;

	// Tracing the ray	
	for(int i = it; i >= 0; --i){
		if (lens[i].isStop) {
			float deltaZ = lens[i].z - hit.z;
			T = Normalize(T);
			float t = deltaZ / T.z;
			hit = hit + t * T;
			if (hit.x * hit.x + hit.y * hit.y > lens[i].aperture * lens[i].aperture) return 0.f;
		}
		else {
			float n2 = (i == 0) ? 1 :lens[i-1].n;
			if (!Propagating(hit, T, lens[i], n2)) return 0.f;
		}
	}
	ray->o = hit;
	ray->d = Normalize(T);
	ray->mint = Hither;
	ray->maxt = (Yon - Hither) / ray->d.z;
	ray->time = Lerp(sample.time, shutterOpen, shutterClose);
	CameraToWorld(*ray, ray);
	// Set exposure weight
	float weight = Dot(Normalize(hit - Pcamera), Vector(0, 0, 1));
	weight *= weight / abs(distance);
	weight *= weight * (lens[0].aperture * lens[0].aperture * M_PI);
	return weight;
}

bool RealisticCamera::ParseSpec(const string &file){
	// open file
	const char* spec = file.c_str();
	ifstream in(spec);
	// handle open error
	if (!in) return false;
	// define size of each reading line
	int size = 256;
	char* line = new char[size];
	in.getline(line, size);
	// define the parameter we need
	float radius = 0, axpos = 0, n = 0, aperture = 0;
	while (!in.eof()){
		if(line[0] != '#'){
			Len len; 
			sscanf(line, "%f%f%f%f", &radius, &axpos, &n, &aperture);
			len.isStop = (n == 0);
			len.z = distance;
			len.radius = radius;
			len.n = (n == 0) ? 1 : n;
			len.aperture = aperture * 0.5;
			distance -= axpos;
			lens.push_back(len);
		};
		in.getline(line, size);
	}
	in.close();
	return true;
}

bool RealisticCamera::Propagating(Point& O, Vector& D, Len surface, float n2) const{
	Vector d = Normalize(D);
	Point C = Point(0.f, 0.f, surface.z - surface.radius);
	Vector OC = O - C;
	float b = Dot(OC, d);
	float c = OC.LengthSquared() - surface.radius * surface.radius;
	float determine = b * b - c;
	float t = 0;
	if (determine < 0) {
		return false;
	}
	else {
		float root = sqrtf(determine);
		t = (surface.radius > 0) ? (-b + root) : (-b - root);
	}
	O = O + t * d;
	if (surface.aperture * surface.aperture < O.y * O.y + O.x * O.x) return false;
	Vector N = (surface.radius > 0.f) ? Normalize(C - O) : Normalize(O - C);

	// Heckber's Method
	float n_ratio = surface.n / n2;
	float c1 = -Dot(d, N);
	float c2 = 1.f - n_ratio * n_ratio * (1.f - c1 * c1);
	if (c2 <= 0.f) return false;
	else c2 = sqrtf(c2);
	D = n_ratio * d + (n_ratio * c1 - c2) * N;

	return true;
}

RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {
	// Extract common camera parameters from \use{ParamSet}
	float hither = params.FindOneFloat("hither", -1);
	float yon = params.FindOneFloat("yon", -1);
	float shutteropen = params.FindOneFloat("shutteropen", -1);
	float shutterclose = params.FindOneFloat("shutterclose", -1);

	// Realistic camera-specific parameters
	string specfile = params.FindOneString("specfile", "");
	float filmdistance = params.FindOneFloat("filmdistance", 70.0); // about 70 mm default to film
 	float fstop = params.FindOneFloat("aperture_diameter", 1.0);	
	float filmdiag = params.FindOneFloat("filmdiag", 35.0);

	Assert(hither != -1 && yon != -1 && shutteropen != -1 &&
		shutterclose != -1 && filmdistance!= -1);
	if (specfile == "") {
	    Severe( "No lens spec file supplied!\n" );
	}
	return new RealisticCamera(cam2world, hither, yon,
				   shutteropen, shutterclose, filmdistance, fstop, 
				   specfile, filmdiag, film);
}
