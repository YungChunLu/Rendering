// lights/MedianCutEnvironmentLight.cpp*
#include "stdafx.h"
#include "lights/MedianCutEnvironmentLight.h"
#include "sh.h"
#include "montecarlo.h"
#include "paramset.h"
#include "imageio.h"

// MedianCutEnvironmentLight Utility Classes
struct infiniteAreaCube {
    // MedianCutEnvironmentLight Public Methods
    infiniteAreaCube(const MedianCutEnvironmentLight *l, const Scene *s,
                     float t, bool cv, float pe)
        : light(l), scene(s), time(t), pEpsilon(pe), computeVis(cv) { }
    Spectrum operator()(int, int, const Point &p, const Vector &w) {
        Ray ray(p, w, pEpsilon, INFINITY, time);
        if (!computeVis || !scene->IntersectP(ray))
            return light->Le(RayDifferential(ray));
        return 0.f;
    }
    const MedianCutEnvironmentLight *light;
    const Scene *scene;
    float time, pEpsilon;
    bool computeVis;
};

// MedianCutEnvironmentLight Method Definitions
MedianCutEnvironmentLight::~MedianCutEnvironmentLight() {
    delete distribution;
    delete radianceMap;
}


MedianCutEnvironmentLight::MedianCutEnvironmentLight(const Transform &light2world,
        const Spectrum &L, int ns, const string &texmap)
    : Light(light2world, ns) {
    int width = 0, height = 0;
    RGBSpectrum *texels = NULL;
    // Read texel data from _texmap_ into _texels_
    if (texmap != "") {
        texels = ReadImage(texmap, &width, &height);
        if (texels)
            for (int i = 0; i < width * height; ++i)
                texels[i] *= L.ToRGBSpectrum();
    }
    if (!texels) {
        width = height = 1;
        texels = new RGBSpectrum[1];
        texels[0] = L.ToRGBSpectrum();
    }
    radianceMap = new MIPMap<RGBSpectrum>(width, height, texels);

	float SolidAngle = ((2.f * M_PI) / (width - 1)) * (M_PI / (1.f * (height - 1)));

    // Compute scalar-valued image _img_ from environment map
    float filter = 1.f / max(width, height);
    float *img = new float[width*height];
    for (int v = 0; v < height; ++v) {
        float vp = (float)v / (float)height;
        float sinTheta = sinf(M_PI * float(v+.5f)/float(height));
        for (int u = 0; u < width; ++u) {
            float up = (float)u / (float)width;
            img[u+v*width] = radianceMap->Lookup(up, vp, filter).y();
            img[u+v*width] *= sinTheta;
			texels[u+v*width] *= (SolidAngle * sinTheta); 
        }
    }

	// Setting the table
	int bound = width * height;	
	float *table = new float[bound];
	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j){
			int me = i * width + j, up = me - width, left = me - 1;
			if (i != 0 && j != 0)
				table[me] = table[up] + table [left] - table[up - 1] + img[me] * SolidAngle;
			else
				if (i == 0 && j == 0)
					table[me] = img[me] * SolidAngle;
					
				else if (i == 0)
					table[me] = table[left] + img[me] * SolidAngle;
				else
					table[me] = table[up] + img[me] * SolidAngle;
		}
	}

	float u_unit =  1.f / (width - 1), v_unit = 1.f / (height - 1);

	// Add the entire probe image
	int leftcorner[2] = {0, 0};
	int l = width - 1, w = height - 1;
	float center[2] = {0.5 * l * u_unit, 0.5 * w * v_unit};
	float halfenergy = 0.5 * table[(leftcorner[0] + l) + (leftcorner[1] + w) * width];
	vector<Region> regions;
	regions.push_back(Region(leftcorner, l, w, center, halfenergy));

	// Partition
	int partition_times = 6;
	for(int n = 0; n < partition_times; ++n){
		vector<Region> region_buffer;
		for (vector<Region>::iterator it = regions.begin(); it != regions.end(); ++it){
			Region region = *it;
			int step, partition_point, partition_point_x, partition_point_y;
			float partition_energy, E1 = 0.f, E2 = 0.f, E3 = 0.f;

			if (*(region.LeftCorner) != 0 && *(region.LeftCorner+1) != 0){
				E1 = table[(*(region.LeftCorner) - 1 + (*(region.LeftCorner+1) - 1) * width)];
				E2 = table[(*(region.LeftCorner) + region.L + (*(region.LeftCorner+1) - 1) * width)];
				E3 = table[(*(region.LeftCorner) - 1 + (*(region.LeftCorner+1) + region.W) * width)];}
			else if ((*(region.LeftCorner) == 0) && (*(region.LeftCorner+1) != 0)){
				E2 = table[(region.L + (*(region.LeftCorner+1) - 1) * width)];}
			else if ((*(region.LeftCorner) != 0) && (*(region.LeftCorner+1) == 0)){
				E3 = table[(*(region.LeftCorner) - 1 + region.W * width)];}

			if (region.IsLatitude){
				step = width;
				partition_point_x = *(region.LeftCorner) + region.L;
				partition_point_y = int(*(region.LeftCorner+1) + 0.5 * region.W);
				if (*(region.LeftCorner) != 0) E3 = table[*(region.LeftCorner) - 1 + partition_point_y * width];}
			else {
				step = 1;
				partition_point_x = int(*(region.LeftCorner) + 0.5 * region.L);
				partition_point_y = *(region.LeftCorner+1) + region.W;
				if (*(region.LeftCorner+1) != 0) E2 = table[partition_point_x + (*(region.LeftCorner+1) - 1) * width];}

			partition_point = partition_point_x + partition_point_y * width;
			partition_energy = table[partition_point] - E2 - E3 + E1;

			if (partition_energy > region.HalfEnergy) step *= -1;

			while(1){
				int next_index = partition_point + step;
				if (region.IsLatitude && *(region.LeftCorner) != 0)
					E3 = table[next_index - region.L - 1];
				else if (!region.IsLatitude && *(region.LeftCorner+1) != 0)
					E2 = table[next_index - (region.W + 1) * width];
				float next_energy = table[next_index] - E2 - E3 + E1;
				// Check if we find the right partition point
				if ((partition_energy - region.HalfEnergy) * (next_energy - region.HalfEnergy) < 0) break;
				partition_point = next_index;
				partition_energy = next_energy;
			}

			int u = partition_point % width;
			int v = (partition_point - u) / width;
			float halfenergy1 = partition_energy * 0.5, halfenergy2 = region.HalfEnergy - halfenergy1;
			int leftcorner1[2] = {*(region.LeftCorner), *(region.LeftCorner+1)};
			int leftcorner2[2] = {0, 0};
			int new_L1 = region.L, new_L2 = region.L, new_W1 = region.W, new_W2 = region.W;

			if (region.IsLatitude){
				new_W1 = v - *(region.LeftCorner+1);
				new_W2 = region.W - new_W1 - 1;
				leftcorner2[0] = *(region.LeftCorner);
				leftcorner2[1] = v + 1;}
			else {
				new_L1 = u - *(region.LeftCorner);
				new_L2 = region.L - new_L1 - 1;
				leftcorner2[0] = u + 1;
				leftcorner2[1] = *(region.LeftCorner+1);}

			float center1[2] = {(leftcorner1[0] * 1.f + 0.5 * new_L1) * u_unit, (leftcorner1[1] * 1.f + 0.5 * new_W1) * v_unit},
				  center2[2] = {(leftcorner2[0] * 1.f + 0.5 * new_L2) * u_unit, (leftcorner2[1] * 1.f + 0.5 * new_W2) * v_unit};
			
			region_buffer.push_back(Region(leftcorner1, new_L1, new_W1, center1, halfenergy1));
			region_buffer.push_back(Region(leftcorner2, new_L2, new_W2, center2, halfenergy2));
		}
		regions.clear();
		regions = region_buffer;
	}

	// Assign VPLs
	PDF = 1.f / regions.size();
	for (vector<Region>::iterator it = regions.begin(); it != regions.end(); ++it){
		Region region = *it;
		RGBSpectrum spectrum = RGBSpectrum(0.f);
		int index = *(region.LeftCorner+1) * width + *(region.LeftCorner);
		for(int v = 0; v <= region.W; ++v){
			for(int u = 0; u <= region.L; ++u){
				spectrum += texels[index + (u + v * width)];
			}
		}
		VPLs.push_back(VPL(region.Center, spectrum));
	}

	// Clean useless data
	delete[] texels;
	delete table;

    // Compute sampling distributions for rows and columns of image
    distribution = new Distribution2D(img, width, height);
    delete[] img;
}


Spectrum MedianCutEnvironmentLight::Power(const Scene *scene) const {
    Point worldCenter;
    float worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    return M_PI * worldRadius * worldRadius *
        Spectrum(radianceMap->Lookup(.5f, .5f, .5f), SPECTRUM_ILLUMINANT);
}


Spectrum MedianCutEnvironmentLight::Le(const RayDifferential &r) const {
    Vector wh = Normalize(WorldToLight(r.d));
    float s = SphericalPhi(wh) * INV_TWOPI;
    float t = SphericalTheta(wh) * INV_PI;
    return Spectrum(radianceMap->Lookup(s, t), SPECTRUM_ILLUMINANT);
}


void MedianCutEnvironmentLight::SHProject(const Point &p, float pEpsilon,
        int lmax, const Scene *scene, bool computeLightVis,
        float time, RNG &rng, Spectrum *coeffs) const {
    // Project _InfiniteAreaLight_ to SH using Monte Carlo if visibility needed
    if (computeLightVis) {
        Light::SHProject(p, pEpsilon, lmax, scene, computeLightVis,
                         time, rng, coeffs);
        return;
    }
    for (int i = 0; i < SHTerms(lmax); ++i)
        coeffs[i] = 0.f;
    int ntheta = radianceMap->Height(), nphi = radianceMap->Width();
    if (min(ntheta, nphi) > 50) {
        // Project _InfiniteAreaLight_ to SH from lat-long representation

        // Precompute $\theta$ and $\phi$ values for lat-long map projection
        float *buf = new float[2*ntheta + 2*nphi];
        float *bufp = buf;
        float *sintheta = bufp;  bufp += ntheta;
        float *costheta = bufp;  bufp += ntheta;
        float *sinphi = bufp;    bufp += nphi;
        float *cosphi = bufp;
        for (int theta = 0; theta < ntheta; ++theta) {
            sintheta[theta] = sinf((theta + .5f)/ntheta * M_PI);
            costheta[theta] = cosf((theta + .5f)/ntheta * M_PI);
        }
        for (int phi = 0; phi < nphi; ++phi) {
            sinphi[phi] = sinf((phi + .5f)/nphi * 2.f * M_PI);
            cosphi[phi] = cosf((phi + .5f)/nphi * 2.f * M_PI);
        }
        float *Ylm = ALLOCA(float, SHTerms(lmax));
        for (int theta = 0; theta < ntheta; ++theta) {
            for (int phi = 0; phi < nphi; ++phi) {
                // Add _InfiniteAreaLight_ texel's contribution to SH coefficients
                Vector w = Vector(sintheta[theta] * cosphi[phi],
                                  sintheta[theta] * sinphi[phi],
                                  costheta[theta]);
                w = Normalize(LightToWorld(w));
                Spectrum Le = Spectrum(radianceMap->Texel(0, phi, theta),
                                       SPECTRUM_ILLUMINANT);
                SHEvaluate(w, lmax, Ylm);
                for (int i = 0; i < SHTerms(lmax); ++i)
                    coeffs[i] += Le * Ylm[i] * sintheta[theta] *
                        (M_PI / ntheta) * (2.f * M_PI / nphi);
            }
        }

        // Free memory used for lat-long theta and phi values
        delete[] buf;
    }
    else {
        // Project _InfiniteAreaLight_ to SH from cube map sampling
        SHProjectCube(infiniteAreaCube(this, scene, time, computeLightVis,
                                       pEpsilon),
                      p, 200, lmax, coeffs);
    }
}


MedianCutEnvironmentLight *CreateMedianCutEnvironmentLight(const Transform &light2world,
        const ParamSet &paramSet) {
    Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    string texmap = paramSet.FindOneFilename("mapname", "");
    int nSamples = paramSet.FindOneInt("nsamples", 1);
    if (PbrtOptions.quickRender) nSamples = max(1, nSamples / 4);
    return new MedianCutEnvironmentLight(light2world, L * sc, nSamples, texmap);
}


Spectrum MedianCutEnvironmentLight::Sample_L(const Point &p, float pEpsilon,
        const LightSample &ls, float time, Vector *wi, float *pdf,
        VisibilityTester *visibility) const {
    PBRT_INFINITE_LIGHT_STARTED_SAMPLE();
	VPL vpl = VPLs[Floor2Int(ls.uComponent * VPLs.size())];
    float theta = *(vpl.Position+1) * M_PI, phi = *(vpl.Position) * 2.f * M_PI;
    float costheta = cosf(theta), sintheta = sinf(theta);
    float sinphi = sinf(phi), cosphi = cosf(phi);
    *wi = LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
                              costheta));

    *pdf = PDF;
    Spectrum Ls = Spectrum(vpl.Spectrum, SPECTRUM_ILLUMINANT);

    // Return radiance value 
    visibility->SetRay(p, pEpsilon, *wi, time);

    PBRT_INFINITE_LIGHT_FINISHED_SAMPLE();
    return Ls;
}


float MedianCutEnvironmentLight::Pdf(const Point &, const Vector &w) const {
    PBRT_INFINITE_LIGHT_STARTED_PDF();
    Vector wi = WorldToLight(w);
    float theta = SphericalTheta(wi), phi = SphericalPhi(wi);
    float sintheta = sinf(theta);
    if (sintheta == 0.f) return 0.f;
    float p = distribution->Pdf(phi * INV_TWOPI, theta * INV_PI) /
           (2.f * M_PI * M_PI * sintheta);
    PBRT_INFINITE_LIGHT_FINISHED_PDF();
    return p;
}


Spectrum MedianCutEnvironmentLight::Sample_L(const Scene *scene,
        const LightSample &ls, float u1, float u2, float time,
        Ray *ray, Normal *Ns, float *pdf) const {
    PBRT_INFINITE_LIGHT_STARTED_SAMPLE();
    // Compute direction for infinite light sample ray

    // Find $(u,v)$ sample coordinates in infinite light texture
    float uv[2], mapPdf;
    distribution->SampleContinuous(ls.uPos[0], ls.uPos[1], uv, &mapPdf);
    if (mapPdf == 0.f) return Spectrum(0.f);

    float theta = uv[1] * M_PI, phi = uv[0] * 2.f * M_PI;
    float costheta = cosf(theta), sintheta = sinf(theta);
    float sinphi = sinf(phi), cosphi = cosf(phi);
    Vector d = -LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
                                    costheta));
    *Ns = (Normal)d;

    // Compute origin for infinite light sample ray
    Point worldCenter;
    float worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    Vector v1, v2;
    CoordinateSystem(-d, &v1, &v2);
    float d1, d2;
    ConcentricSampleDisk(u1, u2, &d1, &d2);
    Point Pdisk = worldCenter + worldRadius * (d1 * v1 + d2 * v2);
    *ray = Ray(Pdisk + worldRadius * -d, d, 0., INFINITY, time);

    // Compute _InfiniteAreaLight_ ray PDF
    float directionPdf = mapPdf / (2.f * M_PI * M_PI * sintheta);
    float areaPdf = 1.f / (M_PI * worldRadius * worldRadius);
    *pdf = directionPdf * areaPdf;
    if (sintheta == 0.f) *pdf = 0.f;
    Spectrum Ls = (radianceMap->Lookup(uv[0], uv[1]), SPECTRUM_ILLUMINANT);
    PBRT_INFINITE_LIGHT_FINISHED_SAMPLE();
    return Ls;
}


