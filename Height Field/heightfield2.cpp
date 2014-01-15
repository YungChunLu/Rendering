
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// shapes/heightfield2.cpp*
#include "stdafx.h"
#include "shapes/heightfield2.h"
#include "shapes/trianglemesh.h"
#include "paramset.h"

// Heightfield2 Method Definitions
Heightfield2::Heightfield2(const Transform *o2w, const Transform *w2o,
        bool ro, int x, int y, const float *zs)
    : Shape(o2w, w2o, ro) {
    nx = x;
    ny = y;
	width = Vector(1.f / (nx-1), 1.f / (ny-1), 0);
    z = new float[nx*ny];
    memcpy(z, zs, nx*ny*sizeof(float));
	nors = new Normal[nx * ny];
	IntialNormal();
}


Heightfield2::~Heightfield2() {
    delete[] z;
	delete[] nors;
}


BBox Heightfield2::ObjectBound() const {
    float minz = z[0], maxz = z[0];
    for (int i = 1; i < nx*ny; ++i) {
        if (z[i] < minz) minz = z[i];
        if (z[i] > maxz) maxz = z[i];
    }
    return BBox(Point(0,0,minz), Point(1,1,maxz));
}


bool Heightfield2::Intersect(const Ray &r, float *tHit, float *rayEpsilon, 
		DifferentialGeometry *dg) const {
	// Transform _Ray_ to object space
    Ray ray;
    (*WorldToObject)(r, &ray);
	// Check ray against overall grid bounds
	BBox bounds = ObjectBound();
	float rayT;
	if (bounds.Inside(ray(ray.mint))) 
		rayT = ray.mint;
	else if (!bounds.IntersectP(ray, &rayT))
		return false;
	Point gridIntersect = ray(rayT);

	/* 
	1. Pos : The coordinates of the voxel currently being considered
	2. NextCrossingT : The parametric go along the ray where it makes its next crossing into
		another voxel in each of the x and y directions
	3. Step : The change in the current voxel coordinates after a step in each direction
	4. DeltaT : The distance along the ray between voxels in each direction relative to ray
	5. Out : The coordinates of the voxel after the last one the ray passes through when it exits the grid
	*/
	
	// Setup 2D DDA for ray
	int nVoxel[2] = {nx-1, ny-1};

	float NextCrossingT[2], DeltaT[2];
	int Step[2], Out[2], Pos[2];
	for (int axis = 0; axis < 2; ++axis){
		// Compute current voxel for axis
		Pos[axis] = Clamp(Float2Int(gridIntersect[axis] * nVoxel[axis]), 0, nVoxel[axis]-1);
		if (ray.d[axis] >= 0) {
			// Handle ray with positive direction for voxel stepping
			NextCrossingT[axis] = rayT + ((Pos[axis] + 1) * width[axis] - gridIntersect[axis]) / ray.d[axis];
			DeltaT[axis] = width[axis] / ray.d[axis];
			Step[axis] = 1;
			Out[axis] = nVoxel[axis];
		}else {
			// Handle ray with negative direction for voxel stepping
			NextCrossingT[axis] = rayT + (Pos[axis] * width[axis] - gridIntersect[axis]) / ray.d[axis];
			DeltaT[axis] = -width[axis] / ray.d[axis];
			Step[axis] = -1;
			Out[axis] = -1;
		}
	}
	// Walk ray through voxel grid
	Point* triangle = new Point[3];
	for (;;) {
		// Check for intersection in current voxel
		triangle[0] = Point(Pos[0] * width[0], Pos[1] * width[1], z[Pos[0] + Pos[1] * nx]);
		triangle[1] = Point((Pos[0] + 1) * width[0], Pos[1] * width[1], z[Pos[0] + 1 + Pos[1] * nx]);
		triangle[2] = Point(Pos[0] * width[0], (Pos[1] + 1) * width[1], z[Pos[0] + (Pos[1] + 1) * nx]);
		if (VoxelIntersect(r, tHit, rayEpsilon, triangle, dg))
			return true;
		triangle[0] = Point((Pos[0] + 1) * width[0], (Pos[1] + 1) * width[1], z[Pos[0] + 1 + (Pos[1] + 1) * nx]);
		if (VoxelIntersect(r, tHit, rayEpsilon, triangle, dg))
			return true;
		// Find stepAxis for stepping to next voxel
		int stepAxis = NextCrossingT[0] < NextCrossingT[1] ? 0 : 1;
		if (ray.maxt < NextCrossingT[stepAxis])
			break;
		Pos[stepAxis] += Step[stepAxis];
        if (Pos[stepAxis] == Out[stepAxis])
			break;
		NextCrossingT[stepAxis] += DeltaT[stepAxis];
	}
	return false;
};

bool Heightfield2::IntersectP(const Ray &r) const {
	// Transform _Ray_ to object space
    Ray ray;
    (*WorldToObject)(r, &ray);
	// Check ray against overall grid bounds
	BBox bounds = ObjectBound();
	float rayT;
	if (bounds.Inside(ray(ray.mint))) 
		rayT = ray.mint;
	else if (!bounds.IntersectP(ray, &rayT))
		return false;
	Point gridIntersect = ray(rayT);

	// Setup 2D DDA for ray
	int nVoxel[2] = {nx-1, ny-1};

	float NextCrossingT[2], DeltaT[2];
	int Step[2], Out[2], Pos[2];
	for (int axis = 0; axis < 2; ++axis){
		// Compute current voxel for axis
		Pos[axis] = Clamp(Float2Int(gridIntersect[axis] * nVoxel[axis]), 0, nVoxel[axis] - 1);
		if (ray.d[axis] >= 0) {
			// Handle ray with positive direction for voxel stepping
			NextCrossingT[axis] = rayT + ((Pos[axis] + 1) * width[axis] - gridIntersect[axis]) / ray.d[axis];
			DeltaT[axis] = width[axis] / ray.d[axis];
			Step[axis] = 1;
			Out[axis] = nVoxel[axis];
		}else {
			// Handle ray with negative direction for voxel stepping
			NextCrossingT[axis] = rayT + (Pos[axis] * width[axis] - gridIntersect[axis]) / ray.d[axis];
			DeltaT[axis] = -width[axis] / ray.d[axis];
			Step[axis] = -1;
			Out[axis] = -1;
		}
	}
	// Walk ray through voxel grid
	Point* triangle = new Point[3];
	for (;;) {
		// Check for intersection in current voxel
		triangle[0] = Point(Pos[0] * width[0], Pos[1] * width[1], z[Pos[0] + Pos[1] * nx]);
		triangle[1] = Point((Pos[0] + 1) * width[0], Pos[1] * width[1], z[Pos[0] + 1 + Pos[1] * nx]);
		triangle[2] = Point(Pos[0] * width[0], (Pos[1] + 1) * width[1], z[Pos[0] + (Pos[1] + 1) * nx]);
		if (VoxelIntersect(r, triangle)) 
			return true;
		triangle[0] = Point((Pos[0] + 1) * width[0], (Pos[1] + 1) * width[1], z[Pos[0] + 1 + (Pos[1] + 1) * nx]);
		if (VoxelIntersect(r, triangle)) 
			return true;
		// Find stepAxis for stepping to next voxel
		int stepAxis = (NextCrossingT[0] < NextCrossingT[1]) ? 0 : 1;
		if (ray.maxt < NextCrossingT[stepAxis])
			break;
		Pos[stepAxis] += Step[stepAxis];
        if (Pos[stepAxis] == Out[stepAxis])
			break;
		NextCrossingT[stepAxis] += DeltaT[stepAxis];
	}
	return false;
}

void Heightfield2::GetShadingGeometry(const Transform &obj2world,
        const DifferentialGeometry &dg,
		DifferentialGeometry *dgShading) const {
			int x = Clamp(Float2Int(dg.u * (nx -1)), 0, nx - 1), y = Clamp(Float2Int(dg.v * (ny -1)), 0, ny - 1);					
			Point p2 = Point((x + 1) * width[0], y * width[1], z[x + nx * y + 1]);
			Point p3 = Point(x * width[0], (y + 1) * width[1], z[x + nx * (y + 1)]);			
			Point p4 = Point(dg.u, p3.y - (dg.u - p3.x), 0);
			Point p1 = (dg.v <= p4.y) ? Point(x * width[0], y * width[1], z[x + nx * y]) : Point((x + 1) * width[0], (y + 1) * width[1], z[x + nx * (y + 1) + 1]);
			
			// Compute the normal at the hit point
			int q = (dg.v <= p4.y) ? (x + nx * y) : (x + 1 + nx * (y + 1));
			Normal normals[3] = {nors[q], nors[x + nx * y + 1], nors[x + nx * (y + 1)]};

			// Compute normal with barycentric method
			//float b[3];
			//float A[2][2] =
			//	{ { p2.x - p1.x, p3.x - p1.x },
			//	{ p2.y - p1.y, p3.y - p1.y } };
			//float C[2] = { dg.u - p1.x, dg.v - p1.y };
			//if (!SolveLinearSystem2x2(A, C, &b[1], &b[2]))
			//	// Handle degenerate parametric mapping
			//	 b[0] = b[1] = b[2] = 1.f/3.f;
			//else
			//	 b[0] = 1.f - b[1] - b[2];
			//Normal hitNormal = (*ObjectToWorld)(Normalize(b[0] * normals[0] + b[1] * normals[1] + b[2] * normals[2]));

			// Compute normal with Pong interpolation
			Point a = Point(p2.x, p2.y,0), b = Point(p3.x, p3.y,0), c = Point(p4.x, p4.y, 0);
			float invL = 1 / (a - b).Length();
			Normal n1 = ((c - a).Length() * normals[2] + (c - b).Length() * normals[1]) * invL,
				   n2 = ((dg.v < p4.y)) ? (fabs(dg.u - p1.x) * normals[1] + fabs(p2.x - dg.u) * normals[0]) * (nx - 1) : (fabs(dg.u - p3.x) * normals[0] + fabs(p1.x - dg.u) * normals[2]) * (nx - 1);
			float d1 = fabs(p4.y - dg.v),
				  d2 = fabs(dg.v - p1.y),
				  invd = 1.f / (d1 + d2);
			Normal hitNormal = (*ObjectToWorld)(Normalize((d2 * Normalize(n1) + d1 * Normalize(n2)) * invd));
				
			// Compute the differential normal at hit point
			Normal dndu, dndv;
			float du1 = p1.x - p3.x;
			float du2 = p2.x - p3.x;
			float dv1 = p1.y - p3.y;
			float dv2 = p2.y - p3.y;
			Normal dn1 = normals[0] - normals[2], dn2 = normals[1] - normals[2];
			float determinant = du1 * dv2 - dv1 * du2;
			if (determinant == 0.f) {
				// Handle zero determinant for triangle partial derivative matrix
				dndu = dndv = Normal(0, 0, 0);
			}
			 else {
				float invdet = 1.f / determinant;
				dndu = ( dv2 * dn1 - dv1 * dn2) * invdet;
				dndv = (-du2 * dn1 + du1 * dn2) * invdet;
			}

			Vector ss = Normalize(dg.dpdu);
			Vector ts = Cross(ss, hitNormal);
			if (ts.LengthSquared() > 0.f) {
				ts = Normalize(ts);
				ss = Cross(ts, hitNormal);
			}
			else
				CoordinateSystem((Vector)hitNormal, &ss, &ts);

			*dgShading = DifferentialGeometry(dg.p, ss, ts,
				(*ObjectToWorld)(dndu), (*ObjectToWorld)(dndv), dg.u, dg.v, dg.shape);
			dgShading->dudx = dg.dudx;  dgShading->dvdx = dg.dvdx;
			dgShading->dudy = dg.dudy;  dgShading->dvdy = dg.dvdy;
			dgShading->dpdx = dg.dpdx;  dgShading->dpdy = dg.dpdy;
}

bool Heightfield2::VoxelIntersect (const Ray &ray, float *tHit, float *rayEpsilon, Point * triangle, DifferentialGeometry *dg) const {
	const Point p1 = (*ObjectToWorld)(triangle[0]);
    const Point p2 = (*ObjectToWorld)(triangle[1]);
    const Point p3 = (*ObjectToWorld)(triangle[2]);
	Vector e1 = p2 - p1;
    Vector e2 = p3 - p1;
    Vector s1 = Cross(ray.d, e2);
    float divisor = Dot(s1, e1);
    
    if (divisor == 0.)
        return false;
    float invDivisor = 1.f / divisor;
	// Compute first barycentric coordinate
    Vector s = ray.o - p1;
    float b1 = Dot(s, s1) * invDivisor;
    if (b1 < 0. || b1 > 1.)
        return false;
    // Compute second barycentric coordinate
    Vector s2 = Cross(s, e1);
    float b2 = Dot(ray.d, s2) * invDivisor;
    if (b2 < 0. || b1 + b2 > 1.)
        return false;

    // Compute _t_ to intersection point
    float t = Dot(e2, s2) * invDivisor;
    if (t < ray.mint || t > ray.maxt)
        return false;
	
	// Compute triangle partial derivatives
    Vector dpdu, dpdv;
	
	// Compute deltas for triangle partial derivatives
    float du1 = triangle[0].x - triangle[2].x;
    float du2 = triangle[1].x - triangle[2].x;
    float dv1 = triangle[0].y - triangle[2].y;
    float dv2 = triangle[1].y - triangle[2].y;
    Vector dp1 = p1 - p3, dp2 = p2 - p3;
    float determinant = du1 * dv2 - dv1 * du2;
    if (determinant == 0.f) {
        // Handle zero determinant for triangle partial derivative matrix
        CoordinateSystem(Normalize(Cross(e2, e1)), &dpdu, &dpdv);
    }
    else {
        float invdet = 1.f / determinant;
        dpdu = ( dv2 * dp1 - dv1 * dp2) * invdet;
        dpdv = (-du2 * dp1 + du1 * dp2) * invdet;
    }

	// Interpolate $(u,v)$ triangle parametric coordinates
    float b0 = 1 - b1 - b2;
    float tu = b0*triangle[0].x + b1*triangle[1].x + b2*triangle[2].x;
    float tv = b0*triangle[0].y + b1*triangle[1].y + b2*triangle[2].y;

    // Fill in _DifferentialGeometry_ from triangle hit
    *dg = DifferentialGeometry(ray(t), dpdu, dpdv,
                               Normal(0,0,0), Normal(0,0,0),
                               tu, tv, this);
	
	*tHit = t;
    *rayEpsilon = 1e-3f * *tHit;
	return true;
}

bool Heightfield2::VoxelIntersect (const Ray &ray, Point * triangle) const {
	const Point p1 = (*ObjectToWorld)(triangle[0]);
    const Point p2 = (*ObjectToWorld)(triangle[1]);
    const Point p3 = (*ObjectToWorld)(triangle[2]);
	Vector e1 = p2 - p1;
    Vector e2 = p3 - p1;
    Vector s1 = Cross(ray.d, e2);
    float divisor = Dot(s1, e1);
    
    if (divisor == 0.)
        return false;
    float invDivisor = 1.f / divisor;
	// Compute first barycentric coordinate
    Vector s = ray.o - p1;
    float b1 = Dot(s, s1) * invDivisor;
    if (b1 < 0. || b1 > 1.)
        return false;
    // Compute second barycentric coordinate
    Vector s2 = Cross(s, e1);
    float b2 = Dot(ray.d, s2) * invDivisor;
    if (b2 < 0. || b1 + b2 > 1.)
        return false;

    // Compute _t_ to intersection point
    float t = Dot(e2, s2) * invDivisor;
    if (t < ray.mint || t > ray.maxt)
        return false;
	return true;
}

void Heightfield2::IntialNormal () const{
	// First Method of assigning normals
	// Two different searching order of finding normal
	const int r1[12] = {-1, 0, -1, 1, 0, 1, 1, 0, 1, -1, 0, -1};
	const int r2[12] = {1, 0, 1, -1, 0, -1, -1, 0, -1, 1, 0, 1};
	const int *v = r1;
	for (int j = 0; j < ny; ++j){
		for (int i = 0; i < nx; ++i){
			if (((i == nx - 1) && (j > 0 && j <= ny -1)) || (((j == ny - 1) && (i > 0 && i < nx - 1))))
				v = r2;
			else
				v = r1;

			vector<Point> points;
			Normal result = Normal();
			// The point at which the normal we want to find stays
			Point p = Point(i * width[0], j * width[1], z[i + nx * j]);

			for (int k = 0; k < 12; k += 2){
				int px = i + v[k], py = j + v[k+1], pz = px + nx * py;
				if ((px >= 0) && (px <= (nx - 1)) && (py >= 0) && (py <= (ny - 1))){
					Point point = Point(px * width[0], py * width[1], z[pz]);
					points.push_back(point);
				}
			}
			for (vector<Point>::iterator it = points.begin(); (it + 1) != points.end(); it++){
				Vector v1 = *it - p, v2 = *(it+1) - p;
				result += Normalize(Normal(Cross(v1, v2)));
			}
			if (points.size() == 6) {
				vector<Point>::iterator it = points.begin();
				Vector v2 = *it - p, v1 = *(it+5) - p;
				result += Normalize(Normal(Cross(v1, v2)));
				nors[i + nx * j] = Normalize(result / 6);
			}
			else nors[i + nx * j] = Normalize(result / (points.size() - 1));
		}
	}

	// Second Method of assigning normals
	//for (int j = 0; j < ny; ++j){
	//	for (int i = 0; i < nx; ++i){
	//		Point p = Point(i * width[0], j * width[1], z[i + nx * j]);
	//		Normal result = Normal();
	//		if ((i == 0) && (j == 0)){
	//			result = Normalize(Normal((0.5*width[1]*z[i+1+nx*j]), (0.5*width[0]*z[i+nx*(j+1)]), -1*width[0]*width[1]));
	//		}
	//		else if ((i == nx-1) && (j == 0)){
	//			result = Normalize(Normal((-0.5*width[1]*z[i-1+nx*j]), (0.5*width[0]*z[i+nx*(j+1)]), -1*width[0]*width[1]));
	//	}
	//		else if ((i == 0) && (j == ny-1)){
	//			result = Normalize(Normal((0.5*width[1]*z[i+1+nx*j]), (-0.5*width[0]*z[i+nx*(j-1)]), -1*width[0]*width[1]));
	//		}
	//		else if ((i == nx-1) && (j == ny-1)){
	//			result = Normalize(Normal((-0.5*width[1]*z[i-1+nx*j]), (-0.5*width[0]*z[i+nx*(j-1)]), -1*width[0]*width[1]));
	//		}
	//		else if ((i > 0) && (i < nx -1) && (j == 0)){
	//			result = Normalize(Normal((0.5*width[1]*(z[i+1+nx*j]-z[i-1+nx*j])), (0.5*width[0]*z[i+nx*(j+1)]), -1*width[0]*width[1]));
	//		}
	//		else if ((i > 0) && (i < nx -1) && (j == ny - 1)){
	//			result = Normalize(Normal((0.5*width[1]*(z[i+1+nx*j]-z[i-1+nx*j])), (-0.5*width[0]*z[i+nx*(j-1)]), -1*width[0]*width[1]));
	//		}
	//		else if ((j > 0) && (j < ny -1) && (i == 0)){
	//			result = Normalize(Normal((0.5*width[1]*z[i+1+nx*j]), (0.5*width[0]*(z[i+nx*(j+1)] - z[i+nx*(j-1)])), -1*width[0]*width[1]));
	//		}
	//		else if ((j > 0) && (j < ny -1) && (i == nx -1)){
	//			result = Normalize(Normal((-0.5*width[1]*z[i-1+nx*j]), (0.5*width[0]*(z[i+nx*(j+1)] - z[i+nx*(j-1)])), -1*width[0]*width[1]));
	//		}
	//		else {
	//			result = Normalize(Normal((0.5*width[1]*(z[i+1+nx*j]-z[i-1+nx*j])), (0.5*width[0]*(z[i+nx*(j+1)] - z[i+nx*(j-1)])), -1*width[0]*width[1]));
	//		}
	//		nors[i + nx * j] = result;
	//	}
	//}
}

Heightfield2 *CreateHeightfield2Shape(const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params) {
    int nu = params.FindOneInt("nu", -1);
    int nv = params.FindOneInt("nv", -1);
    int nitems;
    const float *Pz = params.FindFloat("Pz", &nitems);
    Assert(nitems == nu*nv);
    Assert(nu != -1 && nv != -1 && Pz != NULL);
    return new Heightfield2(o2w, w2o, reverseOrientation, nu, nv, Pz);
}
