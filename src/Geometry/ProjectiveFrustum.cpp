/* Copyright Jukka Jylänki

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License. */

/** @file ProjectiveFrustum.cpp
	@author Jukka Jylänki
	@brief Implementation for the ProjectiveFrustum geometry object. */
#include "ProjectiveFrustum.h"
#include "AABB.h"
#include "Circle.h"
#include "../Math/MathFunc.h"
#include "Plane.h"
#include "Line.h"
#include "OBB.h"
#include "Polyhedron.h"
#include "Polygon.h"
#include "Ray.h"
#include "Sphere.h"
#include "Capsule.h"
#include "Triangle.h"
#include "LineSegment.h"
#include "PBVolume.h"
#include "../Math/float2.h"
#include "../Math/float3x3.h"
#include "../Math/float3x4.h"
#include "../Math/float4x4.h"
#include "../Math/float4.h"
#include "../Math/Quat.h"
#include "../Algorithm/Random/LCG.h"
#include "../Algorithm/GJK.h"

#ifdef MATH_ENABLE_STL_SUPPORT
#include <iostream>
#endif

#if defined(MATH_TINYXML_INTEROP) && defined(MATH_CONTAINERLIB_SUPPORT)
#include "Container/UString.h"
#endif

#if defined(MATH_SIMD) && defined(MATH_AUTOMATIC_SSE)
#include "../Math/float4_sse.h"
#include "../Math/float4_neon.h"
#endif

MATH_BEGIN_NAMESPACE

Plane ProjectiveFrustum::NearPlane() const
{
	return GetPlane(0);
}

Plane ProjectiveFrustum::FarPlane() const
{
	return GetPlane(1);
}

Plane ProjectiveFrustum::LeftPlane() const
{
	return GetPlane(2);
}

Plane ProjectiveFrustum::RightPlane() const
{
	return GetPlane(3);
}

Plane ProjectiveFrustum::TopPlane() const
{
	return GetPlane(4);
}

Plane ProjectiveFrustum::BottomPlane() const
{
	return GetPlane(5);
}

Plane ProjectiveFrustum::GetPlane(int faceIndex) const
{
	assume(0 <= faceIndex && faceIndex <= 5);
	float4 r = inverse.Row(faceIndex>>1);
	if(faceIndex & 1) r = -r;
	r -= inverse.Row(3);
	return r;
}

vec ProjectiveFrustum::CenterPoint() const
{
	float4 p = matrix.Col(3);
	return p.Float3Part().Div(p[3]);
}

/// Generates one of the eight corner points of this Frustum.
/** @param cornerIndex The index of the corner point to generate, in the range [0, 7].
	The points are returned in the order 0: ---, 1: --+, 2: -+-, 3: -++, 4: +--, 5: +-+, 6: ++-, 7: +++.
	(corresponding the XYZ axis directions). */
vec ProjectiveFrustum::CornerPoint(int cornerIndex) const
{
	float4 p = matrix.Col(3);
	p += (cornerIndex & 4) ? matrix.Col(0) : - matrix.Col(0);
	p += (cornerIndex & 2) ? matrix.Col(1) : - matrix.Col(1);
	p += (cornerIndex & 1) ? matrix.Col(2) : - matrix.Col(2);
	return p.Float3Part().Div(p[3]);
}

float ProjectiveFrustum::Volume() const
{
	float4 r3   = matrix.Row(3);
	float4 r3_2 = r3.Mul(r3);
	float4 r3_4 = r3_2.Mul(r3_2);
	float sub_2 = r3_2[3]-r3_2[0]-r3_2[1]-r3_2[2];
	float sub_4 = r3_4[3]-r3_4[0]-r3_4[1]-r3_4[2];
/*
std::cout << "det="<<std::abs(matrix.Determinant4()) << std::endl;
std::cout << "r3_1="<<r3 << std::endl;
std::cout << "r3_2="<<r3_2 << std::endl;
std::cout << "r3_4="<<r3_4 << std::endl;
std::cout << "num="<<2*sub_4+sub_2*sub_2 << std::endl;
std::cout << "den="<<(r3[3]-r3[0]-r3[1]-r3[2]) *
		(r3[3]+r3[0]-r3[1]-r3[2]) *
		(r3[3]-r3[0]+r3[1]-r3[2]) *
		(r3[3]+r3[0]+r3[1]-r3[2]) *
		(r3[3]-r3[0]-r3[1]+r3[2]) *
		(r3[3]+r3[0]-r3[1]+r3[2]) *
		(r3[3]-r3[0]+r3[1]+r3[2]) *
		(r3[3]+r3[0]+r3[1]+r3[2]) << std::endl;
*/
	return 8.0*std::abs(matrix.Determinant4())*(2*sub_4+sub_2*sub_2) / (3.0 *
		(r3[3]-r3[0]-r3[1]-r3[2]) *
		(r3[3]+r3[0]-r3[1]-r3[2]) *
		(r3[3]-r3[0]+r3[1]-r3[2]) *
		(r3[3]+r3[0]+r3[1]-r3[2]) *
		(r3[3]-r3[0]-r3[1]+r3[2]) *
		(r3[3]+r3[0]-r3[1]+r3[2]) *
		(r3[3]-r3[0]+r3[1]+r3[2]) *
		(r3[3]+r3[0]+r3[1]+r3[2]) );
}

bool ProjectiveFrustum::Contains(const vec &point) const
{
	float4 p = inverse*float4(point,1);
	vec q = p.Float3Part().Div(p[3]);
	return q.x>=-1 && q.x<=1 && q.y>=-1 && q.y<=1 && q.z>=-1 && q.z<=1;
}

float ProjectiveFrustum::AlgebraicDistance(const vec &point) const
{
	float4 p = inverse*float4(point,1);
	vec q = p.Float3Part().Div(p.w);
	return q.Abs().MaxElement()-1;
}

void ProjectiveFrustum::Transform(const float3x3 &transform)
{
	matrix = matrix * transform;
}

void ProjectiveFrustum::Transform(const float3x4 &transform)
{
	matrix = matrix * transform;
}

void ProjectiveFrustum::Transform(const float4x4 &transform)
{
	matrix = matrix * transform;
}

void ProjectiveFrustum::Transform(const Quat &transform)
{
	matrix = matrix * transform;
}

void ProjectiveFrustum::GetPlanes(Plane *outArray) const
{
	assume(outArray);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (!outArray)
		return;
#endif
	for(int i = 0; i < 6; ++i)
		outArray[i] = GetPlane(i);
}

PBVolume<6> ProjectiveFrustum::ToPBVolume() const
{
	PBVolume<6> frustumVolume;
	frustumVolume.p[0] = NearPlane();
	frustumVolume.p[1] = LeftPlane();
	frustumVolume.p[2] = RightPlane();
	frustumVolume.p[3] = TopPlane();
	frustumVolume.p[4] = BottomPlane();
	frustumVolume.p[5] = FarPlane();

	return frustumVolume;
}

ProjectiveFrustum ProjectiveFrustum::OptimalEnclosingProjectiveFrustum(const vec *pointArray, int numPoints)
{
	// Precomputation: Generate the convex hull of the input point set. This is because
	// we need vertex-edge-face connectivity information about the convex hull shape, and
	// this also allows discarding all points in the interior of the input hull, which
	// are irrelevant.
	Polyhedron convexHull = Polyhedron::ConvexHull(pointArray, numPoints);
	if (!pointArray || convexHull.v.size() == 0)
	{
		ProjectiveFrustum minFrustum(float4x4::nan);
		return minFrustum;
	}
	return OptimalEnclosingProjectiveFrustum(convexHull);
}


ProjectiveFrustum ProjectiveFrustum::OptimalEnclosingProjectiveFrustum(const Polyhedron &convexHull)
{
	ProjectiveFrustum minFrustum(float4x4::nan);
	float minVolume = FLOAT_INF;

	// Handle degenerate planar cases up front.
	if (convexHull.v.size() <= 3 || convexHull.f.size() <= 1)
	{
		// TODO
		LOGW("Convex hull is degenerate and has only %d vertices/%d faces!", (int)convexHull.v.size(), (int)convexHull.f.size());
		return minFrustum;
	}

	std::vector<float4> facePlanes;
	for(int i = 0; i < convexHull.NumFaces(); ++i) // O(F)
	{
		if (convexHull.f[i].v.size() < 3)
		{
			LOGE("Input convex hull contains a degenerate face %d with only %d vertices! Cannot process this!",
				i, (int)convexHull.f[i].v.size());
			return minFrustum;
		}
		Plane plane = convexHull.FacePlane(i);
		facePlanes.push_back(float4(plane.normal,-plane.d));

		assert((EqualAbs(0,facePlanes.back().Dot(float4(convexHull.Vertex(convexHull.f[i].v[0]),1)))));
		assert((EqualAbs(0,facePlanes.back().Dot(float4(convexHull.Vertex(convexHull.f[i].v[1]),1)))));
		assert((EqualAbs(0,facePlanes.back().Dot(float4(convexHull.Vertex(convexHull.f[i].v[2]),1)))));
	}

	std::vector< std::pair< int, int > > edges = convexHull.EdgeIndices();
	float4x4 n;
	float4x4 nv(float4x4::nan);
	float4x4 inverse;
	const int F = convexHull.NumFaces();
	const int E = convexHull.NumEdges();
	const int V = convexHull.NumVertices();
	bool ok = true;

	std::cout << std::endl;
	std::cout << F << " faces" << std::endl;
	std::cout << E << " edges" << std::endl;
	std::cout << V << " vertices" << std::endl;


	for(int f0 = 0; f0 < F; ++f0) // O(F), face on X=-1
	{
		std::cout << f0 << "..." << std::flush;
		n.Row(0) = facePlanes[f0];
		for(int f1 = f0+1; f1 < F; ++f1) // O(F), face on X=1
		{
			ok = true;
			for(int i=0;i<9 && ok;++i)
				ok = (convexHull.f[f0].v[i%3] != convexHull.f[f1].v[i/3]);
			if(!ok) continue;
			n.Row(1) = facePlanes[f1];
			for(int f2 = 0; f2 < F; ++f2) // O(F), face on Y=-1
			{
				if(f2==f0 || f2==f1) continue;
				n.Row(2) = facePlanes[f2];
				for(int f3 = 0; f3 < F; ++f3) // O(F), face on Z=-1
				{
					if(f3==f0 || f3==f1 || f3==f2) continue;
					n.Row(3) = facePlanes[f3];
					for(int e = 0; e < E; ++e) // O(E), edge on Y=1 : (m3-m1).v12 = 0
					{
						ok = true;
						for(int i=0;i<3 && ok;++i)
							ok = (edges[e].first != convexHull.f[f2].v[i]) && (edges[e].second != convexHull.f[f2].v[i]) ;
						if(!ok) continue;

						float4 v1(convexHull.Vertex(edges[e].first),1);
						float4 v2(convexHull.Vertex(edges[e].second),1);
						float4 nv1 = n*v1;
						float4 nv2 = n*v2;
						nv.Row(1) = nv1.Mul(float4(1,1,-1,0));
						nv.Row(2) = nv2.Mul(float4(1,1,-1,0));
						
						for(int v = 0; v < V; ++v) // O(V), edge on Z=1 : (m3-m2).v0 = 0
						{
							ok = true;
							for(int i=0;i<3 && ok;++i)
								ok = (v != convexHull.f[f3].v[i]);
							if(!ok) continue;
/*
							std::cout << "=X-= " << convexHull.f[f0].v[0] << ", " << convexHull.f[f0].v[1] << ", " << convexHull.f[f0].v[2] << ", " ;
							std::cout << "=X+= " << convexHull.f[f1].v[0] << ", " << convexHull.f[f1].v[1] << ", " << convexHull.f[f1].v[2] << ", " ;
							std::cout << "=Y-= " << convexHull.f[f2].v[0] << ", " << convexHull.f[f2].v[1] << ", " << convexHull.f[f2].v[2] << ", " ;
							std::cout << "=Z-= " << convexHull.f[f3].v[0] << ", " << convexHull.f[f3].v[1] << ", " << convexHull.f[f3].v[2] << ", " ;
							std::cout << "=Y+= " << edges[e].first << ", " << edges[e].second << ", " ;
							std::cout << "=Z+= " << v << std::endl ;
*/
							float4 v0(convexHull.Vertex(v),1);
							float4 nv0 = n*v0;
							nv.Row(0) = nv0.Mul(float4(1,1,0,-1));

							float4 lambda(-nv.Minor(3,0),nv.Minor(3,1),-nv.Minor(3,2),nv.Minor(3,3));
							lambda.Normalize();
/*
							if(!EqualAbs(0,lambda.Dot(nv.Row(0)))) std::cerr << "lambda.Dot(nv.Row(0)) : " << lambda.Dot(nv.Row(0)) << lambda << nv.Row(0) << std::endl;
							if(!EqualAbs(0,lambda.Dot(nv.Row(1)))) std::cerr << "lambda.Dot(nv.Row(1)) : " << lambda.Dot(nv.Row(1)) << lambda << nv.Row(1)<< std::endl;
							if(!EqualAbs(0,lambda.Dot(nv.Row(2)))) std::cerr << "lambda.Dot(nv.Row(2)) : " << lambda.Dot(nv.Row(2)) << lambda << nv.Row(2)<< std::endl;
*/
							inverse.Row(3) = lambda.x*n.Row(0)+lambda.y*n.Row(1);
							inverse.Row(0) = lambda.x*n.Row(0)-lambda.y*n.Row(1);
							inverse.Row(1) = 2*lambda.z*n.Row(2)-inverse.Row(3);
							inverse.Row(2) = 2*lambda.w*n.Row(3)-inverse.Row(3);
/*
							std::cout << "X-: " << inverse*float4(convexHull.Vertex(convexHull.f[f0].v[0]),1) << std::endl;
							std::cout << "X-: " << inverse*float4(convexHull.Vertex(convexHull.f[f0].v[1]),1) << std::endl;
							std::cout << "X-: " << inverse*float4(convexHull.Vertex(convexHull.f[f0].v[2]),1) << std::endl;
							std::cout << "X+: " << inverse*float4(convexHull.Vertex(convexHull.f[f1].v[0]),1) << std::endl;
							std::cout << "X+: " << inverse*float4(convexHull.Vertex(convexHull.f[f1].v[1]),1) << std::endl;
							std::cout << "X+: " << inverse*float4(convexHull.Vertex(convexHull.f[f1].v[2]),1) << std::endl;
							std::cout << "Y-: " << inverse*float4(convexHull.Vertex(convexHull.f[f2].v[0]),1) << std::endl;
							std::cout << "Y-: " << inverse*float4(convexHull.Vertex(convexHull.f[f2].v[1]),1) << std::endl;
							std::cout << "Y-: " << inverse*float4(convexHull.Vertex(convexHull.f[f2].v[2]),1) << std::endl;
							std::cout << "Z-: " << inverse*float4(convexHull.Vertex(convexHull.f[f3].v[0]),1) << std::endl;
							std::cout << "Z-: " << inverse*float4(convexHull.Vertex(convexHull.f[f3].v[1]),1) << std::endl;
							std::cout << "Z-: " << inverse*float4(convexHull.Vertex(convexHull.f[f3].v[2]),1) << std::endl;
							std::cout << "Y+: " << inverse*v1 << std::endl;
							std::cout << "Y+: " << inverse*v2 << std::endl;
							std::cout << "Z+: " << inverse*v0 << std::endl;
*/
/*
							std::cout << "matrix: " << matrix << std::endl;
							std::cout << "lambda: " << lambda << std::endl;
							std::cout << "n: " << n << std::endl;
*/
/*
							std::cout << "nv0: " << nv0 << std::endl;
							std::cout << "nv1: " << nv1 << std::endl;
							std::cout << "nv2: " << nv2 << std::endl;
							std::cout << "mv0: " << matrix*v0 << std::endl;
							std::cout << "mv1: " << matrix*v1 << std::endl;
							std::cout << "mv2: " << matrix*v2 << std::endl;
							assert(EqualAbs(0,n.Row(0).Dot(float4(convexHull.Vertex(convexHull.f[f0].v[0]),1))));
							assert(EqualAbs(0,n.Row(1).Dot(float4(convexHull.Vertex(convexHull.f[f1].v[0]),1))));
							assert(EqualAbs(0,n.Row(2).Dot(float4(convexHull.Vertex(convexHull.f[f2].v[0]),1))));
							assert(EqualAbs(0,n.Row(3).Dot(float4(convexHull.Vertex(convexHull.f[f3].v[0]),1))));

							std::cout << "f0: "<< (matrix.Row(3)-matrix.Row(0)).Dot(float4(convexHull.Vertex(convexHull.f[f0].v[0]),1)) << std::endl;
							std::cout << "f0: "<< (matrix.Row(3)-matrix.Row(0)).Dot(float4(convexHull.Vertex(convexHull.f[f0].v[1]),1)) << std::endl;
							std::cout << "f0: "<< (matrix.Row(3)-matrix.Row(0)).Dot(float4(convexHull.Vertex(convexHull.f[f0].v[2]),1)) << std::endl;
							std::cout << (matrix.Row(3)-matrix.Row(2)).Dot(v0) << std::endl;
							std::cout << (matrix.Row(3)-matrix.Row(1)).Dot(v1) << std::endl;
							std::cout << (matrix.Row(3)-matrix.Row(1)).Dot(v2) << std::endl;
*/
							ProjectiveFrustum frustum(inverse.Inverted(),inverse);
							//for(int i=0; i<8; ++i) std::cout << frustum.CornerPoint(i)<< std::endl;

							// check inside
							float eps = 1e-3;
							if(frustum.AlgebraicDistance(v0.xyz()/v0.w)>eps) continue;
							if(frustum.AlgebraicDistance(v1.xyz()/v1.w)>eps) continue;
							if(frustum.AlgebraicDistance(v2.xyz()/v2.w)>eps) continue;

							float volume = frustum.Volume();
							if(volume < minVolume)
							{
								minVolume=volume;
								minFrustum=frustum;
							}
						}
					}
				}
			}
		}
	}

	return  minFrustum;
}

ProjectiveFrustum operator *(const float3x3 &transform, const ProjectiveFrustum &frustum)
{
	ProjectiveFrustum f(frustum);
	f.Transform(transform);
	return f;
}

ProjectiveFrustum operator *(const float3x4 &transform, const ProjectiveFrustum &frustum)
{
	ProjectiveFrustum f(frustum);
	f.Transform(transform);
	return f;
}

ProjectiveFrustum operator *(const float4x4 &transform, const ProjectiveFrustum &frustum)
{
	ProjectiveFrustum f(frustum);
	f.Transform(transform);
	return f;
}

ProjectiveFrustum operator *(const Quat &transform, const ProjectiveFrustum &frustum)
{
	ProjectiveFrustum f(frustum);
	f.Transform(transform);
	return f;
}

MATH_END_NAMESPACE
