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

/** @file ProjectiveFrustum.h
	@author Mathieu Brédif
	@brief The Projective Frustum geometry object. */
#pragma once

#include "../MathGeoLibFwd.h"
#include "../Math/float2.h"
#include "../Math/float3.h"
#include "../Math/float3x4.h"
#include "../Math/float4x4.h"
#include "Ray.h"

#ifdef MATH_TINYXML_INTEROP
#include "Config/tinyxml/tinyxml.h"
#endif

MATH_BEGIN_NAMESPACE

class ProjectiveFrustum
{
private:
	float4x4 matrix;
	float4x4 inverse;

public:
	ProjectiveFrustum(const float4x4 &matrix_) : matrix(matrix_), inverse(matrix_.Inverted()) { };
	ProjectiveFrustum(const float4x4 &matrix_,const float4x4 &inverse_) : matrix(matrix_), inverse(inverse_) { };

	/// Returns the number of line segment edges that this Frustum is made up of, which is always 12.
	/** This function is used in template-based algorithms to provide an unified API for iterating over the features of a Polyhedron. */
	int NumEdges() const { return 12; }

	/// Computes the plane equation of the near plane of this Frustum.
	/** The normal vector of the returned plane points outwards from the volume inside the frustum, i.e. towards the eye point
		(towards -front). This means the negative half-space of the Frustum is the space inside the Frustum.
		@see front, FarPlane(), LeftPlane(), RightPlane(), TopPlane(), BottomPlane(), GetPlane(), GetPlanes(). */
	Plane NearPlane() const;

	/// Computes the plane equation of the far plane of this Frustum. [similarOverload: NearPlane]
	/** The normal vector of the returned plane points outwards from the volume inside the frustum, i.e. away from the eye point.
		(towards front). This means the negative half-space of the Frustum is the space inside the Frustum.
		@see front, FarPlane(), LeftPlane(), RightPlane(), TopPlane(), BottomPlane(), GetPlane(), GetPlanes(). */
	Plane FarPlane() const;

	/// Returns the plane equation of the specified side of this Frustum.
	/** The normal vector of the returned plane points outwards from the volume inside the frustum.
		This means the negative half-space of the Frustum is the space inside the Frustum.
		[indexTitle: Left/Right/Top/BottomPlane]
		@see NearPlane(), FarPlane(), GetPlane(), GetPlanes(). */
	Plane LeftPlane() const;
	Plane RightPlane() const; ///< [similarOverload: LeftPlane] [hideIndex]
	Plane TopPlane() const; ///< [similarOverload: LeftPlane] [hideIndex]
	Plane BottomPlane() const; ///< [similarOverload: LeftPlane] [hideIndex]

	/// Returns the specified plane of this frustum.
	/** The normal vector of the returned plane points outwards from the volume inside the frustum.
		@param faceIndex A number in the range [0,5], which returns the plane at the selected index from
			the array { near, far, left, right, top, bottom }.
		@see GetPlanes(), NearPlane(), FarPlane(), LeftPlane(), RightPlane(), TopPlane(), BottomPlane(). */
	Plane GetPlane(int faceIndex) const;

	/// Returns all six planes of this Frustum.
	/** The planes will be output in the order { near, far, left, right, top, bottom }.
		@param outArray [out] A pointer to an array of at least 6 elements. This pointer will receive the planes of this Frustum.
			This pointer may not be null.
		@see GetPlane(), NearPlane(), FarPlane(), LeftPlane(), RightPlane(), TopPlane(), BottomPlane(). */
	void GetPlanes(Plane *outArray) const;

	vec CenterPoint() const;

	/// Generates one of the eight corner points of this Frustum.
	/** @param cornerIndex The index of the corner point to generate, in the range [0, 7].
		 The points are returned in the order 0: ---, 1: --+, 2: -+-, 3: -++, 4: +--, 5: +-+, 6: ++-, 7: +++.
		 (corresponding the XYZ axis directions). */
	vec CornerPoint(int cornerIndex) const;

	/// Returns all eight corner points of this array.
	/** @param outPointArray [out] A pointer to an array of at least 8 elements. This pointer will receive the corner vertices
			of this Frustum. This pointer may not be null. */
	void GetCornerPoints(vec *outPointArray) const;

	/// returns the forward matrix defining this Frustum.
	/** @return A 4x4 matrix. */
	float4x4 Matrix() const { return matrix; }
	/// returns the backward matrix defining this Frustum.
	/** @return A 4x4 matrix. */
	float4x4 Inverse() const { return inverse; }

	/// Computes the volume of this Frustum.
	float Volume() const;

	bool Contains(const vec &point) const;

	/// >0 outside, 0 on boundary, <0 inside
	float AlgebraicDistance(const vec &point) const;

	/// Applies a transformation to this Frustum.
	/** @param transform The transformation to apply to this Frustum. This transformation must be
		affine, and must contain an orthogonal set of column vectors (may not contain shear or projection).
		The transformation can only contain uniform scale, and may not contain mirroring.
		@see Translate(), Scale(), classes float3x3, float3x4, float4x4, Quat. */
	void Transform(const float3x3 &transform);
	void Transform(const float3x4 &transform);
	void Transform(const float4x4 &transform);
	void Transform(const Quat &transform);

	/// Converts this Frustum to a PBVolume.
	/** This function returns a plane-bounded volume representation of this Frustum. The conversion is exact, meaning that the
		returned PBVolume<6> represents exactly the same set of points that this Frustum does.
		@see ToPolyhedron(). */
	PBVolume<6> ToPBVolume() const;

	static ProjectiveFrustum OptimalEnclosingProjectiveFrustum(const vec *pointArray, int numPoints);
	static ProjectiveFrustum OptimalEnclosingProjectiveFrustum(const Polyhedron &convexPolyhedron);
};

ProjectiveFrustum operator *(const float3x3 &transform, const ProjectiveFrustum &frustum);
ProjectiveFrustum operator *(const float3x4 &transform, const ProjectiveFrustum &frustum);
ProjectiveFrustum operator *(const float4x4 &transform, const ProjectiveFrustum &frustum);
ProjectiveFrustum operator *(const Quat &transform, const ProjectiveFrustum &frustum);

MATH_END_NAMESPACE
