#include <stdio.h>
#include <stdlib.h>

#include "../src/MathGeoLib.h"
#include "../src/Math/myassert.h"
#include "TestRunner.h"
#include "TestData.h"
#include "../src/Geometry/PBVolume.h"
#include "ObjectGenerators.h"

MATH_IGNORE_UNUSED_VARS_WARNING

using namespace TestData;


float4x4 testMatrix[] = {
	float4x4( 2, 0, 0, 10,	0, 1, 0, 10,	0, 0, 3, 5,	0, 0, 0, 1 ),
	float4x4( 2, 0, 0, 10,	0, 1, 0, 10,	0, 0, 3, 5,	0, 0, 0, 1 )
};

float volume[] = { 48, 48 };

UNIQUE_TEST(ProjectiveFrustum_Planes)
{
	for(int m=0; m*(sizeof(*testMatrix))<sizeof(testMatrix); ++m)
	{
		ProjectiveFrustum f(testMatrix[m]);
		for(int i = 0; i < 9; ++i)
		{
			vec pt;
			if (i == 8)
				pt = f.CenterPoint();
			else
				pt = f.CornerPoint(i);

			assert(f.NearPlane().SignedDistance(pt) < 1e-3f);
			assert(f.FarPlane().SignedDistance(pt) < 1e-3f);
			assert(f.LeftPlane().SignedDistance(pt) < 1e-3f);
			assert(f.RightPlane().SignedDistance(pt) < 1e-3f);
			assert(f.TopPlane().SignedDistance(pt) < 1e-3f);
			assert(f.BottomPlane().SignedDistance(pt) < 1e-3f);
			assert(f.Contains(pt));
		}
	}
}

UNIQUE_TEST(ProjectiveFrustum_Corners)
{
	for(int m=0; m*(sizeof(*testMatrix))<sizeof(testMatrix); ++m)
	{
		ProjectiveFrustum f(testMatrix[m]);
	}
}

UNIQUE_TEST(ProjectiveFrustum_Volume)
{
	for(int m=0; m*(sizeof(*testMatrix))<sizeof(testMatrix); ++m)
	{
		ProjectiveFrustum f(testMatrix[m]);
		assert(EqualRel(f.Volume(), volume[m]));
	}
}

RANDOMIZED_TEST(ProjectiveFrustum_OptimalEnclosingProjectiveFrustum)
{
	const int N = 20;
	vec points[N];
	for(int i = 0; i<N; ++i)
	{
		points[i] = vec::RandomSphere(rng,vec(0,0,0),10);
	}

	ProjectiveFrustum f = ProjectiveFrustum::OptimalEnclosingProjectiveFrustum(points,N);
	for(int i = 0; i<N; ++i)
	{
		assert(f.Contains(points[i]));
	}
}

RANDOMIZED_TEST(ProjectiveFrustum_Volume_Random)
{
	float4x4 matrix = float4x4::RandomGeneral(rng,-10,10);
	matrix.At(3,3) = matrix.Row(3).Abs().SumOfElements(); // ensure m33 > |m30|+|m31|+|m32|
	ProjectiveFrustum f(matrix);

	PBVolume<6> pbvol = f.ToPBVolume();
	Polyhedron ph = pbvol.ToPolyhedron();
	assert(EqualRel(f.Volume(), ph.Volume(),1e-2f));
	
}

