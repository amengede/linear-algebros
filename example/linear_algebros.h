#include <intrin.h>

#ifndef __linear_algebros_h__
#define __linear_algebros_h__

typedef struct vec3 {
	union {
		__m128 vector;
		float data[4];
	};
} vec3;

typedef struct vec4 {
	union {
		__m128 vector;
		float data[4];
	};
} vec4;

typedef struct mat4 {
	union {
		__m128 column[4];
		vec4 column_vector[4];
		float data[16];
	};
} mat4;

/*-------- Conversions        ----------*/

#define pi 3.14159265359f

float linalgDeg2Rad(float angle);

/*-------- Vec3 Operations    ----------*/

vec3 linalgMakeVec3(float x, float y, float z);

float linalgDotVec3(vec3 a, vec3 b);

vec3 linalgCross(vec3 a, vec3 b);

vec3 linalgNormalizeVec3(vec3 a);

vec3 linalgSubVec3(vec3 a, vec3 b);

vec3 linalgAddVec3(vec3 a, vec3 b);

vec3 linalgMulVec3(vec3 a, float scalar);

/*-------- Matrix4 Operations ----------*/

mat4 linalgMakeIdentity4();

mat4 linalgMakePerspectiveProjection(float fovy, float aspect, float near, float far);

mat4 linalgMakeLookAt(vec3 eye, vec3 target, vec3 up);

mat4 linalgMakeTranslation(vec3 translation);

mat4 linalgMakeXRotation(float angle);

mat4 linalgMakeYRotation(float angle);

mat4 linalgMakeZRotation(float angle);

vec4 linalgMulMat4Vec4(mat4 m, vec4 v);

mat4 linalgMulMat4Mat4(mat4 m1, mat4 m2);

#endif // !__linear_algebros_h__
