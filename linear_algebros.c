#include "linear_algebros.h"
#include <math.h>

/*-------- Conversions        ----------*/

float linalgDeg2Rad(float angle) {
	return angle * pi / 180.0f;
}

/*-------- Vec3 Operations    ----------*/

vec3 linalgMakeVec3(float x, float y, float z) {

	vec3 result;

	result.data[0] = x;
	result.data[1] = y;
	result.data[2] = z;
	result.data[3] = 0;

	return result;
}

float linalgDotVec3(vec3 a, vec3 b) {
	return a.data[0] * b.data[0] + a.data[1] * b.data[1] + a.data[2] * b.data[2];
}

vec3 linalgCross(vec3 a, vec3 b) {
	vec3 result;

	result.data[0] = a.data[1] * b.data[2] - a.data[2] * b.data[1];
	result.data[1] = a.data[2] * b.data[0] - a.data[0] * b.data[2];
	result.data[2] = a.data[0] * b.data[1] - a.data[1] * b.data[0];
	result.data[3] = 0;

	return result;
}

vec3 linalgNormalizeVec3(vec3 a) {

	float invMagnitude = 1.0f / sqrtf(a.data[0] * a.data[0] + a.data[1] * a.data[1] + a.data[2] * a.data[2]);

	vec3 result;

	result.data[0] = a.data[0] * invMagnitude;
	result.data[1] = a.data[1] * invMagnitude;
	result.data[2] = a.data[2] * invMagnitude;
	result.data[3] = 0.0f;

	return result;
}

vec3 linalgSubVec3(vec3 a, vec3 b) {
	
	vec3 result;

	result.vector = _mm_sub_ps(a.vector, b.vector);

	return result;
}

vec3 linalgAddVec3(vec3 a, vec3 b) {

	vec3 result;

	result.vector = _mm_add_ps(a.vector, b.vector);

	return result;
}

vec3 linalgMulVec3(vec3 a, float scalar) {

	vec3 result;

	result.vector = _mm_mul_ps(a.vector, _mm_set1_ps(scalar));

	return result;
}

/*-------- Matrix4 Operations ----------*/

mat4 linalgMakeIdentity4() {
	mat4 result;

	result.column[0] = _mm_setr_ps(1, 0, 0, 0);
	result.column[1] = _mm_setr_ps(0, 1, 0, 0);
	result.column[2] = _mm_setr_ps(0, 0, 1, 0);
	result.column[3] = _mm_setr_ps(0, 0, 0, 1);

	return result;
}

mat4 linalgMakePerspectiveProjection(
	float fovy, float aspect, float near, float far) {

	float yMax = near * tanf(linalgDeg2Rad(fovy / 2));
	float xMax = yMax * aspect;

	/*
	
		The matrix is:

		[E 0  A 0]
		[0 F  B 0]
		[0 0  C D]
		[0 0 -1 0]

		Given by:

		float left{ -xMax };
		float right{ -xMax };
		float top{ -yMax };
		float bottom{ yMax };

		float A{ (right + left) / (right - left) };
		float B{ (top + bottom) / (top - bottom) };
		float C{ -(far + near) / (far - near) };
		float D{ -2.0f * far * near / (far - near) };
		float E{ 2.0f * near / (right - left) };
		float F{ 2.0f * near / (top - bottom) };

		(In practice this simplifies out quite a bit though.)
	*/
	float C = -(far + near) / (far - near);
	float D = -2.0f * far * near / (far - near);
	float E = near / xMax;
	float F = near / yMax;

	mat4 result;

	result.column[0] = _mm_setr_ps(E, 0, 0,  0);
	result.column[1] = _mm_setr_ps(0, F, 0,  0);
	result.column[2] = _mm_setr_ps(0, 0, C, -1);
	result.column[3] = _mm_setr_ps(0, 0, D,  0);

	return result;
}

mat4 linalgMakeLookAt(vec3 eye, vec3 target, vec3 up) {

	vec3 forwards = linalgNormalizeVec3(linalgSubVec3(target, eye));
	vec3 right = linalgNormalizeVec3(linalgCross(forwards, up));
	up = linalgNormalizeVec3(linalgCross(right, forwards));

	forwards = linalgMulVec3(forwards, -1);

	mat4 result;

	result.column[0] = _mm_setr_ps(right.data[0], up.data[0], forwards.data[0], 0);
	result.column[1] = _mm_setr_ps(right.data[1], up.data[1], forwards.data[1], 0);
	result.column[2] = _mm_setr_ps(right.data[2], up.data[2], forwards.data[2], 0);
	result.column[3] = _mm_setr_ps(-linalgDotVec3(right, eye), -linalgDotVec3(up, eye), -linalgDotVec3(forwards, eye), 1);

	return result;
}

mat4 linalgMakeTranslation(vec3 translation) {

	mat4 result;

	result.column[0] = _mm_setr_ps(1, 0, 0, 0);
	result.column[1] = _mm_setr_ps(0, 1, 0, 0);
	result.column[2] = _mm_setr_ps(0, 0, 1, 0);
	result.column[3] = _mm_setr_ps(translation.data[0], translation.data[1], translation.data[2], 1);

	return result;
}

mat4 linalgMakeXRotation(float angle) {

	angle = linalgDeg2Rad(angle);
	float cT = cosf(angle);
	float sT = sinf(angle);

	mat4 result;

	result.column[0] = _mm_setr_ps(1,  0,   0, 0);
	result.column[1] = _mm_setr_ps(0, cT, -sT, 0);
	result.column[2] = _mm_setr_ps(0, sT,  cT, 0);
	result.column[3] = _mm_setr_ps(0,  0,   0, 1);

	return result;
}

mat4 linalgMakeYRotation(float angle) {

	angle = linalgDeg2Rad(angle);
	float cT = cosf(angle);
	float sT = sinf(angle);

	mat4 result;

	result.column[0] = _mm_setr_ps( cT, 0, sT, 0);
	result.column[1] = _mm_setr_ps(  0, 1,  0, 0);
	result.column[2] = _mm_setr_ps(-sT, 0, cT, 0);
	result.column[3] = _mm_setr_ps(  0, 0,  0, 1);

	return result;
}

mat4 linalgMakeZRotation(float angle) {

	angle = linalgDeg2Rad(angle);
	float cT = cosf(angle);
	float sT = sinf(angle);

	mat4 result;

	result.column[0] = _mm_setr_ps(cT, -sT, 0, 0);
	result.column[1] = _mm_setr_ps(sT,  cT, 0, 0);
	result.column[2] = _mm_setr_ps( 0,   0, 1, 0);
	result.column[3] = _mm_setr_ps( 0,   0, 0, 1);

	return result;
}

vec4 linalgMulMat4Vec4(mat4 m, vec4 v) {

	vec4 result;

	result.vector = _mm_fmadd_ps(_mm_set1_ps(v.data[0]), m.column[0],
					_mm_fmadd_ps(_mm_set1_ps(v.data[1]), m.column[1],
					_mm_fmadd_ps(_mm_set1_ps(v.data[2]), m.column[2],
					_mm_mul_ps(_mm_set1_ps(v.data[3]), m.column[3])
					)
				)
	);

	return result;
}

mat4 linalgMulMat4Mat4(mat4 m1, mat4 m2) {

	mat4 result;

	result.column_vector[0] = linalgMulMat4Vec4(m2, m1.column_vector[0]);
	result.column_vector[1] = linalgMulMat4Vec4(m2, m1.column_vector[1]);
	result.column_vector[2] = linalgMulMat4Vec4(m2, m1.column_vector[2]);
	result.column_vector[3] = linalgMulMat4Vec4(m2, m1.column_vector[3]);

	return result;
}