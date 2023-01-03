#include "linear_algebros.h"
#include <math.h>

/*-------- Conversions        ----------*/

float linalgDeg2Rad(float angle) {
	return angle * pi / 180.0f;
}

float linalgRad2Deg(float angle) {
	return angle * 180.0f / pi;
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

	result.vector = _mm_mul_ps(a.vector, _mm_set1_ps(invMagnitude));

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

float linalgAngleBetweenVectors3(vec3 a, vec3 b) {

	float denominator = sqrtf(linalgDotVec3(a, a) * linalgDotVec3(b, b));

	float dot = linalgDotVec3(a, b);

	return linalgRad2Deg(acosf(dot / denominator));
}

vec3 linalgProject(vec3 incoming, vec3 basis) {

	return linalgMulVec3(basis, linalgDotVec3(incoming, basis) / linalgDotVec3(basis, basis));
}

vec3 linalgReject(vec3 incoming, vec3 basis) {
	
	return linalgSubVec3(incoming, linalgProject(incoming, basis));
}

vec3 linalgReflect(vec3 incident, vec3 normal) {

	vec3 result;
	// reflected = incident − 2(incident.normal)normal
	result.vector = _mm_fmadd_ps(
		normal.vector, 
		_mm_set1_ps(-2.0f * linalgDotVec3(incident, normal)),
		incident.vector
	);

	return result;
}

vec3 linalgLerpVec3(vec3 a, vec3 b, float t) {

	vec3 result;

	result.vector = _mm_fmadd_ps(
		_mm_sub_ps(b.vector, a.vector),
		_mm_set1_ps(t),
		a.vector
	);

	return result;
}

vec3 linalgSlerpVec3(vec3 a, vec3 b, float t) {

	if (t < 0.1f) {
		return linalgLerpVec3(a, b, t);
	}

	float angle = linalgAngleBetweenVectors3(a, b);

	float denominator = sinf(linalgDeg2Rad(angle));

	vec3 result;

	result.vector = _mm_fmadd_ps(
		a.vector,
		_mm_set1_ps(sinf(1 - t) * angle / denominator),
		_mm_mul_ps(b.vector, _mm_set1_ps(sinf(t * angle) / denominator))
	);

	return result;
}

vec3 linalgNlerpVec3(vec3 a, vec3 b, float t) {

	return linalgNormalizeVec3(linalgLerpVec3(a, b, t));
}

bool linalgCloseVec3(vec3 a, vec3 b) {
	
	vec3 displacement = linalgSubVec3(a, b);

	return linalgDotVec3(displacement, displacement) < eps;
}

/*-------- Vector4 Operations ----------*/

float linalgDotVec4(vec4 a, vec4 b) {
	return a.data[0] * b.data[0] + a.data[1] * b.data[1] + a.data[2] * b.data[2] + a.data[3] * b.data[3];
}

vec4 linalgNormalizeVec4(vec4 a) {

	float invMagnitude = 1.0f / sqrtf(
		a.data[0] * a.data[0] + a.data[1] * a.data[1] + a.data[2] * a.data[2] + a.data[3] * a.data[3]
	);

	vec4 result;

	result.vector = _mm_mul_ps(a.vector, _mm_set1_ps(invMagnitude));

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

mat4 linalgAddMat4(mat4 m1, mat4 m2) {

	mat4 m3;
	m3.chunk[0] = _mm256_add_ps(m1.chunk[0], m2.chunk[0]);
	m3.chunk[1] = _mm256_add_ps(m1.chunk[1], m2.chunk[1]);

	return m3;
}

mat4 linalgMulMat4Scalar(mat4 matrix, float scalar) {

	mat4 m3;
	__m256 scale = _mm256_set1_ps(scalar);
	m3.chunk[0] = _mm256_mul_ps(matrix.chunk[0], scale);
	m3.chunk[1] = _mm256_mul_ps(matrix.chunk[1], scale);

	return m3;
}

mat4 linalgLerpMat4(mat4 m1, mat4 m2, float t) {

	mat4 m3;
	__m256 scale = _mm256_set1_ps(t);

	m3.chunk[0] = _mm256_fmadd_ps(
		_mm256_sub_ps(m2.chunk[0], m1.chunk[0]),
		scale,
		m1.chunk[0]
	);

	m3.chunk[1] = _mm256_fmadd_ps(
		_mm256_sub_ps(m2.chunk[1], m1.chunk[1]),
		scale,
		m1.chunk[1]
	);

	return m3;
}

mat4 linalgTranspose(mat4 matrix) {

	mat4 transposed;

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			transposed.data[i + 4 * j] = matrix.data[j + 4 * i];
		}
	}

	return transposed;
}

mat4 linalgTransformInverse(mat4 matrix) {

	//Get the scale factors
	float a = sqrtf(linalgDotVec4(matrix.column_vector[0], matrix.column_vector[0]));
	float b = sqrtf(linalgDotVec4(matrix.column_vector[1], matrix.column_vector[1]));
	float c = sqrtf(linalgDotVec4(matrix.column_vector[2], matrix.column_vector[2]));

	//Get the rotation vectors, apply inverse scaling
	vec4 X = linalgNormalizeVec4(matrix.column_vector[0]);
	X.vector = _mm_mul_ps(X.vector, _mm_set1_ps(1 / a));
	vec4 Y = linalgNormalizeVec4(matrix.column_vector[1]);
	Y.vector = _mm_mul_ps(Y.vector, _mm_set1_ps(1 / b));
	vec4 Z = linalgNormalizeVec4(matrix.column_vector[2]);
	Z.vector = _mm_mul_ps(Z.vector, _mm_set1_ps(1 / c));
	vec4 T = linalgNormalizeVec4(matrix.column_vector[3]);
	T.vector = _mm_mul_ps(T.vector, _mm_set1_ps(-1));

	mat4 inverse;

	//Column adjustments
	inverse.column[0] = _mm_setr_ps(X.data[0], Y.data[0], Z.data[0], 0);
	inverse.column[1] = _mm_setr_ps(X.data[1], Y.data[1], Z.data[1], 0);
	inverse.column[2] = _mm_setr_ps(X.data[2], Y.data[2], Z.data[2], 0);
	inverse.column[3] = _mm_setr_ps(linalgDotVec4(X, T), linalgDotVec4(Y, T), linalgDotVec4(Z, T), 1);

	return inverse;
}

/*-------- Quaternion Operations ----------*/

quat linalgMakeQuaternionFromComponents(float x, float y, float z, float w) {
	quat q;

	q.vector = _mm_setr_ps(x, y, z, w);

	return q;
}

quat linalgMakeQuaternionFromRotation(float angle, vec3 axis) {
	quat q;

	axis = linalgNormalizeVec3(axis);
	float s = sinf(linalgDeg2Rad(angle / 2));
	float c = cosf(linalgDeg2Rad(angle / 2));

	q.vector = _mm_mul_ps(axis.vector, _mm_set1_ps(s));
	q.data[3] = c;

	return q;
}

quat linalgMakeRotationFromVec2Vec(vec3 a, vec3 b) {

	a = linalgNormalizeVec3(a);
	b = linalgNormalizeVec3(b);

	// a and b might be antiparallel
	if (linalgClose(a, linalgMulVec3(b, -1.0f))) {

		//we want to rotate around a to get to b,
		//pick the least dominant component as the rotation direction
		vec3 ortho = linalgMakeVec3(1, 0, 0);
		if (fabsf(a.data[1]) < fabs(a.data[0])) {
			ortho = linalgMakeVec3(0, 1, 0);
		}
		if (fabsf(a.data[2]) < min(fabs(a.data[0]), fabs(a.data[1]))) {
			ortho = linalgMakeVec3(0, 0, 1);
		}
		vec3 axis = linalgNormalizeVec3(linalgCross(a, ortho));
		return linalgMakeQuaternionFromComponents(
			axis.data[0], 
			axis.data[1], 
			axis.data[2], 
			0.0f
		);
	}

	//Construct the regular quaternion
	vec3 halfVec = linalgNormalizeVec3(linalgAddVec3(a, b));
	vec3 axis = linalgCross(a, halfVec);
	return linalgMakeQuaternionFromComponents(
		axis.data[0], 
		axis.data[1], 
		axis.data[2], 
		linalgDotVec3(a, halfVec)
	);
}

vec3 linalgGetAxisFromQuaternion(quat q) {
	return linalgNormalizeVec3(linalgMakeVec3(q.data[0], q.data[1], q.data[2]));
}

float linalgGetAngleFromQuaternion(quat q) {
	return linalgRad2Deg(2.0f * acosf(q.data[3]));
}

quat linalgAddQuat(quat q1, quat q2) {
	quat result;

	result.vector = _mm_add_ps(q1.vector, q2.vector);

	return result;
}

quat linalgSubQuat(quat q1, quat q2) {
	quat result;

	result.vector = _mm_sub_ps(q1.vector, q2.vector);

	return result;
}

quat linalgMulQuat(quat q, float scalar) {
	quat result;

	result.vector = _mm_mul_ps(q.vector, _mm_set1_ps(scalar));

	return result;
}

float linalgDotQuat(quat q1, quat q2) {
	return q1.data[0] * q2.data[0]
		+ q1.data[1] * q2.data[1]
		+ q1.data[2] * q2.data[2]
		+ q1.data[3] * q2.data[3];
}

bool linalgCloseQuat(quat q1, quat q2) {

	quat displacement = linalgSubQuat(q2, q1);

	return linalgDotQuat(displacement, displacement) < eps;
}

bool linalgQuatSameOrientation(quat q1, quat q2) {

	quat displacement = linalgSubQuat(q1, q2);

	if (linalgDotQuat(displacement, displacement) < eps) {
		return true;
	}

	displacement = linalgAddQuat(q1, q2);

	if (linalgDotQuat(displacement, displacement) < eps) {
		return true;
	}

	return false;
}

quat linalgNormalizeQuat(quat q) {

	float scalar = 1 / sqrtf(linalgDotQuat(q, q));

	return linalgMulQuat(q, scalar);
}

quat linalgGetConjQuat(quat q) {

	return linalgMakeQuaternionFromComponents(
		-q.data[0],
		-q.data[1],
		-q.data[2],
		q.data[3]
	);
}

quat linalgInvQuat(quat q) {

	return linalgMulQuat(linalgGetConjQuat(q), 1 / linalgDotQuat(q, q));
}