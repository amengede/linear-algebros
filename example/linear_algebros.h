#include <intrin.h>
#include <stdbool.h>

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
		__m256 chunk[2];
		__m128 column[4];
		vec4 column_vector[4];
		float data[16];
	};
} mat4;

typedef struct quaternion {
	union {
		__m128 vector;
		float data[4];
	};
} quat;

/*-------- Conversions        ----------*/

#define pi 3.14159265359f
#define eps 1e-10f

/**
	Convert from degrees to radians.

	\param angle the angle to convert (in degrees)
	\returns the angle converted to radians
*/
float linalgDeg2Rad(float angle);

/**
	Convert from radians to degrees.

	\param angle the angle to convert (in radians)
	\returns the angle converted to degrees
*/
float linalgRad2Deg(float angle);

/*-------- Vec3 Operations    ----------*/

/**
	Construct a vec3 from individual floating point components.
*/
vec3 linalgMakeVec3(float x, float y, float z);

/**
	Compute the dot product between two vec3s,
	as with any dot product, the order doesn't matter.

	\param a the first vector
	\param b the second vector
	\returns the dot product: a.b
*/
float linalgDotVec3(vec3 a, vec3 b);

/**
	Compute the cross product between two vec3s

	\param a the first vector
	\param b the second vector
	\returns a new vector storing the cross product: c = axb
*/
vec3 linalgCross(vec3 a, vec3 b);

/**
	\param a the vector to normalize, cannot have zero magnitude
	\returns a new vec3, being parallel with the input a, having unit length.
*/
vec3 linalgNormalizeVec3(vec3 a);

/**
	Compute a vector subtraction.

	\param a the first vector
	\param b the second vector
	\returns a new vector storing the difference: c = a - b
*/
vec3 linalgSubVec3(vec3 a, vec3 b);

/**
	Compute a vector addition.

	\param a the first vector
	\param b the second vector
	\returns a new vector storing the sum: c = a + b
*/
vec3 linalgAddVec3(vec3 a, vec3 b);

/**
	Compute a vector scalar multiplication.

	\param a the vector
	\param scalar the scalar
	\returns a new vector storing the scaled vector: c = scalar * a
*/
vec3 linalgMulVec3(vec3 a, float scalar);

/**
	Get the angle between two vectors.

	\param a the first vector, cannot have zero magnitude
	\param b the second vector, cannot have zero magnitude
	\returns the angle between the vectors a & b, in degrees
*/
float linalgAngleBetweenVectors3(vec3 a, vec3 b);

/**
	Get the projection of one vector onto another.
	Any vector v can be decomposed with regard to another vector u:

		v	= v(parallel with u) + v(perpendicular with u)
			= projection(v onto u) + rejection(v onto u)

	\param incoming the vector to be projected
	\param basis the vector onto which to be projected, cannot have zero magnitude
	\returns a new vector, parallel with basis, storing the vector projection of incoming onto basis
*/
vec3 linalgProject(vec3 incoming, vec3 basis);

/**
	Get the rejection of one vector onto another.
	Any vector v can be decomposed with regard to another vector u:

		v	= v(parallel with u) + v(perpendicular with u)
			= projection(v onto u) + rejection(v onto u)

	\param incoming the vector to be rejected
	\param basis the vector to do the rejecting, cannot have zero magnitude
	\returns a new vector, orthogonal to basis, storing the vector rejection of incoming from basis
*/
vec3 linalgReject(vec3 incoming, vec3 basis);

/**
	Compute a vector reflection.

	\param incident a direction vector incident to (pointing towards) the point of impact.
	\param normal the normal vector about which to reflect. Must have unit length.
	\returns a new vector representing the direction after reflecting.
*/
vec3 linalgReflect(vec3 incident, vec3 normal);

/**
	Linearly interpolate between two vectors.

	\param a the first vector
	\param b the second vector
	\param t the interpolation parameter. Typically between 0 and 1, though this isn't enforced
	\returns a new vector, being a linear interpolation between a and b.
*/
vec3 linalgLerpVec3(vec3 a, vec3 b, float t);

/**
	Spherical Linear interpolation between two vectors.
	lerp will take a straight line between vectors, on the other hand,
	slerp interpolates angle-wise, in a rotational sense.

	\param a the first vector, should be normalized.
	\param b the second vector, should be normalized.
	\param t the interpolation parameter. Typically between 0 and 1, though this isn't enforced
	\returns a new vector, being a linear interpolation between a and b.
*/
vec3 linalgSlerpVec3(vec3 a, vec3 b, float t);

/**
	Normalized Linear interpolation between two vectors.
	Normalizing the result of lerp will approximate slerp.

	\param a the first vector, should be normalized.
	\param b the second vector, should be normalized.
	\param t the interpolation parameter. Typically between 0 and 1, though this isn't enforced
	\returns a new vector, being a linear interpolation between a and b.
*/
vec3 linalgNlerpVec3(vec3 a, vec3 b, float t);

/**
	Indicates whether two vectors are within epsilon of one another.
*/
bool linalgCloseVec3(vec3 a, vec3 b);

/*-------- Vector4 Operations ----------*/

/*
	Under the hood, vec3 is really a vec4 with zeroed out w component,
	however as convention has it right now, it's being treated as a seperate type.
*/

/**
	\returns a normalized copy of the given vector.
*/
vec4 linalgNormalizeVec4(vec4 a);

/**
	\returns the dot product result = a.b
*/
float linalgDotVec4(vec4 a, vec4 b);

/*-------- Matrix4 Operations ----------*/

/**
	\returns a new 4x4 matrix storing the identity transform.
*/
mat4 linalgMakeIdentity4();

/**
	Make a perspective projection matrix.

	\param fovy the field of view angle of the frustrum (in degrees)
	\param aspect the aspect ratio width/height
	\param near the near view distance of the frustrum
	\param far the far distance of the frustrum
	\returns a new mat4 representing the perspective projection transform
*/
mat4 linalgMakePerspectiveProjection(float fovy, float aspect, float near, float far);

/**
	Make a view matrix (translates and rotates the world around the 
	given reference point)

	\param eye the position of the viewer
	\param target the position the viewer is looking at
	\param up the up direction from the viewer's point of reference
	\returns a new mat4 representing the view transform
*/
mat4 linalgMakeLookAt(vec3 eye, vec3 target, vec3 up);

/**
	Make a translation transform matrix.

	\param translation the displacement to apply
	\returns a new mat4 representing the transform
*/
mat4 linalgMakeTranslation(vec3 translation);

/**
	Make a rotation around the x-axis.

	\param angle the angle to rotate by (in degrees)
	\returns a new mat4 representing the transform
*/
mat4 linalgMakeXRotation(float angle);

/**
	Make a rotation around the y-axis.

	\param angle the angle to rotate by (in degrees)
	\returns a new mat4 representing the transform
*/
mat4 linalgMakeYRotation(float angle);

/**
	Make a rotation around the z-axis.

	\param angle the angle to rotate by (in degrees)
	\returns a new mat4 representing the transform
*/
mat4 linalgMakeZRotation(float angle);

/**
	Transform a vector by a matrix.

	\param m the matrix to apply
	\param v the vector to transform
	\returns a new vec4 representing the matrix multiplication: result = m*v
*/
vec4 linalgMulMat4Vec4(mat4 m, vec4 v);

/**
	Multiply two matrices

	\param m1 the original matrix
	\param m2 the new matrix to multiply onto m1
	\returns a new mat4 representing the matrix multiplication: m3 = m2*m1
*/
mat4 linalgMulMat4Mat4(mat4 m1, mat4 m2);

/**
	\returns the matrix sum m3 = m1 + m2
*/
mat4 linalgAddMat4(mat4 m1, mat4 m2);

/**
	\returns the scalar multiplication result = scalar * matrix
*/
mat4 linalgMulMat4Scalar(mat4 matrix, float scalar);

/**
	Blend (linearly interpolate) two matrices.

	\param m1 the start matrix (t = 0)
	\param m2 the end matrix (t = 1)
	\param t the interpolation parameter
	\returns the result m3 = m1 + t * (m2 - m1)
*/
mat4 linalgLerpMat4(mat4 m1, mat4 m2, float t);

/**
	\returns a transposed copy of the given matrix
*/
mat4 linalgTranspose(mat4 matrix);

/**
	Compute a transform matrix inverse.

	General matrix inverses are computationally intense, however
	most transform matrices can be expressed in the form:

	M = (aX | bY | cZ | I)

	where:

		a: x axis scaling factor
		X: forwards rotation vector
		b: y axis scaling factor
		Y: right rotation vector
		c: z axis scale factor
		Z: up rotation vector
		I: [0 0 0 1]

	Matrices in this form have a closed form inverse.

	Source: https://lxjk.github.io/2017/09/03/Fast-4x4-Matrix-Inverse-with-SSE-SIMD-Explained.html

	\param matrix the matrix to invert
	\returns the inverse
*/
mat4 linalgTransformInverse(mat4 matrix);

/*-------- Quaternion Operations ----------*/

/**
	\returns a quaternion made from individual components.
*/
quat linalgMakeQuaternionFromComponents(float x, float y, float z, float w);

/**
	Make a quaternion from a rotation operation.

	\param angle the rotation angle (in degrees)
	\param axis the axis of rotation
	\returns the corresponding quaternion
*/
quat linalgMakeQuaternionFromRotation(float angle, vec3 axis);

/**
	Make a quaternion tracking a rotation from vector a to vector b.
*/
quat linalgMakeRotationFromVec2Vec(vec3 a, vec3 b);

/**
	\returns the quaternion's axis of rotation
*/
vec3 linalgGetAxisFromQuaternion(quat q);

/**
	\returns the quaternion's angle, in degrees
*/
float linalgGetAngleFromQuaternion(quat q);

/**
	\returns the sum of two quaternions
*/
quat linalgAddQuat(quat q1, quat q2);

/**
	\returns the difference of two quaternions
*/
quat linalgSubQuat(quat q1, quat q2);

/**
	\returns a scaled copy of the quaternion
*/
quat linalgMulQuat(quat q, float scalar);

/**
	\returns the dot product of two quaternions
*/
float linalgDotQuat(quat q1, quat q2);

/**
	\returns whether two quaternions are sufficiently close.
*/
bool linalgCloseQuat(quat q1, quat q2);

/**
	It's possible for two quaternions to be the same, but mathematically different.
	Ie. a rotation in the opposite angle through the opposite axis is actually still the same.

	\returns whether two quaternions have the same orientation.
*/
bool linalgQuatSameOrientation(quat q1, quat q2);

/**
	\returns a normalized quaternion
*/
quat linalgNormalizeQuat(quat q);

/**
	\returns the conjugate quaternion, a rotation of the same angle, around the opposite axis.
*/
quat linalgGetConjQuat(quat q);

/**
	\returns the inverse of the given quaternion.
*/
quat linalgInvQuat(quat q);

#endif // !__linear_algebros_h__
