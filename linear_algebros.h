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
#define eps 1e-30f

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

	\param a the first vector
	\param b the second vector
	\returns the angle between the vectors a & b, in degrees
*/
float angleBetweenVectors3(vec3 a, vec3 b);

/**
	Get the projection of one vector onto another.
	Any vector v can be decomposed with regard to another vector u:

		v	= v(parallel with u) + v(perpendicular with u)
			= projection(v onto u) + rejection(v onto u)

	\param incoming the vector to be projected
	\param basis the vector onto which to be projected
	\returns a new vector, parallel with basis, storing the vector projection of incoming onto basis
*/
vec3 project(vec3 incoming, vec3 basis);

/**
	Get the rejection of one vector onto another.
	Any vector v can be decomposed with regard to another vector u:

		v	= v(parallel with u) + v(perpendicular with u)
			= projection(v onto u) + rejection(v onto u)

	\param incoming the vector to be rejected
	\param basis the vector to do the rejecting
	\returns a new vector, orthogonal to basis, storing the vector rejection of incoming from basis
*/
vec3 reject(vec3 incoming, vec3 basis);

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

#endif // !__linear_algebros_h__
