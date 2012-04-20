#ifndef QUATERNION_H
#define QUATERNION_H

//3rd party
#ifdef _WIN32
#include "armadillo.h"
#else
#include <armadillo>
#endif

#include "math.h"

template <typename Real>
class Quaternion
{
public:
	// A quaternion is q = w + x*i + y*j + z*k where (w,x,y,z) is not
	// necessarily a unit-length vector in 4D.

	// Construction.
	Quaternion ();  // uninitialized
	Quaternion (Real w, Real x, Real y, Real z);
	Quaternion (const Quaternion& q);

	// Quaternion for the input rotation matrix.
	//*********Quaternion (const Matrix3<Real>& rot);

	// Quaternion for the rotation of the axis-angle pair.
	//*********Quaternion (const Vector3<Real>& axis, Real angle);

	// Quaternion for the rotation matrix with specified columns.
	//*********Quaternion (const Vector3<Real> rotColumn[3]);

	// Coordinate access as an array:  0 = w, 1 = x, 2 = y, 3 = z.
	//*********inline operator const Real* () const;
	//*********inline operator Real* ();
	//*********inline Real operator[] (int i) const;
	//*********inline Real& operator[] (int i);
	//*********inline Real W () const;
	//*********inline Real& W ();
	//*********inline Real X () const;
	//*********inline Real& X ();
	//*********inline Real Y () const;
	//*********inline Real& Y ();
	//*********inline Real Z () const;
	//*********inline Real& Z ();

	// Assignment.
	 inline Quaternion& operator= (const Quaternion& q);
	 inline Quaternion operator+ (const Quaternion& q) const;
	// Comparison (for use by STL containers).
	//*********inline bool operator== (const Quaternion& q) const;
	//*********inline bool operator!= (const Quaternion& q) const;
	//*********inline bool operator<  (const Quaternion& q) const;
	//*********inline bool operator<= (const Quaternion& q) const;
	//*********inline bool operator>  (const Quaternion& q) const;
	//*********inline bool operator>= (const Quaternion& q) const;

	// Arithmetic operations.
	//*********inline Quaternion operator+ (const Quaternion& q) const;
	//*********inline Quaternion operator- (const Quaternion& q) const;
	inline Quaternion operator* (const Quaternion& q) const;
	inline Quaternion operator* (Real scalar) const;
	//*********inline Quaternion operator/ (Real scalar) const;
	//*********inline Quaternion operator- () const;

	//arma::Col<double>
	inline Quaternion operator* (arma::vec &) const;
	inline friend Quaternion<Real> operator* (arma::vec &, Quaternion &);


	//*********friend Quaternion<Real> operator* (Real scalar,	const Quaternion<Real>& q)
	//{
	//	return q*scalar;
	//}

	// Arithmetic updates.
	//*********inline Quaternion& operator+= (const Quaternion& q);
	//*********inline Quaternion& operator-= (const Quaternion& q);
	//*********inline Quaternion& operator*= (Real scalar);
	//*********inline Quaternion& operator/= (Real scalar);

	// Conversion between quaternions, matrices, and axis-angle.
	void FromRotationMatrix (const arma::mat& rot);
	void ToRotationMatrix (arma::mat& rot) const;
	//void FromRotationMatrix (const arma::vec rotColumn[3]);
	//void ToRotationMatrix (Vector3<Real> rotColumn[3]) const;

	//*********void FromAxisAngle (const Vector3<Real>& axis, Real angle);
	//*********void ToAxisAngle (Vector3<Real>& axis, Real& angle) const;

	// Functions of a quaternion.
	//*********inline Real Length () const;  // length of 4-tuple
	//*********inline Real SquaredLength () const;  // squared length of 4-tuple
	//*********inline Real Dot (const Quaternion& q) const;  // dot product of 4-tuples
	//*********inline Real Normalize (Real epsilon = Math<Real>::ZERO_TOLERANCE);
	//*********Quaternion Inverse () const;  // apply to non-zero quaternion
	//*********Quaternion Conjugate () const;  // negate x, y, and z terms
	//*********Quaternion Exp () const;  // apply to quaternion with w = 0
	//*********Quaternion Log () const;  // apply to unit-length quaternion

	// Rotation of a vector by a quaternion.
	//*********Vector3<Real> Rotate (const Vector3<Real>& vec) const;

	// Spherical linear interpolation.
	//*********Quaternion& Slerp (Real t, const Quaternion& p, const Quaternion& q);
	//*********Quaternion& SlerpExtraSpins (Real t, const Quaternion& p,const Quaternion& q, int extraSpins);

	// Intermediate terms for spherical quadratic interpolation.
	Quaternion& Intermediate (const Quaternion& q0, const Quaternion& q1, const Quaternion& q2);

	// Spherical quadratic interpolation.
	Quaternion& Squad (Real t, const Quaternion& q0, const Quaternion& a0,
		const Quaternion& a1, const Quaternion& q1);

	// Compute a quaternion that rotates unit-length vector V1 to unit-length
	// vector V2.  The rotation is about the axis perpendicular to both V1 and
	// V2, with angle of that between V1 and V2.  If V1 and V2 are parallel,
	// any axis of rotation will do, such as the permutation (z2,x2,y2), where
	// V2 = (x2,y2,z2).
	//*********Quaternion& Align (const Vector3<Real>& v1, const Vector3<Real>& v2);

	// Decompose a quaternion into q = q_twist * q_swing, where q is 'this'
	// quaternion.  If V1 is the input axis and V2 is the rotation of V1 by
	// q, q_swing represents the rotation about the axis perpendicular to
	// V1 and V2 (see Quaternion::Align), and q_twist is a rotation about V1.
	//*********void DecomposeTwistTimesSwing (const Vector3<Real>& v1,Quaternion& twist, Quaternion& swing);

	// Decompose a quaternion into q = q_swing * q_twist, where q is 'this'
	// quaternion.  If V1 is the input axis and V2 is the rotation of V1 by
	// q, q_swing represents the rotation about the axis perpendicular to
	// V1 and V2 (see Quaternion::Align), and q_twist is a rotation about V1.
	//*********void DecomposeSwingTimesTwist (const Vector3<Real>& v1,Quaternion& swing, Quaternion& twist);

	// *** Find closest quaternions with unconstrained angles.

	// Closest quaternion of the form (cx + sx*i).
	Quaternion GetClosestX () const;

	// Closest quaternion of the form (cy + sy*j).
	Quaternion GetClosestY () const;

	// Closest quaternion of the form (cz + sz*k).
	Quaternion GetClosestZ () const;

	// Closest quaternion of the form (cx + sx*i)*(cy + sy*j).
	Quaternion GetClosestXY () const;

	// Closest quaternion of the form (cy + sy*j)*(cx + sx*i).
	Quaternion GetClosestYX () const;

	// Closest quaternion of the form (cz + sz*k)*(cx + sx*i).
	Quaternion GetClosestZX () const;

	// Closest quaternion of the form (cx + sx*i)*(cz + sz*k).
	Quaternion GetClosestXZ () const;

	// Closest quaternion of the form (cy + sy*j)*(cz + sz*k).
	Quaternion GetClosestYZ () const;

	// Closest quaternion of the form (cz + sz*k)*(cy + sy*j).
	Quaternion GetClosestZY () const;

	// Factor to (cx + sx*i)*(cy + sy*j)*(cz + sz*k).
	void FactorXYZ (Real& cx, Real& sx, Real& cy, Real& sy, Real& cz,
		Real& sz);

	// Factor to (cx + sx*i)*(cz + sz*k)*(cy + sy*j).
	void FactorXZY (Real& cx, Real& sx, Real& cz, Real& sz, Real& cy,
		Real& sy);

	// Factor to (cy + sy*j)*(cz + sz*k)*(cx + sx*i).
	void FactorYZX (Real& cy, Real& sy, Real& cz, Real& sz, Real& cx,
		Real& sx);

	// Factor to (cy + sy*j)*(cx + sx*i)*(cz + sz*k).
	void FactorYXZ (Real& cy, Real& sy, Real& cx, Real& sx, Real& cz,
		Real& sz);

	// Factor to (cz + sz*k)*(cx + sx*i)*(cy + sy*j).
	void FactorZXY (Real& cz, Real& sz, Real& cx, Real& sx, Real& cy,
		Real& sy);

	// Factor to (cz + sz*k)*(cy + sy*j)*(cx + sx*i).
	void FactorZYX (Real& cz, Real& sz, Real& cy, Real& sy, Real& cx,
		Real& sx);

	// *** Find closest quaternions with constrained angles.
	class Constraints
	{
	public:
		Constraints ();  // uninitialized
		Constraints (Real minAngle, Real maxAngle);
		void SetAngles (Real minAngle, Real maxAngle);
		bool IsValid (Real x, Real y) const;

		Real MinAngle;       // in [-PI/2,PI/2]
		Real MaxAngle;       // in [m_fMinAngle/2,PI/2]
		Real CosMinAngle;    // = cos(m_fMinAngle)
		Real SinMinAngle;    // = sin(m_fMinAngle)
		Real CosMaxAngle;    // = cos(m_fMaxAngle)
		Real SinMaxAngle;    // = sin(m_fMaxAngle)
		Real DiffCosMaxMin;  // = cos(m_fMaxAngle) - cos(m_fMinAngle)
		Real DiffSinMaxMin;  // = sin(m_fMaxAngle) - sin(m_fMinAngle)
		Real CosAvrAngle;    // = cos((m_fMinAngle + m_fMaxAngle)/2)
		Real SinAvrAngle;    // = sin((m_fMinAngle + mM_faxAngle)/2)
	};

	// Closest constrained quaternion of the form (cx + sx*i).
	Quaternion GetClosestX (const Constraints& xCon) const;

	// Closest constrained quaternion of the form (cy + sy*j).
	Quaternion GetClosestY (const Constraints& yCon) const;

	// Closest constrained quaternion of the form (cz + sz*k).
	Quaternion GetClosestZ (const Constraints& zCon) const;

	// Closest constrained quaternion of the form (cx + sx*i)*(cy + sy*j).
	Quaternion GetClosestXY (const Constraints& xCon,
		const Constraints& yCon) const;

	// Closest constrained quaternion of the form (cy + sy*j)*(cx + sx*i).
	Quaternion GetClosestYX (const Constraints& yCon,
		const Constraints& xCon) const;

	// Closest constrained quaternion of the form (cz + sz*k)*(cx + sx*i).
	Quaternion GetClosestZX (const Constraints& zCon,
		const Constraints& xCon) const;

	// Closest constrained quaternion of the form (cx + sx*i)*(cz + sz*k).
	Quaternion GetClosestXZ (const Constraints& xCon,
		const Constraints& zCon) const;

	// Closest constrained quaternion of the form (cz + sz*k)*(cy + sy*j).
	Quaternion GetClosestZY (const Constraints& zCon,
		const Constraints& yCon) const;

	// Closest constrained quaternion of the form (cy + sy*j)*(cz + sz*k).
	Quaternion GetClosestYZ (const Constraints& yCon,
		const Constraints& zCon) const;

	// Special quaternions.
	static const Quaternion ZERO;
	static const Quaternion IDENTITY;

private:
	// Closest unconstrained quaternion of the form:
	//   (cx + sx*i) when axis = 1,
	//   (cy + sy*j) when axis = 2,
	//   (cz + sz*k) when axis = 3
	Quaternion GetClosest (int axis) const;

	// Closest constrained quaternion of the form:
	//   (cx + sx*i) when axis = 1,
	//   (cy + sy*j) when axis = 2,
	//   (cz + sz*k) when axis = 3
	Quaternion GetClosest (int axis, const Constraints& con) const;

	// Order of storage is (w,x,y,z).
	Real data[4];
};

template<class Real>
Quaternion<Real>::Quaternion()
{
	data[0] = data[1] = data[2] = data[3] = 0.0;
}
template<class Real>
Quaternion<Real>::Quaternion(Real w, Real x, Real y, Real z)
{
	data[0] = w;
	data[1] = x;
	data[2] = y;
	data[3] = z;
}
template<class Real>
Quaternion<Real>::Quaternion(const Quaternion& q)
{
	data[0] = q.data[0];
	data[1] = q.data[1];
	data[2] = q.data[2];
	data[3] = q.data[3];
}

template<class Real>
Quaternion<Real>& Quaternion<Real>::operator= (const Quaternion& q)
{
	if(this==&q)
		return *this;

	data[0] = q.data[0];
	data[1] = q.data[1];
	data[2] = q.data[2];
	data[3] = q.data[3];
	return *this;
}
template <typename Real>
inline Quaternion<Real> Quaternion<Real>::operator* (Real scalar) const
{
	Quaternion result;
	for (int i = 0; i < 4; ++i)
	{
		result.data[i] = scalar*data[i];
	}
	return result;
}

template <class Real>
void Quaternion<Real>::FromRotationMatrix (const arma::mat& rot)
{
	// Algorithm in Ken Shoemake's article in 1987 SIGGRAPH course notes
	// article "Quaternion Calculus and Fast Animation".

	//KOLLA ROW/COLUM-MAJOR LOLEX

	const int next[3] = { 1, 2, 0 };

	Real trace = rot(0,0) + rot(1,1) + rot(2,2);
	Real root;

	if (trace > (Real)0)
	{
		// |w| > 1/2, may as well choose w > 1/2
		root = std::sqrt(trace + (Real)1);  // 2w
		data[0] = ((Real)0.5)*root;
		root = ((Real)0.5)/root;  // 1/(4w)
		data[1] = (rot(2,1) - rot(1,2))*root;
		data[2] = (rot(0,2) - rot(2,0))*root;
		data[3] = (rot(1,0) - rot(0,1))*root;
	}
	else
	{
		// |w| <= 1/2
		int i = 0;
		if (rot(1,1) > rot(0,0))
		{
			i = 1;
		}
		if (rot(2,2) > rot(i,i))
		{
			i = 2;
		}
		int j = next[i];
		int k = next[j];

		root = std::sqrt(rot(i,i) - rot(j,j) - rot(k,k) + (Real)1);
		Real* quat[3] = { &data[1], &data[2], &data[3] };
		*quat[i] = ((Real)0.5)*root;
		root = ((Real)0.5)/root;
		data[0] = (rot(k,j) - rot(j,k))*root;
		*quat[j] = (rot(j,i) + rot(i,j))*root;
		*quat[k] = (rot(k,i) + rot(i,k))*root;
	}
}

template <typename Real>
void Quaternion<Real>::ToRotationMatrix (arma::mat& rot) const
{
	//KOLLA ROW/COLUM-MAJOR LOLEX
	Real twoX  = ((Real)2)*data[1];
	Real twoY  = ((Real)2)*data[2];
	Real twoZ  = ((Real)2)*data[3];
	Real twoWX = twoX*data[0];
	Real twoWY = twoY*data[0];
	Real twoWZ = twoZ*data[0];
	Real twoXX = twoX*data[1];
	Real twoXY = twoY*data[1];
	Real twoXZ = twoZ*data[1];
	Real twoYY = twoY*data[2];
	Real twoYZ = twoZ*data[2];
	Real twoZZ = twoZ*data[3];

	rot(0,0) = (Real)1 - (twoYY + twoZZ);
	rot(0,1) = twoXY - twoWZ;
	rot(0,2) = twoXZ + twoWY;
	rot(1,0) = twoXY + twoWZ;
	rot(1,1) = (Real)1 - (twoXX + twoZZ);
	rot(1,2) = twoYZ - twoWX;
	rot(2,0) = twoXZ - twoWY;
	rot(2,1) = twoYZ + twoWX;
	rot(2,2) = (Real)1 - (twoXX + twoYY);
}
template <class Real>
Quaternion<Real> Quaternion<Real>::operator* (arma::vec & v) const
{
	//Create a q of the vector set w = 0;
	Quaternion<Real> tmp( 0.0,v(0),v(1), v(2));

	return (*this)*tmp;

}

template <typename Real>
inline Quaternion<Real> Quaternion<Real>::operator+ (const Quaternion& q)
	const
{
	Quaternion result;
	for (int i = 0; i < 4; ++i)
	{
		result.data[i] = data[i] + q.data[i];
	}
	return result;
}

template <class Real>
inline Quaternion<Real> operator* (arma::vec & v, Quaternion<Real> & q)
{
	//Create a q of the vector set w = 0;
	Quaternion<Real> tmp( 0.0,v(0),v(1), v(2));
	return tmp*q;
}


template <typename Real>
inline Quaternion<Real> Quaternion<Real>::operator* (const Quaternion& q) const
{
	// NOTE:  Multiplication is not generally commutative, so in most
	// cases p*q != q*p.

	Quaternion result;

	result.data[0] =
		data[0]*q.data[0] -
		data[1]*q.data[1] -
		data[2]*q.data[2] -
		data[3]*q.data[3];

	result.data[1] =
		data[0]*q.data[1] +
		data[1]*q.data[0] +
		data[2]*q.data[3] -
		data[3]*q.data[2];

	result.data[2] =
		data[0]*q.data[2] +
		data[2]*q.data[0] +
		data[3]*q.data[1] -
		data[1]*q.data[3];

	result.data[3] =
		data[0]*q.data[3] +
		data[3]*q.data[0] +
		data[1]*q.data[2] -
		data[2]*q.data[1];

	return result;
}



#endif
