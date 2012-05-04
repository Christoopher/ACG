#ifndef COLLISION_DETECTION_H
#define COLLISION_DETECTION_H

#include <vector>
#include <cmath>


#include "OpenGLViewer.h"
#include "RigidBody.h"
#include "Contact.h"
#include "OpenMPDefines.h"


//3rd party
#ifdef _WIN32
#include "armadillo.h"
#else
#include <armadillo>
#endif

#include "omp.h"

typedef std::pair<int, int> RigidBodyIndexPair;

class MinkowskiSet
{
public:
	/*	
	Stores the MK point and the indicies of the the two body vrts 
	MPoint = bodyA.vrts[i] - bodyA.vrts[j] 
	*/
	struct MPoint
	{
		MPoint() : Ai(-1), Bj(-1) //-1 of not set
		{
			p.zeros(3,1);
		}
		arma::vec p; //The Minkowski point
		int Ai,Bj; //Vertex indices of body A and B

		bool operator==(const MPoint & rhs)
		{
			arma::vec d = p - rhs.p;
			if(arma::dot(d,d) < 10e-6)
				return true;
			return false;
		}
	};


	MinkowskiSet(RigidBody &a, RigidBody &b) : _a(a), _b(b)
	{
		//Number of unique pairs
		_points.resize(a.getNumVrts()*b.getNumVrts());
		_stride = b.getNumVrts();
	}

	MPoint *
	getStartingPoint()
	{
		//Do some smart picking of the first point =)

		//e.g closest point to origin
		int minIndex = -1;
		float minDist = 100000;
		for(int i = 0; i < _a.getNumVrts(); ++i)
		{
			for(int j = 0; j < _b.getNumVrts(); ++j)
			{
				int offset = i+_stride*j;
				_points[offset].Ai = i;
				_points[offset].Bj = j;
				_points[offset].p = _a.getWorldVerts()[i] - _b.getWorldVerts()[j];
				float d = arma::norm(_points[offset].p,2);
				if(d < minDist)
				{
					minDist = d;
					minIndex = offset;
				}
			}
		}
		return &_points[minIndex];		
	}

	/*  Returns the farthest away point in the direction d  */
	MPoint *
	support(arma::vec & d)
	{
		int i = _supportMax(d,_a);
		int j = _supportMin(d,_b);
		int offset = i+_stride*j;
		//If not set, calc MPoint
		if(_points[offset].Ai < 0)
		{
			_points[offset].Ai = i;
			_points[offset].Bj = j;
			_points[offset].p = _a.getWorldVerts()[i] - _b.getWorldVerts()[j];
		}

		return &_points[offset];
	}


private:		
	//MPoints: 
	std::vector<MPoint> _points;		
	RigidBody & _a,_b;
	int _stride;
	/* Find the furthest away support vertex in a body in direction d. Returns the index of the vert */
	int
	_supportMax(arma::vec & d, RigidBody & body)
	{
		float max = -100000.0;
		float dist = 0.0;
		int index;
		for(int i = 0; i < body.getNumVrts(); ++i)
		{
			dist = arma::dot(d,body.getWorldVerts()[i]);
			if(dist > max)
			{
				max = dist;	
				index = i;
			}
		}

		return index;
	}

	/* Find the closest support vertex in a body in direction d. Returns the index of the vert */
	int
	_supportMin(arma::vec & d, RigidBody & body)
	{
		float min = 100000.0;
		float dist = 0.0;
		int index;
		for(int i = 0; i < body.getNumVrts(); ++i)
		{
			dist = arma::dot(d,body.getWorldVerts()[i]);
			if(dist < min)
			{
				min = dist;	
				index = i;
			}
		}
		return index;
	}

};

/*
	The simplex is cinstructed by the GJK algorithm.
	It keeps points from the minkowski difference set.
	The simplex is later passed to the EPA algorithm.

	@todo:
		-Might hold extra information about closest point on each n-simplex.
		-Search direction for each n-simplex	
*/
class Simplex
{
public:
	Simplex() : _rank(0) { _points.resize(4); }
	/* convenience method to fetch an MPoint */
	MinkowskiSet::MPoint &
	operator[](int i) { return *_points[i]; }

	std::vector<MinkowskiSet::MPoint *> &
	getImpl() { return _points; }

	/* convenience method tp fetch the last inserted point */
	MinkowskiSet::MPoint &
	last() { return *_points[_rank-1]; }

	int
	rank() { return _rank; }

	/*	
	Inserts a new point in the simplex.
	Return false if point already exists
	Assumes rank <= 3 i.e the last point we insert
	creates a tetra-simplex
	*/
	bool
	insert(MinkowskiSet::MPoint * p)
	{
		for(int i = 0; i < _rank; ++i)
		{
			if(*p == *_points[i])
				return false;
		}

		/*	The point was not previouslt inserted
			So we add it to the simplex
		*/		
		_points[_rank] = p;
		++_rank;

		return true;
	}

		
	/* deletion of one point (lazy) */
	void
	remove(int i)
	{
		//Swap point 'i' to 'last'
		for(int j = i; j < _rank-1; ++j)
			std::swap(_points[j],_points[j+1]);
		--_rank;
	}

	/* only keep point 'i' */
	void
	keep(int i)
	{
		//Swap down point 'i' to 0
		for(int j = i; j > 0; --j)
			std::swap(_points[j],_points[j-1]);
		_rank = 1; //One point left
	}

private:
	//MinkowskiSet::MPoint * _points[4]; //ptr's to MPoints in the minkowski diff
	std::vector<MinkowskiSet::MPoint *> _points; //Might be better to use a vector
	int _rank;
};

class GJK
{
public:
	friend class EPA;

	GJK(RigidBody & a, RigidBody & b) 
		: _bodyA(a) , _bodyB(b), _mk(a,b), _intersects(false), THRESHOLD(10e-9)
	{
		_run();
	}

	bool
	intersects() { return _intersects; }

	Simplex * 
	simplex() { return &_s; }
	

private:
	void
	_run()
	{
		MinkowskiSet::MPoint * A = _mk.getStartingPoint();
		_s.insert(A);
		arma::vec d(-A->p);
		while(true)
		{
			A = _mk.support(d);
			if(arma::dot(A->p,d) < THRESHOLD)
				return;

			if(!_s.insert(A))
			{
				//@todo: Deal with trying to insert an already inserted point
			}
			else
			{
				if(_doSimplex(d))
				{
					_intersects = true;
					return;
				}
			}
		}

	}

	enum _simplextype
	{
		LINE = 2, TRIANGLE = 3, TETRAHEDRON = 4 
	};

	/* returns true if done */
	bool
	_doSimplex(arma::vec &d)
	{
		arma::vec Ao = -_s.last().p;	
		switch(_s.rank())
		{
#pragma region LINE
		case LINE: //Having S = [B] + new point A = line AB
			{
				arma::vec AB = _s[0].p - _s.last().p;
				//If origin lies between point A and B
				if(dot(AB, -_s.last().p) > 0)
				{
					d = cross(cross(AB, Ao), AB);
				}
				else //If origin is beyond A in the direction AB
				{
					d = Ao; 
					_s.keep(1); //Only keep A in the simplex
				}	
			}
			break;
#pragma endregion
#pragma region TRIANGLE
		case TRIANGLE: //Having S = [C B] + new point A = triangle ABC
			{
				arma::vec AB = _s[1].p - _s.last().p;
				arma::vec AC = _s[0].p - _s.last().p;
				arma::vec ABC = cross(AB, AC);
				if(dot(cross(ABC, AC), Ao) > 0)
				{
					if(dot(AC, Ao) > 0) //keep [C A] from [C B A]
					{
						d = cross(cross(AC, Ao), AC);
						//
						_s.remove(1);
					}
					else //STAR
					{
						if(dot(AB, Ao) > 0) //keep [B A] from [C B A]
						{
							d = cross(cross(AB, Ao), AB);

							_s.remove(0);
						}
						else //Only keep [ A ] from [ C B A ]
						{
							d = Ao;
							_s.keep(2);
						}
					}
				}
				else
				{
					if(dot(cross(AB, ABC), Ao) > 0) 
					{
						//STAR
						if(dot(AB, Ao) > 0) //keep [B A] from [C B A]
						{
							d = cross(cross(AB, Ao), AB);

							_s.remove(0);
						}
						else //Only keep [ A ] from [ C B A ]
						{
							d = Ao;
							_s.keep(2);
						}
					}
					else //d is either ABC or -ABC. The simplex is [A B C]
					{
						if(dot(ABC, Ao) > 0)
							d = ABC;
						else
							d = -ABC;
					}
				}	
			}
			break;
#pragma endregion
#pragma region TETRAHEDRON
		case TETRAHEDRON: //Having S = [D C B] + new point A = tetrahedron ABCD
			{
				arma::vec ABC = cross(_s[2].p - _s.last().p, _s[1].p - _s.last().p); //cross(AB, AC) -> triangle ABC
				if(dot(ABC, Ao) > 0) //If on the outside of ABC
				{
					//Remove D
					_s.remove(0);
					return _doSimplex(d);
				}
				else //Inside tetra or outside of ADB or ACD
				{
					arma::vec ADB = cross(_s[0].p - _s.last().p, _s[2].p - _s.last().p); //cross(AD, AB) -> triangle ADB
					if(dot(ADB, Ao) > 0) //Outside ADB
					{
						//Remove C
						_s.remove(1);
						return _doSimplex(d);
					}
					else //Outside ACD or inside tetra
					{
						arma::vec ACD = cross(_s[1].p - _s.last().p, _s[0].p - _s.last().p); //cross(AC, AD) -> triangle ACD
						if(dot(ACD, Ao) > 0) //Outside of ACD
						{
							//remove B;
							_s.remove(2);
							return _doSimplex(d);
						}
						else //inside tetrahedron
							return true; //WE ARE DONE :)
					}
				}
			}
			break;
#pragma endregion
		default:
			std::cout << "There was a 'doSimplex'-error \n";
			break;
		}

		return false;
	}

	RigidBody & _bodyA, _bodyB;
	Simplex _s;
	MinkowskiSet _mk;
	bool _intersects;
	const float THRESHOLD;
};

struct Plane
{
	Plane(MinkowskiSet::MPoint * A, MinkowskiSet::MPoint * B, MinkowskiSet::MPoint * C)
	{
		isPoint[0] = A;
		isPoint[1] = B;
		isPoint[2] = C;

		// create the normal to the current simplex face
		arma::vec n = arma::cross(C->p - A->p, B->p - A->p);

		//Check so that it is pointing away from the origin
		if(arma::dot(n, A->p) < 0) //A could also be B or C (just one point on the plane)
			n = -n;

		// calculate the distance from the origin to the this plane
		n = n/arma::norm(n,2);
		p = dot(A->p, n) * n;
		dist = arma::norm(p,2);

	}

	arma::vec p;
	double dist;
	MinkowskiSet::MPoint * isPoint[3];
};

class EPA
{
public:
	EPA(GJK & gjk) : _gjk(gjk)  
	{
		_run();
	}

	Plane &
	plane() { return _planes[_cpi]; }

private:
	void
	_run()
	{
		MinkowskiSet::MPoint * A;

		//create planes for the starting simplex and insert to planes
		Plane abc(&_gjk._s[0], &_gjk._s[1], &_gjk._s[2]);
		Plane abd(&_gjk._s[0], &_gjk._s[1], &_gjk._s[3]);
		Plane acd(&_gjk._s[0], &_gjk._s[2], &_gjk._s[3]);
		Plane bcd(&_gjk._s[1], &_gjk._s[2], &_gjk._s[3]);
		_planes.push_back(abc);
		_planes.push_back(abd);
		_planes.push_back(acd);
		_planes.push_back(bcd);

		while (true) 
		{
			//Find the closest face
			_findClosestPlane();   //Sets _cpi (closestPlaneIndex) and _closestDist


			A = _gjk._mk.support(_planes[_cpi].p);
			double d = arma::dot(A->p, _planes[_cpi].p);
			if (d - _closestDist < 0.00001) 
			{
				//Finished!
				return; 
			} 
			else //We now have the closest plane, which is the plane from where a new tetrahedron will be "extruded" by 					inserting a new point on the non-origin-side of this plane and adding the 3 new planes. A is the new point.
			{ 
				Plane abc(A, _planes[_cpi].isPoint[0], _planes[_cpi].isPoint[1]);
				Plane abd(A, _planes[_cpi].isPoint[0], _planes[_cpi].isPoint[2]);
				Plane acd(A, _planes[_cpi].isPoint[1], _planes[_cpi].isPoint[2]);
				_planes.erase(_planes.begin() + _cpi);
				_planes.push_back(abc);
				_planes.push_back(abd);
				_planes.push_back(acd);
			}
		}
	}

	void
	_findClosestPlane()
	{
		_closestDist = 10000000.0;

		for (int i = 0; i < _planes.size(); ++i) 
		{
			// check the distance against the other distances
			if (_planes[i].dist < _closestDist) 
			{
				// if this face is closer than previous faces we use it
				_closestDist = _planes[i].dist;
				_cpi = i;
				_csNormal = _planes[i].p;
			}
		}
	}

	
	
	double _closestDist;
	int _cpi; //closestPlaneIndex
	GJK & _gjk;
	arma::vec _csNormal;
	std::vector<Plane> _planes;

};

void
makeContact(Plane & plane,RigidBody & a,RigidBody & b, Contact & c)
{
	/* 
	The simplex is ordered  [ A B C P0 ] 
	Where P0 is the vector from the origin in the minkowski difference to the closest
	point on the hull. ABC is the triangle/plane on which the point exists.
	Also ABC and P0 are orthogonal
	dot(ABC,P0) = 0;

	C
	| \  
	|   \ 
	|     \
	|   P   \
	|         \
	A_ _ _ _ _ _B

	*/

	/* Begin with calculating the barycentric coordinates */

	MinkowskiSet::MPoint *A,*B,*C;
	A = plane.isPoint[0];
	B = plane.isPoint[1];
	C = plane.isPoint[2];

	arma::vec n,na,nb,nc;
	n = na = nb = nc = arma::zeros<arma::vec>(3,1);

	/* N = AB x AC = ABC */
	n = arma::cross(B->p - A->p, C->p - A->p);

	/* na = BC x BP = BCP */
	na = arma::cross(C->p - B->p, plane.p - B->p);

	/* nb = CA x CP = CAP */
	nb = arma::cross(A->p - C->p, plane.p - C->p);

	/* nc = AB x AP = ABP */
	nc = arma::cross(B->p - A->p, plane.p - A->p);

	double nnormsq = n(0)*n(0)+n(1)*n(1)+n(2)*n(2);
	double lA = arma::dot(n,na) / nnormsq;
	double lB = arma::dot(n,nb) / nnormsq;
	double lC = arma::dot(n,nc) / nnormsq;

	/* Calculate world coordinates */
	arma::vec aWorld = lA*a.getWorldVerts()[A->Ai] + lB*a.getWorldVerts()[B->Ai] + lC*a.getWorldVerts()[C->Ai];
	arma::vec bWorld = lA*a.getWorldVerts()[A->Bj] + lB*a.getWorldVerts()[B->Bj] + lC*a.getWorldVerts()[C->Bj];
	arma::vec pWorld = aWorld - bWorld; /* The translation vector */
	c.P = bWorld;
	
	/* Determine vertex-face or edge-edge */

	/*	
		VERTEX-FACE CASE

		(1)
		IF A_vertex <-> B_face

		-> A->Ai == B->Ai == C->Ai

		-> lA == lB == lC == 1/3

		(2)
		IF B_vertex <-> A_face

		-> A->Bi == B->Bi == C->Bi

		-> lA == lB == lC == 1/3

	*/

	std::cout << plane.p(0) << "," << plane.p(1) << "," << plane.p(2) << "\n";

	/* A_vertex <-> B_face */
	if(A->Ai == B->Ai && B->Ai == C->Ai) // ===> A->Ai == B->Ai == C->Ai
	{
	
	
	}
	/* B_vertex <-> A_face */
	else if(A->Ai == B->Ai && B->Ai == C->Ai) // ===> A->Bi == B->Bi == C->Bi
	{
	
	
	}
}

void
narrowPhase(RigidBody & bodyA,RigidBody & bodyB, std::vector<Contact> & contacts)
{
	GJK gjk(bodyA,bodyB);

	if(gjk.intersects())
	{
		EPA epa(gjk);				
		Contact c;
		
		makeContact(epa.plane(),bodyA,bodyB, c);

		bodyA.isColliding = true;
		bodyB.isColliding = true;
		play = false;
	}

}




/***********************************************************************************************/
//									  BROAD PHASE
/***********************************************************************************************/
bool broadPhase(RigidBody & a, RigidBody & b)
{
	//Min max lists
	float maxAx,maxAy,maxAz, minAx, minAy, minAz;
	float maxBx,maxBy,maxBz, minBx, minBy, minBz;
	maxAx = maxAy = maxAz = maxBx = maxBy = maxBz = -10e10;
	minAx = minAy = minAz = minBx = minBy = minBz = +10e10;


	//Calculate bounds for body A
	for(int i = 0; i<a.getWorldVerts().size(); ++i)
	{
		//X axis
	    if(maxAx < a.getWorldVerts()[i](0))
			maxAx = a.getWorldVerts()[i](0);
		if(minAx > a.getWorldVerts()[i](0))
			minAx = a.getWorldVerts()[i](0);
		//Y axis
		if(maxAy < a.getWorldVerts()[i](1))
			maxAy = a.getWorldVerts()[i](1);
		if(minAy > a.getWorldVerts()[i](1))
			minAy = a.getWorldVerts()[i](1);
		//Z axis
		if(maxAz < a.getWorldVerts()[i](2))
			maxAz = a.getWorldVerts()[i](2);
		if(minAz > a.getWorldVerts()[i](2))
			minAz = a.getWorldVerts()[i](2);
	}
	//Calculate bounds for body A
	for(int i = 0; i<b.getWorldVerts().size(); ++i)
	{
		//X axis
		if(maxBx < b.getWorldVerts()[i](0) )
			maxBx = b.getWorldVerts()[i](0);
		if(minBx > b.getWorldVerts()[i](0))
			minBx = b.getWorldVerts()[i](0);
		//Y axis
		if(maxBy < b.getWorldVerts()[i](1))
			maxBy = b.getWorldVerts()[i](1);
		if(minBy > b.getWorldVerts()[i](1))
			minBy = b.getWorldVerts()[i](1);
		//Z axis
		if(maxBz < b.getWorldVerts()[i](2))
			maxBz = b.getWorldVerts()[i](2);
		if(minBz > b.getWorldVerts()[i](2))
			minBz = b.getWorldVerts()[i](2);
	}
	
	//Check for x-axis
	if(maxAx < minBx || maxBx < minAx)
		return false;
	//Check for y-axis
	if(maxAy < minBy || maxBy < minAy)
		return false;
	//Check for z-axis
	if(maxAz < minBz || maxBz < minAz)
		return false;

	//If we come here there is a collision
	return true;

}

/***********************************************************************************************/
//================================= COLLISION DETECTION ========================================
/***********************************************************************************************/
void 
collision_detection(std::vector<RigidBody> & bodies, float t, float dt, std::vector<Contact> & contacts)
{
	//std::vector<RigidBody> backupBodies(bodies.begin(), bodies.end());
	float t0 = t;
	float t1 = t + dt;
	float dt_threshold = 0.5;

	std::vector<RigidBodyIndexPair> colPairs;


#if _DEBUG
	for(int i = 0; i < bodies.size(); ++i)
	{
		bodies[i].isColliding = false;
	}
#endif

	//Broad Phase collision (AABB)
	//Possible collisions are stored in possibleContacts vector
	for (int i = 0; i < bodies.size() - 1; ++i)
	{
		for (int j = i + 1; j < bodies.size(); ++j)
		{
			if(broadPhase(bodies[i],bodies[j]))
			{
				colPairs.push_back(RigidBodyIndexPair(i,j));				
			}
		}
	}

#if NDEBUG
	//For all possible contacts check if we have contact
	//if so add the contact thethe contacts vector
	std::vector<std::vector<Contact> > contactsVector;
	contactsVector.resize(MAX_THREADS);
	int chunksize = colPairs.size() - (colPairs.size() % MAX_THREADS) / MAX_THREADS;
#pragma omp parallel for shedule(dynamic, chunksize)
	for (int i = 0; i < colPairs.size(); ++i)
	{
		narrowPhase(bodies[ colPairs[i].first ], bodies[ colPairs[i].second ], contactsVector[omp_get_thread_num()]);
	}

	for(int i = 0; i < contactsVector.size(); ++i)
	{
		contacts.insert(contacts.end(), contactsVector[i].begin(),contactsVector[i].end());
	}
#else
	for (int i = 0; i < colPairs.size(); ++i)
	{
		narrowPhase(bodies[ colPairs[i].first ], bodies[ colPairs[i].second ], contacts);
	}
#endif

}


#endif