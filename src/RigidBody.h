#ifndef RIGIDBODY_H
#define RIGIDBODY_H
#include "UTLog.h"
#include <assert.h>
//3rd party
#ifdef _WIN32
#include "armadillo.h"
#else
#include <armadillo>
#endif

#include "Quaternion.h"

typedef Quaternion<double> Quaterniond;

arma::vec
Force(double dt, arma::vec X,Quaterniond Q,arma::vec P,arma::vec L, arma::mat R,arma::vec V,arma::vec W)
{
    arma::vec F = arma::zeros<arma::vec>(3,1);
    F(1) = -120*9.82;

    return F;
}

arma::vec
Torque(double dt, arma::vec X,Quaterniond Q,arma::vec P,arma::vec L, arma::mat R,arma::vec V,arma::vec W)
{
    arma::vec F = arma::zeros<arma::vec>(3,1);
    return F;
}


class RigidBody
{
public:
	RigidBody() : mass(120.0), inv_mass(1.0/120.0), isColliding(false)
	{
		inertia.eye(3,3); inertia *= mass/12.0;
		inv_inertia = arma::inv(inertia);
        force_fun = &Force;
        torque_fun = &Torque;
		X.zeros(3,1);
		V.zeros(3,1);
		W.zeros(3,1);
		R.eye(4,4);
		L.zeros(3,1);	
		P.zeros(3,1);
     /*   if(force_fun==NULL)
        {
            LOG("force_fun is null!");
            assert(force_fun==NULL);
        }
        else if(torque_fun==NULL)
        {
            LOG("torque_fun is null!");
            assert("torque_fun==NULL");
        }*/
        m_externalForce = force_fun(0.0,X,Q,P,L,R,V,W);
        m_externalTorque = torque_fun(0.0,X,Q,P,L,R,V,W);
        float r = 1.0;
		arma::vec v = arma::zeros<arma::vec>(3,1);
		v(0) = -r; v(1) = -r; v(2) = -r; _vertices.push_back(v);
		v(0) = r; v(1) = -r; v(2) = -r; _vertices.push_back(v);
		v(0) = r; v(1) = -r; v(2) = r; _vertices.push_back(v);
		v(0) = -r; v(1) = -r; v(2) = r; _vertices.push_back(v);
		v(0) = -r; v(1) = r; v(2) = -r; _vertices.push_back(v);	
		v(0) = r; v(1) = r; v(2) = -r; _vertices.push_back(v);
		v(0) = r; v(1) = r; v(2) = r; _vertices.push_back(v);
		v(0) = -r; v(1) = r; v(2) = r; _vertices.push_back(v);
	}
	
	void
        AppendInternalForce (arma::vec& intForce)
        {
        m_internalForce += intForce;
        }

        void
        AppendInternalTorque (arma::vec& intTorque)
        {
        m_internalTorque += intTorque;
        }

	void
	init()
	{
		for (int i = 0; i < _vertices.size(); ++i)
		{
			arma::vec p = R.submat(0,0,2,2)*_vertices[i];
			p += X;
			_verticesWorldCache.push_back(p);
		}
	}



	void reset()
	{
		X = V = W = L = P = arma::zeros<arma::vec>(3,1);
		R = arma::eye(4,4);
	}


	//Constants
	double mass, inv_mass;
	arma::mat inertia, inv_inertia;

	//State Variables
	arma::vec X, P, L; //Pos, lin mom., ang mom.
	Quaterniond Q; //Orientation
	//Derived State VAriables
	arma::mat R; // Orientiation
	arma::vec V, W; //lin vel, ang vel
	
	// external force/torque at current time of simulation
    arma::vec3 m_externalForce, m_externalTorque;
    // motion of bodies, then reset to zero for next pass.
    arma::vec3 m_internalForce, m_internalTorque;
 	
	
	bool isColliding;

	typedef arma::vec (*Function)
	(
		double,
		arma::vec,
		Quaterniond,
		arma::vec,
		arma::vec,
		arma::mat,
		arma::vec,
		arma::vec
	);

	Function torque_fun;
	Function force_fun;

	void 
	update(double t, double dt)
	{
		double halfdt = 0.5 * dt, sixthdt = dt / 6.0;
		double tphalfdt = t + halfdt, tpdt = t + dt;

		//Temp States
		arma::vec XN, PN, LN, VN, WN;
		XN = PN = LN = VN = WN = arma::zeros<arma::vec>(3,1);
		Quaterniond QN;
		arma::mat RN = arma::zeros<arma::mat>(4,4);


		//STEp1
		arma::vec A1DXDT(V);
		arma::vec tmp = (0.5 * W);
        Quaterniond A1DQDT = Q*tmp;//tmp;
        //arma::vec A1DPDT = force_fun(t,X,Q,P,L,R,V,W);
        //arma::vec A1DLDT = torque_fun(t,X,Q,P,L,R,V,W);

        arma::vec A1DPDT = m_externalForce  + m_internalForce;
        arma::vec A1DLDT = m_externalTorque + m_internalTorque;
		XN = X + A1DXDT * halfdt;
		QN = Q + A1DQDT * halfdt;
		PN = P + A1DPDT * halfdt;
		LN = L + A1DLDT * halfdt;
		Convert(QN,PN,LN,RN,VN,WN);

		//new stuff------------------------

        m_internalForce  = arma::zeros<arma::vec>(3);
        m_internalTorque = arma::zeros<arma::vec>(3);
                //----------------------------------

		//Step 2
		arma::vec A2DXDT(VN);
		tmp = 0.5*WN;
		Quaterniond A2DQDT = QN*tmp;
		arma::vec A2DPDT = force_fun(tphalfdt,XN,QN,PN,LN,RN,VN,WN);
		arma::vec A2DLDT = torque_fun(tphalfdt,XN,QN,PN,LN,RN,VN,WN);
		XN = X + A2DXDT * halfdt;
		QN = Q + A2DQDT * halfdt;
		PN = P + A2DPDT * halfdt;
		LN = L + A2DLDT * halfdt;
		Convert(QN,PN,LN,RN,VN,WN);

		//Step 3
		arma::vec A3DXDT(VN);
		tmp = 0.5*WN;
		Quaterniond A3DQDT = QN*tmp;
		arma::vec A3DPDT = force_fun(tphalfdt,XN,QN,PN,LN,RN,VN,WN);
		arma::vec A3DLDT = torque_fun(tphalfdt,XN,QN,PN,LN,RN,VN,WN);
		XN = X + A3DXDT * dt;
		QN = Q + A3DQDT * dt;
		PN = P + A3DPDT * dt;
		LN = L + A3DLDT * dt;
		Convert(QN,PN,LN,RN,VN,WN);

		//Step 4
		arma::vec A4DXDT(VN);
		tmp = 0.5*WN;
		Quaterniond A4DQDT = QN*tmp;
		arma::vec A4DPDT = force_fun(tpdt,XN,QN,PN,LN,RN,VN,WN);
		arma::vec A4DLDT = torque_fun(tpdt,XN,QN,PN,LN,RN,VN,WN);
		X = X + (A1DXDT + (A2DXDT + A3DXDT)*2.0 + A4DXDT) * sixthdt;
        Q = Q + (A1DQDT + (A2DQDT + A3DQDT)*2.0 + A4DQDT) * sixthdt;
		P = P + (A1DPDT + (A2DPDT + A3DPDT)*2.0 + A4DPDT) * sixthdt;
		L = L + (A1DLDT + (A2DLDT + A3DLDT)*2.0 + A4DLDT) * sixthdt;
		Convert(Q,P,L,R,V,W);
		

		//Recalc world transform
		_calcWorldTransform();
	}

	void 
	Convert(const Quaterniond & Q,const arma::vec & P,const arma::vec &L, arma::mat &R, arma::vec &V, arma::vec &W) const
	{
		Q.ToRotationMatrix(R);
		V = inv_mass * P;
		W = R.submat(0,0,2,2)*inv_inertia*arma::trans(R.submat(0,0,2,2))*L;
	}

	std::vector<arma::vec> &
	getWorldVerts()	{ return _verticesWorldCache; }

	int
	getNumVrts() { return _vertices.size(); }

protected:

private:

	void
	_calcWorldTransform()
	{
		for (int i = 0; i < _vertices.size(); ++i)
		{
			arma::vec p = R.submat(0,0,2,2)*_vertices[i];
			p += X;
			_verticesWorldCache[i] = p;
		}
	}

	std::vector<arma::vec> _vertices;
	std::vector<arma::vec> _verticesWorldCache;
};



#endif
