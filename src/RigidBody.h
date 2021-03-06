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
		Force(double dt, arma::vec X,Quaterniond Q,arma::vec P,arma::vec L, arma::mat R,arma::vec V,arma::vec W,double mass)
{
		arma::vec F = arma::zeros<arma::vec>(3,1);
		F(1) = -mass*9.82;

		return F;
}

arma::vec
Torque(double dt, arma::vec X,Quaterniond Q,arma::vec P,arma::vec L, arma::mat R,arma::vec V,arma::vec W, double mass)
{
		arma::vec F = arma::zeros<arma::vec>(3,1);
		return F;
}


class RigidBody
{
public:
		RigidBody() : mass(916.7), inv_mass(1.0/916.7), isColliding(false)
		{
				inertia.eye(3,3); inertia *= mass/6.0; //r=0.5 => s = 1 => s^2 = 1 (inertia = eye()*s^2*mass / 6.0)
				inv_inertia = arma::inv(inertia);
				force_fun = &Force;
				torque_fun = &Torque;
				X.zeros(3,1);
				V.zeros(3,1);
				W.zeros(3,1);
				R.eye(4,4);
				L.zeros(3,1);
				P.zeros(3,1);
				m_internalForce.zeros(3,1);
				m_internalTorque.zeros(3,1);

				isMovable = true;
				m_externalForce = force_fun(0.0,X,Q,P,L,R,V,W,mass);
				m_externalTorque = torque_fun(0.0,X,Q,P,L,R,V,W,mass);
				float r = 0.5;
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
		init(double x_ = 1.0, double y_ = 1.0, double z_ = 1.0,double density =  916.7)
		{
				sx = x_; sy = y_; sz = z_;
				arma::mat mScale; mScale.zeros(3,3);
				mScale(0,0) = sx;
				mScale(1,1) = sy;
				mScale(2,2) = sz;
				if(isMovable)
				{
						mass = density*sx*sy*sz; //Mass can be considered as density here l00ool0.
						inv_mass = 1.0/mass;
						inertia.eye(3,3);
						//Inertia tensor see:
						//http://en.wikipedia.org/wiki/List_of_moment_of_inertia_tensors
						double height2 = sy*sy;
						double width2  = sx*sx;
						double depth2  = sz*sz;
						double massConstant = mass/12.0;
						inertia(0,0) = massConstant * (height2 + depth2);
						inertia(1,1) = massConstant * (width2  + depth2);
						inertia(2,2) = massConstant * (width2  + height2);
						inv_inertia = arma::inv(inertia);
						m_externalForce = force_fun(0.0,X,Q,P,L,R,V,W,mass);
						//The initial torque should still be
				 }
				else
				{
				  inv_mass = 0;
				  //inv_inertia already 0 in constructor.
				}



				for (unsigned int i = 0; i < _vertices.size(); ++i)
				{
						arma::vec p = mScale * _vertices[i];
						//arma::vec p = R.submat(0,0,2,2)*_vertices[i];
						p = R.submat(0,0,2,2)*p;
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
		arma::vec m_externalForce, m_externalTorque;
		// motion of bodies, then reset to zero for next pass.
		arma::vec m_internalForce, m_internalTorque;
		float sx,sy,sz;
		bool isMovable;
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
				arma::vec,
				double
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
				if(arma::norm(W,2) > 0)
					int halleberry = 0;

				//arma::vec tmp = (0.5 * W);
				// OLD STUFF THAT DOES NOT WORK (maybe...)
				/*
				arma::vec A1DXDT(V);
				Quaterniond A1DQDT = Q*tmp;//tmp;
				*/
				
				arma::vec A1DXDT(V);
				Quaterniond tmpQ = Quaterniond(0.0,W(0),W(1),W(2));
				Quaterniond A1DQDT = (tmpQ*Q)*0.5;
				
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
				tmpQ = Quaterniond(0.0,WN(0),WN(1),WN(2));
				Quaterniond A2DQDT = (tmpQ*QN)*0.5;
				arma::vec A2DPDT = force_fun(tphalfdt,XN,QN,PN,LN,RN,VN,WN,mass);
				arma::vec A2DLDT = torque_fun(tphalfdt,XN,QN,PN,LN,RN,VN,WN,mass);
				XN = X + A2DXDT * halfdt;
				QN = Q + A2DQDT * halfdt;
				PN = P + A2DPDT * halfdt;
				LN = L + A2DLDT * halfdt;
				Convert(QN,PN,LN,RN,VN,WN);

				//Step 3
				arma::vec A3DXDT(VN);
				tmpQ = Quaterniond(0.0,WN(0),WN(1),WN(2));
				Quaterniond A3DQDT = (tmpQ*QN)*0.5;
				arma::vec A3DPDT = force_fun(tphalfdt,XN,QN,PN,LN,RN,VN,WN,mass);
				arma::vec A3DLDT = torque_fun(tphalfdt,XN,QN,PN,LN,RN,VN,WN,mass);
				XN = X + A3DXDT * dt;
				QN = Q + A3DQDT * dt;
				PN = P + A3DPDT * dt;
				LN = L + A3DLDT * dt;
				Convert(QN,PN,LN,RN,VN,WN);

				//Step 4
				arma::vec A4DXDT(VN);
				tmpQ = Quaterniond(0.0,WN(0),WN(1),WN(2));
				Quaterniond A4DQDT = (tmpQ*QN)*0.5;
				arma::vec A4DPDT = force_fun(tpdt,XN,QN,PN,LN,RN,VN,WN,mass);
				arma::vec A4DLDT = torque_fun(tpdt,XN,QN,PN,LN,RN,VN,WN,mass);
				X = X + (A1DXDT + (A2DXDT + A3DXDT)*2.0 + A4DXDT) * sixthdt;
				Q = Q + (A1DQDT + (A2DQDT + A3DQDT)*2.0 + A4DQDT) * sixthdt;
				P = P + (A1DPDT + (A2DPDT + A3DPDT)*2.0 + A4DPDT) * sixthdt;
				L = L + (A1DLDT + (A2DLDT + A3DLDT)*2.0 + A4DLDT) * sixthdt;
				Convert(Q,P,L,R,V,W);

				float detR = 0.0;
				if(abs( (detR = arma::det(R)) - 1.0) > 0.1)
				{
					LOG("Rotation matrix error, abs(arma::det(R) - 1.0) > 0.1");
				}

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
			arma::mat mScale; mScale.zeros(3,3);
			mScale(0,0) = sx;
			mScale(1,1) = sy;
			mScale(2,2) = sz;
			for (int i = 0; i < _vertices.size(); ++i)
			{
					arma::vec p = mScale * _vertices[i];
					//arma::vec p = R.submat(0,0,2,2)*_vertices[i];
					p = R.submat(0,0,2,2)*p;

					

					p += X;
					_verticesWorldCache[i] = p;
			}

		}

		std::vector<arma::vec> _vertices;
		std::vector<arma::vec> _verticesWorldCache;
};



#endif