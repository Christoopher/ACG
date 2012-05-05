#include <vector>
#include <cmath>

#include "Contact.h"
#include "CollisionDetection.h"
#include "CollisionResponse.h"

//3rd party
#ifdef _WIN32
#include "armadillo.h"
#else
#include <armadillo>
#endif

const int nrOfRigidBodies = 2;
std::vector<RigidBody> rigidBodyArray;
std::vector<arma::vec> cube;
std::vector<Contact> contacts;

bool first_frame = true;
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

void
physicsSetState()
{
	rigidBodyArray[0].X(0) = rigidBodyArray[0].X(2) = 0;
	rigidBodyArray[0].X(1) = 5.0;
	rigidBodyArray[0].mass = 120;
	

//	double val = 1.0/sqrt(2.0);
	rigidBodyArray[0].R(0,0) = 0.891006524188368; rigidBodyArray[0].R(0,1) = 0.0; rigidBodyArray[0].R(0,2) = -0.453990499739547;
	rigidBodyArray[0].R(1,0) = 0.226995249869773; rigidBodyArray[0].R(1,1) = 0.866025403784439; rigidBodyArray[0].R(1,2) = 0.445503262094184;
	rigidBodyArray[0].R(2,0) =0.39316730585124; rigidBodyArray[0].R(2,1) = -0.5; rigidBodyArray[0].R(2,2) = 0.771634284884801;
	

	rigidBodyArray[0].Q.FromRotationMatrix(rigidBodyArray[0].R);


	rigidBodyArray[0].force_fun = &Force;
	rigidBodyArray[0].torque_fun = &Torque;
	rigidBodyArray[0].init();
	
	/*
	val = -1.0/sqrt(2.0);
	rigidBodyArray[1].R(0,0) = val; rigidBodyArray[1].R(0,1) = -val; rigidBodyArray[1].R(0,2) = 0.0;
	rigidBodyArray[1].R(1,0) = val; rigidBodyArray[1].R(1,1) = val; rigidBodyArray[1].R(1,2) = 0.0;
	rigidBodyArray[1].R(2,0) = 0.0; rigidBodyArray[1].R(2,1) = 0.0; rigidBodyArray[1].R(2,2) = 1.0;
	*/

	rigidBodyArray[1].X(0) = 0; rigidBodyArray[1].X(1) = 0; rigidBodyArray[1].X(2) = 0;
	rigidBodyArray[1].inv_mass = 0;
	rigidBodyArray[1].inv_inertia = arma::zeros<arma::mat>(3,3);
	rigidBodyArray[1].force_fun = &Force;
	rigidBodyArray[1].torque_fun = &Torque;
	rigidBodyArray[1].init();
}

void
physicsInit()
{

	createCubeList(1.0,cube);
	rigidBodyArray.resize(nrOfRigidBodies);

	physicsSetState();
}

void 
physicsReset()
{
	for (int i = 0; i < rigidBodyArray.size(); ++i)
	{
		rigidBodyArray[i].reset();
	}

	physicsSetState();
}

void
physicsTerminate()
{
	//delete [] rigidBodyArray;
}

void 
physicsTick(double t, double dt)
{
	contacts.clear();
	collision_detection(rigidBodyArray, t, dt, contacts);
	setContacts(contacts);
	//collisionResponse(rigidBodyArray,contacts);
	for (int i = 0; i < rigidBodyArray.size(); ++i)
	{
		//std::cout << rigidBodyArray[0].X(1) << "\n";
		rigidBodyArray[i].update(t,dt);
	}
}
