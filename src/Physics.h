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

const int nrOfRigidBodies = 400;
std::vector<RigidBody> rigidBodyArray;
std::vector<arma::vec> cube;
std::vector<Contact> contacts;

void 
rot(arma::mat & R, int axel, float angle)
{
	if(axel == 1) //x
	{
		R(0,0) = 1.0; 		R(0,1) = 0.0; 			R(0,2) = 0.0;
		R(1,0) = 0.0; 		R(1,1) = cos(angle); 	R(1,2) = -sin(angle);
		R(2,0) = 0.0; 		R(2,1) = sin(angle); 	R(2,2) = cos(angle);
	}
	else if(axel == 2) //y
	{
		R(0,0) = cos(angle); 	R(0,1) = 0.0; 		R(0,2) = sin(angle);
		R(1,0) = 0.0; 		R(1,1) = 1.0; 		R(1,2) = 0.0;
		R(2,0) = -sin(angle); R(2,1) = 0.0; 		R(2,2) = cos(angle);
	}
	else if(axel == 3) //z
	{
		R(0,0) = cos(angle); 	R(0,1) = -sin(angle);	R(0,2) = 0.0;
		R(1,0) = sin(angle); 	R(1,1) = cos(angle);	R(1,2) = 0.0;
		R(2,0) = 0.0; 			R(2,1) = 0.0; 			R(2,2) = 1.0;
	}
	else
	{
		std::cout << "That's not possible. Choose x, y or z as rotation axix!" << endl;
	}
}

void
physicsSetState()
{

//	rigidBodyArray[0].X(0) = rigidBodyArray[0].X(2) = 0;
//	rigidBodyArray[0].X(1) = 5.0;	
//	rigidBodyArray[0].mass = 120;
//	
//
////	double val = 1.0/sqrt(2.0);
//	rigidBodyArray[0].R(0,0) = 0.891006524188368; rigidBodyArray[0].R(0,1) = 0.0; rigidBodyArray[0].R(0,2) = -0.453990499739547;
//	rigidBodyArray[0].R(1,0) = 0.226995249869773; rigidBodyArray[0].R(1,1) = 0.866025403784439; rigidBodyArray[0].R(1,2) = 0.445503262094184;
//	rigidBodyArray[0].R(2,0) =0.39316730585124; rigidBodyArray[0].R(2,1) = -0.5; rigidBodyArray[0].R(2,2) = 0.771634284884801;
//
//	rigidBodyArray[0].Q.FromRotationMatrix(rigidBodyArray[0].R);
//
//	rigidBodyArray[0].force_fun = &Force;
//	rigidBodyArray[0].torque_fun = &Torque;
//	rigidBodyArray[0].init(1,1,1);
	srand ( time(NULL) );
	
	float m = 0.0;
	for(int i = 0; i < nrOfRigidBodies-1; ++i)
	{
		rigidBodyArray[i].X(0) = 2*(rand()/(float)RAND_MAX - 0.5)*4; 
		rigidBodyArray[i].X(2) = 2*(rand()/(float)RAND_MAX - 0.5)*4;
		rigidBodyArray[i].X(1) = 3 + m;
		m += 1.6;

		double val = 1.0/sqrt(2.0);
		rigidBodyArray[i].R(0,0) = 0.891006524188368; rigidBodyArray[i].R(0,1) = 0.0; rigidBodyArray[i].R(0,2) = -0.453990499739547;
		rigidBodyArray[i].R(1,0) = 0.226995249869773; rigidBodyArray[i].R(1,1) = 0.866025403784439; rigidBodyArray[i].R(1,2) = 0.445503262094184;
		rigidBodyArray[i].R(2,0) =0.39316730585124; rigidBodyArray[i].R(2,1) = -0.5; rigidBodyArray[i].R(2,2) = 0.771634284884801;

		
		/*rigidBodyArray[i].R(0,0) = 1; rigidBodyArray[i].R(0,1) = 0.0; rigidBodyArray[i].R(0,2) = 0;
		rigidBodyArray[i].R(1,0) = 0; rigidBodyArray[i].R(1,1) = val; rigidBodyArray[i].R(1,2) = val;
		rigidBodyArray[i].R(2,0) = 0; rigidBodyArray[i].R(2,1) = -val; rigidBodyArray[i].R(2,2) = val;*/

		//rot(rigidBodyArray[i].R,1,45);


		rigidBodyArray[i].Q.FromRotationMatrix(rigidBodyArray[i].R);

		rigidBodyArray[i].force_fun = &Force;
		rigidBodyArray[i].torque_fun = &Torque;

		rigidBodyArray[i].init(1.0,1.0,1.0);
	}




	
	/*
	val = -1.0/sqrt(2.0);
	rigidBodyArray[1].R(0,0) = val; rigidBodyArray[1].R(0,1) = -val; rigidBodyArray[1].R(0,2) = 0.0;
	rigidBodyArray[1].R(1,0) = val; rigidBodyArray[1].R(1,1) = val; rigidBodyArray[1].R(1,2) = 0.0;
	rigidBodyArray[1].R(2,0) = 0.0; rigidBodyArray[1].R(2,1) = 0.0; rigidBodyArray[1].R(2,2) = 1.0;
	*/

	/*rigidBodyArray[1].R(0,0) = 0.891006524188368; rigidBodyArray[1].R(0,1) = 0.0; rigidBodyArray[1].R(0,2) = -0.453990499739547;
	rigidBodyArray[1].R(1,0) = 0.226995249869773; rigidBodyArray[1].R(1,1) = 0.866025403784439; rigidBodyArray[1].R(1,2) = 0.445503262094184;
	rigidBodyArray[1].R(2,0) =0.39316730585124; rigidBodyArray[1].R(2,1) = -0.5; rigidBodyArray[1].R(2,2) = 0.771634284884801;
	rigidBodyArray[1].Q.FromRotationMatrix(rigidBodyArray[1].R);*/
	int index = nrOfRigidBodies-1;
	rigidBodyArray[index].X(0) = 0; rigidBodyArray[index].X(1) = 0; rigidBodyArray[index].X(2) = 0;
	rigidBodyArray[index].inv_inertia = arma::zeros<arma::mat>(3,3);
	rigidBodyArray[index].force_fun = &Force;
	rigidBodyArray[index].torque_fun = &Torque;
	rigidBodyArray[index].isMovable = false;
	rigidBodyArray[index].init(20,1,20);
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

	if(contacts.size() > 0)
	  collisionResponse(contacts);

#pragma omp parallel for
	for (int i = 0; i < rigidBodyArray.size(); ++i)
	{
		//std::cout << rigidBodyArray[0].X(1) << "\n";
		rigidBodyArray[i].update(t,dt);
	}
	//step = false; play = false;
}
