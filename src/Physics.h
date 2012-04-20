#include <vector>
#include <cmath>

//3rd party
#ifdef _WIN32
#include "armadillo.h"
#else
#include <armadillo>
#endif


const int nrOfRigidBodies = 2;
RigidBody * rigidBodyArray;
std::vector<arma::vec> cube;

bool first_frame = true;
arma::vec Force	(double dt, arma::vec X,Quaterniond Q,arma::vec P,arma::vec L, arma::mat R,arma::vec V,arma::vec W)
{
	arma::vec F = arma::zeros<arma::vec>(3,1);
	F(1) = -120*9.82;

	return F;
}

arma::vec Torque(double dt, arma::vec X,Quaterniond Q,arma::vec P,arma::vec L, arma::mat R,arma::vec V,arma::vec W)
{
	arma::vec F = arma::zeros<arma::vec>(3,1);
	return F;
}


void physics_set_state()
{
	rigidBodyArray[0].X(0) = rigidBodyArray[0].X(2) = 0;
	rigidBodyArray[0].X(1) = 10.0;
	rigidBodyArray[0].mass = 120;
	

	double val = 1.0/sqrt(2.0);
	rigidBodyArray[0].R(0,0) = 0.891006524188368; rigidBodyArray[0].R(0,1) = 0.0; rigidBodyArray[0].R(0,2) = -0.453990499739547;
	rigidBodyArray[0].R(1,0) = 0.226995249869773; rigidBodyArray[0].R(1,1) = 0.866025403784439; rigidBodyArray[0].R(1,2) = 0.445503262094184;
	rigidBodyArray[0].R(2,0) =0.39316730585124; rigidBodyArray[0].R(2,1) = -0.5; rigidBodyArray[0].R(2,2) = 0.771634284884801;
	

	rigidBodyArray[0].Q.FromRotationMatrix(rigidBodyArray[0].R);
	

	rigidBodyArray[0].force_fun = &Force;
	rigidBodyArray[0].torque_fun = &Torque;

	
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
}


void physics_init()
{

	create_cube_list(1.0,cube);
	rigidBodyArray = new RigidBody[nrOfRigidBodies];
	physics_set_state();
}

void physics_reset()
{
	for (int i = 0; i < nrOfRigidBodies; ++i)
	{
		rigidBodyArray[i].reset();
	}

	physics_set_state();
}


void physics_terminate()
{
	delete [] rigidBodyArray;
}


void resolve_collision(RigidBody &A, RigidBody &B, arma::vec & P, arma::vec & N)
{
	const double e = 0.6;  //energy loss;


	arma::vec rA,rB,kA,kB,uA,uB,impulse;
	rA = rB = kA= kB = uA = uB= impulse = arma::zeros<arma::vec>(3,1);
	

	//Both the point P and X is in world coordinates
	
	rA = P - A.X;
	rB = P - B.X;
	kA = arma::cross(rA,N);
	kB = arma::cross(rB,N);

	uA = A.inv_inertia * kA;
	uB = B.inv_inertia * kB;

	using arma::dot;
	double numer = -(1.0+e) * (dot(N,A.V - B.V) + dot(A.W,kA) - dot(B.W, kB));
	double denom = A.inv_mass + B.inv_mass + dot(kA,uA) + dot(kB,uB);
	double f = numer / denom;
		
	impulse = f * N;

	//Apply impulse	
	A.P += impulse;
	B.P -= impulse;
	A.L += f*kA;
	B.L -= f*kB;

	//update vel/ang

	A.V = A.P*A.inv_mass;
	//B.V = A.P*B.inv_mass; this is what it says in the book but probably wrong!!
	B.V = B.P*B.inv_mass;

	A.W = A.inv_inertia.submat(0,0,2,2) * A.L;
	B.W = B.inv_inertia.submat(0,0,2,2) * B.L;

	//A.W += f * uA;
	//B.W = f * uB; //could be a minus sign

}


void collision_detection()
{
	arma::vec N = arma::zeros<arma::vec>(3,1);
	arma::vec WP = arma::zeros<arma::vec>(3,1);
	N(1) = 1.0;
	
	for (int i = 0; i < cube.size(); ++i)
	{
		WP = cube[i];

		WP = rigidBodyArray[0].R.submat(0,0,2,2) * WP;
		WP += rigidBodyArray[0].X;

		//Vertex i värld
		if(WP(1) <= 0.0) 
		{
			//nu ligger WP under eller på
			double move = abs(WP(1));
			rigidBodyArray[0].X(1) += move;
			WP(1) += move;
			resolve_collision(rigidBodyArray[0],rigidBodyArray[1], WP ,N);
		}
	}


}


void physics_tick(double t, double dt)
{
	for (int i = 0; i < nrOfRigidBodies; ++i)
	{
		//std::cout << rigidBodyArray[0].X(1) << "\n";
		collision_detection();
		rigidBodyArray[i].update(t,dt);
	}
}
