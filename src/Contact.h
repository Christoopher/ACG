#ifndef CONTACT_H
#define CONTACT_H

#include <cmath>

//3rd party
#ifdef _WIN32
#include "armadillo.h"
#else
#include <armadillo>
#endif

#include "RigidBody.h"

struct Contact
{
	Contact() : isVFContact(true), A(NULL), B(NULL)
	{
		P.zeros(3,1);
		N.zeros(3,1);
		EA.zeros(3,1);
		EB.zeros(3,1);
		isVFContact = true;

		pOnA.zeros(3,1);

	}

	RigidBody * A; //body containing vertex
	RigidBody * B; //body containing face
	arma::vec P; //contact point, //point on body b from EPA
	arma::vec N; //outward normal of face
	arma::vec EA; //edge from A
	arma::vec EB; //edge from B
	bool isVFContact; //true if vertex-face, false if edge-edge


	arma::vec pOnA;
};



#endif