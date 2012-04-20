//Standard headers
#include <limits>
#include <stdio.h>
#include <iostream>

//3rd party headers
#ifdef _WIN32
	#include "armadillo.h"
#elif
	#include <armadillo>
#endif


#include "OpenGLViewer.h"
#include "Physics.h"
#include "RigidBody.h"
#include "Quaternion.h"

int main(void)
{
	physics_init();

	OpenGl_initViewer(600, 600);
	double t = 0, dt = 1.0/60.0;
	while(running) {
		if(reset)
			physics_reset();
		reset = false;
		

		if(step || play)
		{
			physics_tick(t,dt);			
		}
			
		//Draw
		OpenGl_drawAndUpdate(running, rigidBodyArray, nrOfRigidBodies);
		t += dt;
	}
	

	TerminateViewer();
	physics_terminate();

	return 0;
}



