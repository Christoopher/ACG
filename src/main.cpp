//Standard headers
#include <limits>
#include <stdio.h>
#include <iostream>

//3rd party headers
#ifdef _WIN32
	#include "armadillo.h"
#else
	#include <armadillo>
#endif

#include "OpenGLViewer.h"
#include "Physics.h"
#include "RigidBody.h"
#include "Quaternion.h"

int main(void)
{
	physicsInit();

	OpenGl_initViewer(600, 600);
	double t = 0, dt = 1.0/1200.0;
	while(running) {
		if(reset)
			physicsReset();
		reset = false;
		

		if(step || play)
		{
			physicsTick(t,dt);	
		}
			
		//Draw
		OpenGl_drawAndUpdate(running, rigidBodyArray);
		t += dt;
	}
	

	TerminateViewer();
	physicsTerminate();

	return 0;
}



