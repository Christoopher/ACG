#ifndef _OPENGL_VIEWER_
#define _OPENGL_VIEWER_

//C/C++ headers
#include <stdlib.h>
#include <stdio.h>
#include <vector>

//Opengl and glew headers

#ifdef _WIN32
#include "glew.h"
#include "wglew.h"
#else
#include "GL/glew.h"
#include "GL/gl.h"
#endif


//3rd party
#ifdef _WIN32
#include "armadillo.h"
#include "glfw.h"
#include "glut.h"
#else
#include <armadillo>
#include "GL/glfw.h"
#include "GL/glut.h"
#endif

#include "RigidBody.h"
#include "Contact.h"

//FUL VARRE
std::vector<Contact> * conts;

//----------------------------------------------------------------------------//
// Variables declaration
//----------------------------------------------------------------------------//

//Window properties and GUI-related
int winw, winh;
int lastmousex, lastmousey, lastwheelpos;
float posDx = 0.0f, posDy = 0.0f, zoom = 0.0f, rotDx = 0.0f, rotDy = 0.0f;
double t0 = 0;
int frames = 0;
bool running = true;
GLfloat edge = 3.0f;
bool step = false, reset = false, showgrid = false, play = false;
bool up_is_down;
bool down_is_down;
bool shiftdown = false;
bool gjkdraw = true;
bool wireframe = false;
bool skip = false;

void setContacts(std::vector<Contact> & c)
{
	conts = &c;
}



//----------------------------------------------------------------------------//
// Print opengl error to console
//----------------------------------------------------------------------------//
void getOpenGLError()
{

	GLenum errCode;
	const GLubyte *errString;

	if ((errCode = glGetError()) != GL_NO_ERROR) {
		errString = gluErrorString(errCode);
		fprintf (stderr, "OpenGL Error: %s\n", errString);
	}
}

float getFPS(void)
{
	double t = glfwGetTime();
	return (float)((double)frames / (t - t0));

}
//----------------------------------------------------------------------------//
// Display fps, winh, winw, zoom and title
//----------------------------------------------------------------------------//
void showFPS(int winw, int winh, float zoom) {

	static char titlestr[200];
	double t, fps;
	// Get current time
	t = glfwGetTime();  // Number of seconds since glfwInit()
	// If one second has passed, or if this is the very first frame
	if( (t - t0) > 1.0 || frames == 0 )
	{
		fps = (double)frames / (t - t0);
		sprintf(titlestr, "FLIP3D, %dx%d pixels, %.1fx zoom, %.1f FPS -> %.1f Mpixels/s",
			winw, winh, zoom, fps, winw*winh*fps*1e-6);
		glfwSetWindowTitle(titlestr);
		t0 = t;
		frames = 0;
	}
	frames ++;
}

//----------------------------------------------------------------------------//
// Cleans up the OpenGL viewer: Destroy window, destroys shaders, fbos, vbos, etc...
//----------------------------------------------------------------------------//
void TerminateViewer()
{
	//Terminate the window
	glfwTerminate();
}

//----------------------------------------------------------------------------//
// Responsible for updating OpenGL due to screen resizing
//----------------------------------------------------------------------------//
void Resize()
{
	glfwGetWindowSize(&winw, &winh);
	glViewport(0, 0, winw, winh);

	glMatrixMode(GL_PROJECTION);

	glLoadIdentity();
	gluPerspective(65.0, (float)winw / winh, 0.1, 1000);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

}

void createCubeList(float r, std::vector<arma::vec> &cube) {

	arma::vec v = arma::zeros<arma::vec>(3,1);
	v(0) = -r; v(1) = -r; v(2) = -r; cube.push_back(v);
	v(0) = r; v(1) = -r; v(2) = -r; cube.push_back(v);
	v(0) = r; v(1) = -r; v(2) = r; cube.push_back(v);
	v(0) = -r; v(1) = -r; v(2) = r; cube.push_back(v);
	v(0) = -r; v(1) = r; v(2) = -r; cube.push_back(v);
	v(0) = r; v(1) = r; v(2) = -r; cube.push_back(v);
	v(0) = r; v(1) = r; v(2) = r; cube.push_back(v);
	v(0) = -r; v(1) = r; v(2) = r; cube.push_back(v);
}

void drawColorCube(float r,float R, float G, float B) {

	glColor3f(R,G,B);
	glBegin(GL_QUADS);
	// +X face
	//glColor3f(1,0,0);
	glNormal3f(1,0,0);
	glVertex3f(r,r,-r);
	glVertex3f(r,r,r);
	glVertex3f(r,-r,r);
	glVertex3f(r,-r,-r);
	// -X face
	//glColor3f(1,0,0);
	glNormal3f(-1,0,0);
	glVertex3f(-r,r,r);
	glVertex3f(-r,r,-r);
	glVertex3f(-r,-r,-r);
	glVertex3f(-r,-r,r);
	// +Y face
	//glColor3f(1,0,0);
	glNormal3f(0,1,0);
	glVertex3f(-r,r,r);
	glVertex3f(r,r,r);
	glVertex3f(r,r,-r);
	glVertex3f(-r,r,-r);
	// -Y face
	//glColor3f(1,0,0);
	glNormal3f(0,-1,0);
	glVertex3f(-r,-r,r);
	glVertex3f(-r,-r,-r);
	glVertex3f(r,-r,-r);
	glVertex3f(r,-r,r);
	// +Z face
	//glColor3f(1,0,0);
	glNormal3f(0,0,1);
	glVertex3f(-r,r,r);
	glVertex3f(-r,-r,r);
	glVertex3f(r,-r,r);
	glVertex3f(r,r,r);
	// -Z face
	//glColor3f(1,0,0);
	glNormal3f(0,0,-1);
	glVertex3f(-r,r,-r);
	glVertex3f(r,r,-r);
	glVertex3f(r,-r,-r);
	glVertex3f(-r,-r,-r);
	glEnd();
}

void draw_grid(int size)
{
	glColor3f(0.32,0.32,0.32);
	float dx = 1.0;
	float x = -size;
	float z = size;
	glBegin(GL_LINES);
	for (int i = -10; i < 10; ++i)
	{
		glVertex3f(x, 0, z);
		glVertex3f(x, 0, -z);
		x += dx;
	}

	for (int i = -10; i <= 10; ++i)
	{
		glVertex3f(-x, 0, z);
		glVertex3f(x, 0, z);
		z -= dx;
	}

		glEnd();
}

void draw_axis(int size)
{
	//Draw coordinate axis
	glViewport(0,0,100,100);

	glBegin(GL_LINES);
	glColor3f(1.0,0.0,0.0);
	glVertex3f(0,0,0);
	glVertex3f(size,0,0);

	glColor3f(0.0,1.0,0.0);
	glVertex3f(0,0,0);
	glVertex3f(0.0,size,0.0);

	glColor3f(0.0,0.0,1.0);
	glVertex3f(0,0,0);
	glVertex3f(0.0,0.0,size);
	glEnd();



	//Draw coordinate axis
	glViewport(0,0,winw,winh);
}

//----------------------------------------------------------------------------//
// Draws all content
//----------------------------------------------------------------------------//

static float g_lightPos[4] = { 10, 10, +20, 1 };  // Position of light
void OpenGl_drawAndUpdate(bool &running, std::vector<RigidBody> &rb)
{
	if(wireframe)
	{
		glDisable(GL_DEPTH_TEST);
		glDisable(GL_CULL_FACE);
		glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
		glDisable(GL_LIGHTING);
	}
	else
	{
		glEnable(GL_LIGHTING);
		glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
		glLightfv(GL_LIGHT0, GL_POSITION, g_lightPos);
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_CULL_FACE);
	}


	running = !glfwGetKey( GLFW_KEY_ESC ) && glfwGetWindowParam( GLFW_OPENED );

	if(!running)
		return;

	showFPS(winw, winh, zoom);
	Resize(); //Update viewport and projection matrix

	//Clear the buffer color and depth
	glClearColor(0.25f,0.25f,0.25f,1.0f);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);



	//Disables vsync
	//wglSwapIntervalEXT(0);

	// Set up viewing transformation, looking down -Z axis
	glLoadIdentity();
	glTranslatef(0,0,zoom);
	gluLookAt(0, 5, 10, 0, 0, 1, 0, 1, 0);

	// Set up the stationary light

	glTranslatef(posDx,posDy,0);
	glRotatef(-rotDx, 1,0,0);
	glRotatef(-rotDy, 0,1,0);


	// Render the scene
	glPointSize(4);
	glBegin(GL_POINTS);
		glVertex3f(0,0,0);

		// Draw contact points and contact vector
		if(conts != NULL)
			for(int i = 0; i < conts->size(); ++i)
			{
				glColor3f(1.0,1.0,0.0);  //Yellow
				glVertex3f(conts->at(i).P(0),conts->at(i).P(1), conts->at(i).P(2));
				glColor3f(1.0,1.0,1.0);
				glVertex3f(conts->at(i).pOnA(0),conts->at(i).pOnA(1), conts->at(i).pOnA(2));
			}
	glEnd();

	//Draws the collision normal at the point of collision
	glBegin(GL_LINES);
		if(conts != NULL)
		for(int i = 0; i < conts->size(); ++i)
		{
			glColor3f(1.0,1.0,1.0);
			glVertex3f(conts->at(i).P(0),conts->at(i).P(1), conts->at(i).P(2));
			glVertex3f(conts->at(i).N(0)+conts->at(i).P(0),conts->at(i).N(1)+conts->at(i).P(1), conts->at(i).N(2)+conts->at(i).P(2));
		}
	glEnd();
	
	//Draw normalized angular velocity vector
	/*
	glBegin(GL_LINES);		
		for(int i = 0; i < rb.size(); ++i)
		{
			glColor3f(1.0,0.0,1.0);
			glVertex3f(rb[i].X(0),rb[i].X(1),rb[i].X(2));
			arma::vec w = rb[i].W;
			//w = w/norm(w,2);
			glVertex3f(w(0)+rb[i].X(0),w(1)+rb[i].X(1), w(2)+rb[i].X(2));


			glColor3f(1.0,1.0,0.0);
			glVertex3f(rb[i].X(0),rb[i].X(1),rb[i].X(2));
			w = rb[i].V;
			//w = w/norm(w,2);
			glVertex3f(w(0)+rb[i].X(0),w(1)+rb[i].X(1), w(2)+rb[i].X(2));
		}
	glEnd();
	*/

	if(showgrid)
		draw_grid(10.0);

	draw_axis(3.0);

	//for (int i = rb.size()-1; i >= 0; --i)
	for (int i = 0; i < rb.size(); ++i)
	{		
		glPushMatrix();
			glTranslatef(rb[i].X(0),rb[i].X(1), rb[i].X(2));
			glPushMatrix();
			glMultMatrixd(rb[i].R.memptr());
			glScalef(rb[i].sx, rb[i].sy, rb[i].sz);

			//if(rb[i].isColliding)
			//	drawColorCube(0.5,1.0,0.0,0.0);
			//else
				drawColorCube(0.5,0.6,0.6,0.6);

			glPopMatrix();
		glPopMatrix();
		rb[i].isColliding = false;
	}

	glfwSwapBuffers();

}


//----------------------------------------------------------------------------//
// GLFW Keyboard callback
//----------------------------------------------------------------------------//
void GLFWCALL KeyboardFunc( int key, int action )
{
	if(key == 'G')
	{
		if(action == GLFW_PRESS)
			gjkdraw = !gjkdraw;
	}

	if(key == 'M')
	{
		if(action == GLFW_PRESS)
			skip = !skip;

		if(skip)
			std::cout << "skip = true";
		else
			std::cout << "skip = false";
	}

	
	if(key == 'P')
	{
		if(action == GLFW_PRESS)
			play = !play;
	}

	if(key == 'G')
	{
		if(action == GLFW_PRESS)
			showgrid = !showgrid;

	}

	if(key == 'W')
	{
		if(action == GLFW_PRESS)
			wireframe = !wireframe;

	}

	if(key == 'R' && action == GLFW_PRESS)
		reset = true;

	if(key == 'S' && action == GLFW_PRESS)
		step = true;
	if(key == 'S' && action == GLFW_RELEASE)
		step = false;

	if(key == 'C' && action == GLFW_PRESS) //Cam reset
	{
		rotDy = 0;
		rotDx = 0;
		posDx = 0;
		posDy = 0;
		zoom = 0;

	}

	if(key == GLFW_KEY_ESC)
	{
		running = false;
	}

	if (key == GLFW_KEY_UP && action == GLFW_PRESS)
		up_is_down = true;
	if(key == GLFW_KEY_UP && action == GLFW_RELEASE)
		up_is_down = false;

	if(key == GLFW_KEY_DOWN && action == GLFW_PRESS)
		down_is_down = true;
	if(key == GLFW_KEY_DOWN  && action == GLFW_RELEASE)
		down_is_down = false;


	if (key == GLFW_KEY_LSHIFT && action == GLFW_PRESS)
		shiftdown = true;
	if(key == GLFW_KEY_LSHIFT && action == GLFW_RELEASE)
		shiftdown = false;



}


//----------------------------------------------------------------------------//
// GLFW MouseButton callback
//----------------------------------------------------------------------------//
void GLFWCALL MouseButtonFunc( int button, int action )
{

}

//----------------------------------------------------------------------------//
// GLFW MouseWheel callback
//----------------------------------------------------------------------------//
void GLFWCALL MouseWheelFunc( int pos )
{
	if(shiftdown)
		zoom += (pos - lastwheelpos) *0.2;
	else
		zoom += (pos - lastwheelpos) *1.2;

	lastwheelpos = pos;
}

//----------------------------------------------------------------------------//
// GLFW MouseButton callback
//----------------------------------------------------------------------------//
void GLFWCALL MousePosFunc( int x, int y )
{
	if(glfwGetKey(GLFW_KEY_LCTRL) == GLFW_PRESS && glfwGetMouseButton(GLFW_MOUSE_BUTTON_1) == GLFW_PRESS)
	{
		rotDy += (lastmousex - x) * 0.5;
		rotDx += (lastmousey - y) * 0.5;
	}
	else if(glfwGetMouseButton(GLFW_MOUSE_BUTTON_1) == GLFW_PRESS)
	{
		posDx -= (lastmousex - x) * 0.05;
		posDy += (lastmousey - y) * 0.05;
	}
	lastmousex = x;
	lastmousey = y;
}

//----------------------------------------------------------------------------//
// Creates and sets up a window
//----------------------------------------------------------------------------//
void OpenGl_initViewer(int width_, int height_)
{
	winw = width_;
	winh = height_;

	// Initialize GLFW
	if( !glfwInit() )
	{
		exit( EXIT_FAILURE );
	}

	// Open an OpenGL window using glfw
	if( !glfwOpenWindow( winw,winh, 8,8,8,8,24,0, GLFW_WINDOW ) )
	{
		glfwTerminate();
		exit( EXIT_FAILURE );
	}

	//Init glew!
	GLenum err = glewInit();
	if (GLEW_OK != err)
	{
		/* Problem: glewInit failed, something is seriously wrong. */
		fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
	}
	fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));

	//Setup callbacks
	glfwSetKeyCallback(KeyboardFunc);
	glfwSetMouseButtonCallback(MouseButtonFunc);
	glfwSetMousePosCallback(MousePosFunc);
	glfwSetMouseWheelCallback(MouseWheelFunc);

	glfwGetMousePos(&lastmousex, &lastmousey);
	lastwheelpos = glfwGetMouseWheel();

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_CULL_FACE);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_COLOR_MATERIAL);
	glHint(GL_CLIP_VOLUME_CLIPPING_HINT_EXT,GL_FASTEST);



	int depth;
	glGetIntegerv(GL_DEPTH_BITS, &depth);
	std::cout << "GL_DEPTH_BITS: " << depth << "\n";

}













#endif
