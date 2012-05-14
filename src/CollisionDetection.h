#ifndef COLLISION_DETECTION_H
#define COLLISION_DETECTION_H

//#define VISGJK
//#define VISEPA
//#define FIX_COLLISION

#include <vector>
#include <cmath>
#include <cstdlib>


#include "OpenGLViewer.h"
#include "RigidBody.h"
#include "Contact.h"
#include "OpenMPDefines.h"
#include "UTLog.h"


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
        int maxIndex = -1;
        float maxDist = 100000;
        for(int i = 0; i < _a.getNumVrts(); ++i)
        {
            for(int j = 0; j < _b.getNumVrts(); ++j)
            {
                int offset = i+_stride*j;
                _points[offset].Ai = i;
                _points[offset].Bj = j;
                _points[offset].p = _a.getWorldVerts()[i] - _b.getWorldVerts()[j];
                float d = arma::norm(_points[offset].p,2);
                if(d < maxDist)
                {
                    maxDist = d;
                    maxIndex = offset;
                }
            }
        }
        return &_points[maxIndex];
    }

    /*  Returns the farthest away point in the direction d  */
    MPoint *
    support(arma::vec & d)
    {
        int i = _supportMax(d,_a);
        int j = _supportMin(d,_b);
        int offset = i+_stride*j;
        //If not set, calc MPoint
        if(_points[offset].Ai < 0) //-1 if not set
        {
            _points[offset].Ai = i;
            _points[offset].Bj = j;
            _points[offset].p = _a.getWorldVerts()[i] - _b.getWorldVerts()[j];
        }

        return &_points[offset];
    }

    void
    remove(int i)
    {
        std::swap(_points[i],_points.back());
        _points.pop_back();
    }

    std::vector<MPoint> _points;
private:
    //MPoints:
    RigidBody & _a,_b;
    int _stride;
    /* Find the furthest away support vertex in a body in direction d. Returns the index of the vert */
    int
    _supportMax(arma::vec & d, RigidBody & body)
    {
        float max = -100000.0;
        float dist = 0.0;
        int index = -1;
        for(int i = 0; i < body.getNumVrts(); ++i)
        {
            dist = arma::dot(d,body.getWorldVerts()[i]);
            if(dist > max)
            {
                max = dist;
                index = i;
            }
        }

        if(index < 0) //unset if -1
        {
            LOG("Could not find a MAX point in _supportMax: probably on the border");
            abort();

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
        //for(int j = i; j < _rank-1; ++j)
        std::swap(_points[i],_points[_rank-1]);
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

#ifdef VISGJK
    void
    drawGJK(bool &running, MinkowskiSet & mk, Simplex &s)
    {
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_CULL_FACE);
        glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
        glBlendFunc(GL_SRC_ALPHA,GL_ONE);


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
        //glLightfv(GL_LIGHT0, GL_POSITION, g_lightPos);

        glTranslatef(posDx,posDy,0);
        glRotatef(-rotDx, 1,0,0);
        glRotatef(-rotDy, 0,1,0);


        // Render the scene
        glPointSize(4);
        glDisable(GL_LIGHTING);
        glBegin(GL_POINTS);
        glColor3f(0.0,0.0,0.0);
        glPointSize(8);
        glVertex3f(0,0,0);
        glPointSize(4);
        // Draw contact points and contact vector
        glColor3f(0.5,0.5,0.5);
        for(int i = 0; i < mk._points.size(); ++i)
        {
            glVertex3f(mk._points[i].p(0),mk._points[i].p(1),mk._points[i].p(2));
        }
        glEnd();


        // A
        if(A != NULL)
        {
            glBegin(GL_POINTS);
            glColor3f(0.0,0.0,1.0);
            glVertex3f(A->p(0),A->p(1),A->p(2));
            glEnd();
        }

        glBegin(GL_LINES);

        //Direction
        glColor3f(1.0,0.0,1.0);
        arma::vec dir(d);
        dir = dir/norm(dir,2);

        glVertex3f(0,0,0);
        glVertex3f(dir(0),dir(1),dir(2));

        if(s.rank() == 4)
        {
            //normal
            glColor3f(1.0,0.0,0.0);  //RÖD
            dir = normal;
            dir = dir/norm(dir,2);
            dir += s.getImpl()[3]->p;

            glVertex3f(s.getImpl()[3]->p(0),s.getImpl()[3]->p(1),s.getImpl()[3]->p(2));
            glVertex3f(dir(0),dir(1),dir(2));


            //A0
            glColor3f(0.0,1.0,0.0); //GRÖN
            dir = s.getImpl()[3]->p;
            glVertex3f(0,0,0);
            glVertex3f(dir(0),dir(1),dir(2));
        }

        glEnd();


        glEnable(GL_BLEND);
        switch(s.rank())
        {
        case 2:
            {
                glBegin(GL_LINES);
                for(int i = 0; i < s.rank()-1; ++i)
                {
                    glColor3f(1.0,1.0,1.0);
                    glVertex3f(s[i].p(0),s[i].p(1),s[i].p(2));
                    glVertex3f(s[i+1].p(0),s[i+1].p(1),s[i+1].p(2));
                }
                glEnd();
                break;
            }
        case 3:
            {
                glBegin(GL_TRIANGLES);
                glColor4f(0.6,0.6,0.6,0.2);
                glVertex3f(s[0].p(0),s[0].p(1),s[0].p(2));
                glVertex3f(s[1].p(0),s[1].p(1),s[1].p(2));
                glVertex3f(s[2].p(0),s[2].p(1),s[2].p(2));
                glEnd();
                break;
            }
        case 4:
            {
                glBegin(GL_TRIANGLES);
                //ACB
                if(plane == 1)
                    glColor4f(0.2,0.9,0.2,0.2);
                else
                    glColor4f(0.6,0.6,0.6,0.2);
                glVertex3f(s[3].p(0),s[3].p(1),s[3].p(2)); //A 3
                glVertex3f(s[1].p(0),s[1].p(1),s[1].p(2)); //C 1
                glVertex3f(s[2].p(0),s[2].p(1),s[2].p(2)); //B 2

                //ABD
                if(plane == 2)
                    glColor4f(0.2,0.9,0.2,0.2);
                else
                    glColor4f(0.6,0.6,0.6,0.2);

                glVertex3f(s[3].p(0),s[3].p(1),s[3].p(2)); //A 3
                glVertex3f(s[2].p(0),s[2].p(1),s[2].p(2)); //B 2
                glVertex3f(s[0].p(0),s[0].p(1),s[0].p(2)); //D 0

                //ADC
                if(plane == 3)
                    glColor4f(0.2,0.9,0.2,0.2);
                else
                    glColor4f(0.6,0.6,0.6,0.2);

                glVertex3f(s[3].p(0),s[3].p(1),s[3].p(2)); //A 3
                glVertex3f(s[0].p(0),s[0].p(1),s[0].p(2)); //D 0
                glVertex3f(s[1].p(0),s[1].p(1),s[1].p(2)); //C 1

                //BDC
                glColor4f(0.6,0.6,0.6,0.2);
                glVertex3f(s[2].p(0),s[2].p(1),s[2].p(2)); //B 2
                glVertex3f(s[0].p(0),s[0].p(1),s[0].p(2)); //D 0
                glVertex3f(s[1].p(0),s[1].p(1),s[1].p(2)); //C 1

                glEnd();
                break;
            }
        }

        glfwSwapBuffers();
        glDisable(GL_BLEND);
    }
#endif

    void
        _run()
    {
        A = _mk.getStartingPoint();
        _s.insert(A);
        d = -A->p;
        _iterations = 0;
        while(true)
        {
            if(_iterations > 200)
            {
                LOG("Reached 200 iterations in GJK algorithm. aborting...");
                return;
            }
            A = _mk.support(d);

#ifdef VISGJK
            gjkdraw = true;
            while(gjkdraw)
            {
                drawGJK(gjkdraw,_mk,_s);
            }
#endif

            float dotprod=0.0;
            if((dotprod = arma::dot(A->p,d)) < THRESHOLD)
            {
                LOG("No intersection in gjk");
                return;
            }

            if(!_s.insert(A))
            {
                //@todo: Deal with trying to insert an already inserted point
                LOG("inserting an already existing point in GJK");
            }
            else
            {
                if(_doSimplex(d))
                {
                    _intersects = true;

                    if(!_tetraContainsOrigin())
                    {
                        LOG("GJK ended with a tetra not containing the origin");
                    }
#ifdef VISGJK
                    gjkdraw = true;
                    while(gjkdraw)
                    {
                        drawGJK(gjkdraw,_mk,_s);
                    }
#endif

                    return;
                }
            }
#ifdef VISGJK
            gjkdraw = true;
            while(gjkdraw)
            {
                drawGJK(gjkdraw,_mk,_s);
            }
#endif
            ++_iterations;
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
                    return false;
                }
                else //If origin is beyond A in the direction AB
                {
                    d = Ao;
                    _s.keep(1); //Only keep A in the simplex
                    return false;
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
                        return false;
                    }
                    else //STAR
                    {
                        if(dot(AB, Ao) > 0) //keep [B A] from [C B A]
                        {
                            d = cross(cross(AB, Ao), AB);

                            _s.remove(0);
                            return false;
                        }
                        else //Only keep [ A ] from [ C B A ]
                        {
                            d = Ao;
                            _s.keep(2);
                            return false;
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
                            return false;
                        }
                        else //Only keep [ A ] from [ C B A ]
                        {
                            d = Ao;
                            _s.keep(2);
                            return false;
                        }
                    }
                    else //d is either ABC or -ABC. The simplex is [A B C]
                    {
                        if(dot(ABC, Ao) > 0)
                            d = ABC;
                        else
                            d = -ABC;
                        return false;
                    }
                }
            }
            break;
#pragma endregion
#pragma region TETRAHEDRON
        case TETRAHEDRON: //Having S = [D C B A]
            {
                float dotprod = 0.0;
                arma::vec t1 = _s[1].p - _s[3].p;
                arma::vec t2 = _s[2].p - _s[3].p;
                arma::vec ACB = cross(_s[1].p - _s[3].p,_s[2].p - _s[3].p); //cross(AC, AB) -> triangle ACB
                if(arma::dot(ACB,_s[0].p - _s[3].p) > 0) // dot(n,AD) > 0
                {
                    ACB = -ACB;
                    LOG("flipped normal");
                }
#ifdef VISGJK
                plane = 1;
                d = ACB;
                normal = ACB;
                LOG("check if outside ACB");
                gjkdraw = true;
                while(gjkdraw)
                {
                    drawGJK(gjkdraw,_mk,_s);
                }
#endif

                if((dotprod = dot(ACB, Ao)) > 0) //If on the outside of ACB
                {
                    if(_tetraContainsOrigin())
                    {
                        LOG("Origin is contained but gjk is removing D. dotprod =" << dotprod);
                        return true;
                    }
                    //Remove D
                    _s.remove(0);

                    d = ACB;
                    return false;
                    //return _doSimplex(d);
                }
                else //Inside tetra or outside of ABD or ADC
                {
                    /*
                        The normal dir problem might be caused by numerical issues when

                        the point lies almost exactly on the plane in question
                    */
                    arma::vec t1 = _s[2].p - _s[3].p;
                    arma::vec t2 = _s[0].p - _s[3].p;
                    arma::vec ABD = cross(_s[2].p - _s[3].p, _s[0].p - _s[3].p); //cross(AB, AD) -> triangle ABD
                    if(arma::dot(ABD,_s[1].p - _s[3].p) > 0) // dot(n,AC) > 0
                    {
                        ABD = -ABD;
                        LOG("flipped normal");
                    }


#ifdef VISGJK
                    plane = 2;
                    normal = ABD;
                    d = ABD;
                    LOG("check if outside ABD");
                    gjkdraw = true;
                    while(gjkdraw)
                    {
                        drawGJK(gjkdraw,_mk,_s);
                    }
#endif
                    if((dotprod = dot(ABD, Ao)) > 0) //Outside ABD
                    {
                        if(_tetraContainsOrigin())
                        {
                            LOG("Origin is contained but gjk is removing C. dotprod =" << dotprod);

                            return true;
                        }


                        //Remove C
                        _s.remove(1);
                        d = ABD;
                        return false;
                        //return _doSimplex(d);
                    }
                    else //Outside ADC or inside tetra
                    {
                        LOG("check if outside ADC");
                        arma::vec ADC = cross(_s[0].p - _s[3].p, _s[1].p - _s[3].p); //cross(AD, AC) -> triangle ADC
                        if(arma::dot(ADC,_s[2].p - _s[3].p) > 0) // dot(n,AB) > 0
                        {
                            ADC = -ADC;
                            LOG("flipped normal");
                        }
#ifdef VISGJK
                        plane = 3;
                        d = ADC;
                        normal = ADC;
                        gjkdraw = true;
                        while(gjkdraw)
                        {
                            drawGJK(gjkdraw,_mk,_s);
                        }
#endif
                        if((dotprod = dot(ADC, Ao)) > 0) //Outside of ACD
                        {
                            if(_tetraContainsOrigin())
                            {
                                LOG("Origin is contained but gjk is removing B. dotprod =" << dotprod);

                                return true;
                            }
                            //remove B;
                            //LOG("A : " << _s[3].Ai << "," << _s[3].Bj);
                            //LOG("B : " << _s[2].Ai << "," << _s[2].Bj);
                            //LOG("C : " << _s[1].Ai << "," << _s[1].Bj);
                            //LOG("D : " << _s[0].Ai << "," << _s[0].Bj);
                            _s.remove(2);
                            d = ADC;
                            return false;
                            //return _doSimplex(d);
                        }
                        else //inside tetrahedron
                        {
                            return true; //WE ARE DONE :)
                        }
                    }
                }
            }
            break;
#pragma endregion
        default:
            LOG("in doSimplex() rank was not of [1,2,3,4] \n");
            abort();
            break;
        }

        return false;
    }

    bool
    _tetraContainsOrigin()
    {
        if(_s.getImpl()[3] == NULL)
        {
            LOG("trying to run _tetraContainsOrigin() with only 3 points");
            return false;
        }
        //Check if we trully enclosed the origin
        arma::mat44 D0,D1, D2, D3, D4;

        D0(0,0) = _s[0].p(0); D0(0,1) = _s[0].p(1); D0(0,2) = _s[0].p(2); D0(0,3) = 1;
        D0(1,0) = _s[1].p(0); D0(1,1) = _s[1].p(1); D0(1,2) = _s[1].p(2); D0(1,3) = 1;
        D0(2,0) = _s[2].p(0); D0(2,1) = _s[2].p(1); D0(2,2) = _s[2].p(2); D0(2,3) = 1;
        D0(3,0) = _s[3].p(0); D0(3,1) = _s[3].p(1); D0(3,2) = _s[3].p(2); D0(3,3) = 1;

        //BCD
        D1(0,0) = 0.0;		  D1(0,1) = 0.0;		D1(0,2) = 0.0;		  D1(0,3) = 1;
        D1(1,0) = _s[1].p(0); D1(1,1) = _s[1].p(1); D1(1,2) = _s[1].p(2); D1(1,3) = 1;
        D1(2,0) = _s[2].p(0); D1(2,1) = _s[2].p(1); D1(2,2) = _s[2].p(2); D1(2,3) = 1;
        D1(3,0) = _s[3].p(0); D1(3,1) = _s[3].p(1); D1(3,2) = _s[3].p(2); D1(3,3) = 1;

        //ACD
        D2(0,0) = _s[0].p(0); D2(0,1) = _s[0].p(1); D2(0,2) = _s[0].p(2); D2(0,3) = 1;
        D2(1,0) = 0.0;		  D2(1,1) = 0.0;		D2(1,2) = 0.0;		  D2(1,3) = 1;
        D2(2,0) = _s[2].p(0);	D2(2,1) = _s[2].p(1); D2(2,2) = _s[2].p(2); D2(2,3) = 1;
        D2(3,0) = _s[3].p(0); D2(3,1) = _s[3].p(1); D2(3,2) = _s[3].p(2); D2(3,3) = 1;

        //ABD
        D3(0,0) = _s[0].p(0); D3(0,1) = _s[0].p(1); D3(0,2) = _s[0].p(2); D3(0,3) = 1;
        D3(1,0) = _s[1].p(0); D3(1,1) = _s[1].p(1); D3(1,2) = _s[1].p(2); D3(1,3) = 1;
        D3(2,0) = 0.0;		  D3(2,1) = 0.0;		D3(2,2) = 0.0;	      D3(2,3) = 1;
        D3(3,0) = _s[3].p(0); D3(3,1) = _s[3].p(1); D3(3,2) = _s[3].p(2); D3(3,3) = 1;

        //ABC
        D4(0,0) = _s[0].p(0); D4(0,1) = _s[0].p(1); D4(0,2) = _s[0].p(2); D4(0,3) = 1;
        D4(1,0) = _s[1].p(0); D4(1,1) = _s[1].p(1); D4(1,2) = _s[1].p(2); D4(1,3) = 1;
        D4(2,0) = _s[2].p(0); D4(2,1) = _s[2].p(1); D4(2,2) = _s[2].p(2); D4(2,3) = 1;
        D4(3,0) = 0.0;		  D4(3,1) = 0.0;		D4(3,2) =0.0;		  D4(3,3) = 1;



        float detD0 = det(D0,true);
        float detD1 = det(D1,true);
        float detD2 = det(D2,true);
        float detD3 = det(D3,true);
        float detD4 = det(D4,true);




        if( (detD0 <= 0.0 && detD1 <= 0.0 && detD2 <= 0.0 && detD3 <= 0.0 && detD4 <= 0.0) ||
            (detD0 > 0.0 && detD1 > 0.0 && detD2 > 0.0 && detD3 > 0.0 && detD4 > 0.0) )
        {
            if(abs(detD0) < 10e-7) //close to zero
            {
                LOG("tetrahedron is degenerate (the points are coplanar)");
                abort();
            }
        }
        else
        {
            if(abs(detD1) < 10e-7) //close to zero
            {
                LOG("origin is on the border of plane 0");
            }
            else if(abs(detD2) < 10e-7) //close to zero
            {
                LOG("origin is on the border of plane 1");
            }
            else if(abs(detD3) < 10e-7) //close to zero
            {
                LOG("origin is on the border of plane 2");
            }
            else if(abs(detD4) < 10e-7) //close to zero
            {
                LOG("origin is on the border of plane 3");
            }
            else
            {
                return false;
            }

        }

        LOG("D1: " << detD1 << ", D2: " << detD2 << ", D3: " << detD3 << ", D4: " << detD4);
        return true;
    }

    RigidBody & _bodyA, _bodyB;
    MinkowskiSet::MPoint * A;
    Simplex _s;
    MinkowskiSet _mk;
    bool _intersects;
    const float THRESHOLD;
    int _iterations;
    arma::vec d;
#ifdef VISGJK
    arma::vec normal;
    int plane;
#endif

};

struct Plane
{
    Plane(MinkowskiSet::MPoint * A, MinkowskiSet::MPoint * B, MinkowskiSet::MPoint * C)
    {
        points[0] = A;
        points[1] = B;
        points[2] = C;

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
    MinkowskiSet::MPoint * points[3];
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

            if(_planes.size() > 300)
            {
                LOG("Added: " <<  _planes.size() << "planes in EPA algorithm. aborting...");
                return;
            }


            A = _gjk._mk.support(_planes[_cpi].p);
            double d = arma::dot(A->p, _planes[_cpi].p/arma::norm(_planes[_cpi].p,2));
            if (abs(d - _closestDist) < 0.00001)
            {
                //Finished!
                return;
            }
            else //We now have the closest plane, which is the plane from where a new tetrahedron will be "extruded" by	inserting a new point on the non-origin-side of this plane and adding the 3 new planes. A is the new point.
            {
                Plane abc(A, _planes[_cpi].points[0], _planes[_cpi].points[1]);
                Plane abd(A, _planes[_cpi].points[0], _planes[_cpi].points[2]);
                Plane acd(A, _planes[_cpi].points[1], _planes[_cpi].points[2]);
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

    //Set the bodies;
    c.A = &a;
    c.B = &b;

    //Calculate collision point;

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
    A = plane.points[0];
    B = plane.points[1];
    C = plane.points[2];

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


#if (DEBUG)
    double sumL = lA+lB+lC;
    if(sumL < 0.999 || sumL > 1.000999)
        LOG("lA
#endif

        /* Calculate world coordinates */
    arma::vec aWorld = lA*a.getWorldVerts()[A->Ai] + lB*a.getWorldVerts()[B->Ai] + lC*a.getWorldVerts()[C->Ai];
    arma::vec bWorld = lA*b.getWorldVerts()[A->Bj] + lB*b.getWorldVerts()[B->Bj] + lC*b.getWorldVerts()[C->Bj];
    arma::vec pWorld = aWorld - bWorld; /* The translation vector */
    c.P = bWorld;
    c.N = aWorld;


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



    //All points in A are the same
    if(A->Ai == B->Ai && B->Ai == C->Ai) // ===> A->Ai == B->Ai == C->Ai
    {
        //Can be vertex-face / vertex-edge / vertex-vertex
        LOG("vertex-face collision");
        c.isVFContact = true;
        //ACB plane i.e N = AC x AB
        c.N = arma::cross(b.getWorldVerts()[C->Bj] - b.getWorldVerts()[A->Bj], b.getWorldVerts()[B->Bj] - b.getWorldVerts()[A->Bj]);
        c.N = (1.0/arma::norm(c.N,2))*c.N;

        if(arma::dot(c.N,c.P-b.X) < 0.0)
            c.N = -c.N;
    }
    //All points are different in a
    else if(A->Ai != B->Ai && A->Ai != C->Ai  && B->Ai != C->Ai)
    {
        //can be face-face / face-vertex / face-edge
        LOG("face-face collision. NOT SUPPORTED YET. treating it as vertex face");
        c.isVFContact = true;

        //ACB plane i.e N = AC x AB
        c.N = arma::cross(b.getWorldVerts()[C->Bj] - b.getWorldVerts()[A->Bj], b.getWorldVerts()[B->Bj] - b.getWorldVerts()[A->Bj]);
        c.N /= arma::norm(c.N,2);
    }
    //Two points in object A are the same
    else
    {
        //can be edge-vertex / edge-face / edge-vertex
        LOG("edge-edge collision");
        //c.N = pWorld/norm(pWorld,2);
        c.isVFContact = false;

        int i1 = plane.points[0]->Ai, i2;
        if(plane.points[1]->Ai != i1)
            i2 = plane.points[1]->Ai;
        else
            i2 = plane.points[2]->Ai;

        c.EA = a.getWorldVerts()[i2] - a.getWorldVerts()[i1];

        i1 = plane.points[0]->Bj;
        if(plane.points[1]->Bj != i1)
            i2 = plane.points[1]->Bj;
        else
            i2 = plane.points[2]->Bj;
        c.EB = b.getWorldVerts()[i2] - b.getWorldVerts()[i1];

        c.N = arma::cross(c.EA,c.EB);
        c.N /= arma::norm(c.N,2);

        if(arma::dot(c.N,c.P-b.X) < 0) //c.P-b.X point into the object
            c.N *= -1;
    }

}

void
narrowPhase(RigidBody & bodyA,RigidBody & bodyB, std::vector<Contact> & contacts)
{
    LOG("runnin gjk");
    GJK gjk(bodyA,bodyB);

    if(gjk.intersects())
    {
        EPA epa(gjk);
        Contact c;

        makeContact(epa.plane(),bodyA,bodyB, c);
        contacts.push_back(c);

#ifdef FIX_COLLISION
        //Ugly way to separate objects
        if(bodyA.inv_mass == 0)
        {
            bodyB.X(0) += c.N(0);
            bodyB.X(1) += c.N(1);
            bodyB.X(2) += c.N(2);
        }
        else //Every other cases we only move bodyA
        {
            bodyA.X(0) -= c.N(0);
            bodyA.X(1) -= c.N(1);
            bodyA.X(2) -= c.N(2);
        }
#endif

        bodyA.isColliding = true;
        bodyB.isColliding = true;
        play = false;
    }
}

/***********************************************************************************************/
//									  BROAD PHASE
/***********************************************************************************************/
bool
broadPhase(RigidBody & a, RigidBody & b)
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
#pragma omp parallel for schedule(dynamic, chunksize)
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
