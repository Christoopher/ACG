#ifndef COLLISION_RESPONSE_H
#define COLLISION_RESPONSE_H
#include <cstdlib>
#include <cmath>
//3rd party
#ifdef _WIN32
#include "armadillo.h"
#else
#include <armadillo>
#endif
void
createA(std::vector<Contact>& contactArray,arma::mat& A)
{
    for(unsigned int i=0; i < contactArray.size(); ++i)
    {
        Contact ci = contactArray.at(i);

        arma::vec riAlphaCrossN = arma::cross(ci.P-ci.A->X,ci.N);
        arma::vec riBetaCrossN = arma::cross(ci.P-ci.B->X,ci.N);
        /*
        std::cout << "ci.N: " << ci.N << std::endl;
        std::cout << "ci.P: " << ci.P << std::endl;
        std::cout << "ci.A->X: " << ci.A->X << std::endl;
        std::cout << "ci.B->X: " << ci.B->X << std::endl;*/
        std::cout << "ci.P-ci.A->X: " << ci.P-ci.A->X << std::endl;
        std::cout << "ci.P-ci.B->X: " << ci.P-ci.B->X << std::endl;
        /*
        std::cout << "riAlphaCrossN: " << riAlphaCrossN << std::endl;
        std::cout << "riBetaCrossN: " << riBetaCrossN << std::endl;
        */
        for(unsigned int j=0; j < contactArray.size(); ++j)
        {
            Contact cj = contactArray.at(j);

            arma::vec rjAlphaCrossN = arma::cross(cj.P-cj.A->X,cj.N);
            arma::vec rjBetaCrossN = arma::cross(cj.P-cj.B->X,cj.N);

            float dotCommon = arma::dot(ci.N,cj.N);
            A(i,j)=0.0;
            float a = ci.A->inv_mass * dotCommon;

            a += arma::dot(riAlphaCrossN,ci.A->inv_inertia * rjAlphaCrossN);

            A(i,j) +=(ci.A==cj.A)*a;
            A(i,j) -=(ci.A==cj.B)*a;

            //För de sista två ifsen
            float b = ci.B->inv_mass * dotCommon + arma::dot(riBetaCrossN,ci.B->inv_inertia * rjBetaCrossN);
            A(i,j) -= (ci.B==cj.A)*b;
            A(i,j) += (ci.B==cj.B)*b;

            /*
            std::cout << "cj.N: " << cj.N << std::endl;
            std::cout << "cj.P: " << cj.P << std::endl;
            std::cout << "cj.A->X: " << cj.A->X << std::endl;
            std::cout << "cj.B->X: " << cj.B->X << std::endl;
            std::cout << "cj.P-cj.A->X: " << cj.P-ci.A->X << std::endl;
            */
            /*
            std::cout << "rjAlphaCrossN: " << rjAlphaCrossN << std::endl;
            std::cout << "rjBetaCrossN: " << rjBetaCrossN << std::endl;
            std::cout << "dotCommon: " << dotCommon << std::endl;
            std::cout << "ci.A->inv_mass*dotCommon: " << ci.A->inv_mass * dotCommon << std::endl;
            std::cout << "+=dot(riAlphaCrossN,ci.A->inv_inertia * rjAlphaCrossN): " << a << std::endl;
            std::cout << "i: " << i << "j: "<< j << std::endl;
            std::cout << "(ci.A==cj.A)*a: " << (ci.A==cj.A)*a << std::endl;
            std::cout << "(ci.A==cj.B)*a: " << (ci.A==cj.B)*a << std::endl;
            std::cout << "ci.B->inv_mass: " << ci.B->inv_mass << std::endl;
            std::cout << "dot(riBetaCrossN,ci.B->inv_inertia * rjBetaCrossN): " << dot(riBetaCrossN,ci.B->inv_inertia * rjBetaCrossN) << std::endl;
            std::cout << "b: " << b << std::endl;
            std::cout << "(ci.B==cj.A)*b: " << (ci.B==cj.A)*b << std::endl;
            std::cout << "(ci.B==cj.B)*b: " << (ci.B==cj.B)*b << std::endl;
            std::cout << "ci.A->inv_inertia" << ci.A->inv_inertia << std::endl;
            std::cout << "ci.B->inv_inertia" << ci.B->inv_inertia << std::endl;
            */
        }
    }

}
void
computePreImpulseVelocity (std::vector<Contact>& contactArray,
arma::vec& ddot)
{
    arma::vec3 rAi,rBi;
    rAi.set_size(3,1);rBi.set_size(3,1);
    Contact ci;

    for (unsigned int i = 0; i < contactArray.size(); ++i)
    {
        ci = contactArray.at(i);


        rAi = ci.P - ci.A->X;
        rBi  = ci.P - ci.B->X;

        //Eq. 6.98
        ddot(i) = arma::dot(ci.N,(ci.A->V + arma::cross(ci.A->W,rAi)) - (ci.B->V + arma::cross(ci.B->W,rBi)));
    }


    std::cout << "ci.N: " << ci.N << std::endl;
    std::cout << "ci.A->V: " << ci.A->V << std::endl;
    std::cout << "ci.A->W: " << ci.A->W << std::endl;
    std::cout << "ci.B->V: " << ci.B->V << std::endl;
    std::cout << "ci.B->W: " << ci.B->W << std::endl;

    std::cout << "ci.A->V + + arma::cross(ci.A->W,rAi): " << ci.A->V + arma::cross(ci.A->W,rAi) << std::endl;
    std::cout << "(ci.B->V + arma::cross(ci.B->W,rBi)): " << (ci.B->V + arma::cross(ci.B->W,rBi)) << std::endl;
    std::cout << "ddot: " << ddot << std::endl;
   /* */
}

void
lemke(arma::mat& M,arma::vec& q, arma::vec& z)
{
    int n = q.n_rows;
    float zErrorTol = 1e-5;
    float pivTol    = 1e-8;
    int maxIter = std::min(1000,25*n);
    unsigned int err = 0;

    //Check trivial solution
    bool trivial = true;
    for(int i = 0;i<n;++i)
    {
        if(q(i)<0)
        {
            trivial = false;
            break;
        }
    }
    if(trivial)
    {
        z = arma::zeros<arma::vec>(n);
        return;
    }
    z = arma::zeros<arma::vec>(2*n);


    //Determine initial basis
    arma::mat B = -(arma::diagmat(arma::eye(n,n)));
    arma::uvec bas; bas.set_size(n);
      for(int i=0;i<n;++i)
        bas.at(i)= n + i;

    //Determine initial values

      arma::vec x; x.set_size(n,1);

       //   x = -arma::solve(B,q); Ändrat här själv
          x = q;
          /*
          std::cout << "q: " << q << std::endl;
          std::cout << "B: " << B << std::endl;
          std::cout << "-arma::solve(B,q): " << x << std::endl;
          */
    //Check if initial basis provides solution, already done now? since x = q ?!
    /*
    trivial = true;
    for(unsigned int i = 0;i<n;++i)
    {
        if(x(i)<0);
        {
            trivial = false;
            break;
        }
    }
    if(trivial)
    {
        z.rows(n,2*n-1) = x;
        z = z.rows(0,n-1);
        //z.rows(n+1,2*n) = x;
        //z = z.rows(1,n); //wierd? z(bas)=x; z=z(1:n);
        return;
    }
    */
   // arma::uword t = 2*n+1;
    arma::uword t = 2*n; //Correct for armadillo?
    arma::uword entering = t;

    //Determine initial leaving variable
    arma::uword lvindex;
    double tval = x.min(lvindex);
    //std::cout << "x: " << x << std::endl;
    //std::cout << "tval: " << tval << std::endl;
    //std::cout << "lvindex: " << lvindex << std::endl;
    //std::cout << "bas: " << bas << std::endl;
    //std::cout << "t: " << t << std::endl;
    tval = -tval;
    arma::uword leaving = bas(lvindex);
    bas(lvindex)=t;
    //std::cout << "leaving: " << leaving << std::endl;
    //std::cout << "bas(lvindex)=t: " << bas << std::endl;
   // std::cout << "x before addition of scalar: " << x << std::endl;
    x=x+tval; //vector + double samma som i matlab? JA
    //std::cout << "x after addition of scalar: " << x << std::endl;
    x(lvindex)=tval;
    //std::cout << "x(lvindex)=tval: " << x << std::endl;
    //std::cout << "B before changing B.col(lvindex): " << B << std::endl;
    //B.col(lvindex)=-B*arma::ones<arma::vec>(n,1);
    B.col(lvindex) = arma::ones<arma::vec>(n,1);
    //std::cout << "B after changing B.col(lvindex): " << B << std::endl;
    arma::vec Be; Be.set_size(n,1); Be.zeros();
    //Main iterations!
    arma::uvec j;j.set_size(n);j.zeros(n);
    int iter;
    for(iter=0;iter<maxIter;++iter)
    {
        //Check if done
        if(leaving == t) break;
        else if(leaving < n) //leaving <= n
        {
           // std::cout << "bas i loop: " << bas << std::endl;
          //  std::cout << "leaving: " << leaving << std::endl;
           // std::cout << "n: " << n << std::endl;
            entering = n+leaving;
            Be.zeros();
            Be(leaving) = -1.0;//sparse(leaving,1,-1.0,n,1);
          //  std::cout << "Be: " << Be << std::endl;
        }
        else
        {
            entering = leaving-n;
            if(entering >= n)
               std::cout << "entering index out of bounds: " << entering << "n: " << n << std::endl;// LOG("entering index out of bounds: " << entering << "n: " << n);
            Be = M.col(entering);

        }
        arma::vec d;

        d = arma::solve(B,Be);

        /*
        std::cout << "M: " << M << std::endl;
        std::cout << "entering = leaving -n: " << entering << std::endl;
        std::cout << "Be = M.col(entering): " << Be << std::endl;
        std::cout << "B: " << B << std::endl;
        std::cout << "d=solve(B,Be): " << d << std::endl;
        */


        //Find new leaving variable
      //  std::cout << "d: " << d << std::endl;
      //  std::cout << "tol: " << pivTol << std::endl;
        j = arma::find(d>pivTol);
       // std::cout << "j = arma::find(d>pivTol): " << j << std::endl;
        if(j.empty())
        {
            LOG("Inne i unbounded ray!");
          //  LOG("j = arma::find(d>pivTol): " << j);
           // LOG("d: " << d);
            err=2;
            break;
        }

        double theta = arma::min((x.elem(j)+zErrorTol)/d.elem(j)); //double + vector samma som matlab?
        //std::cout << "x.elem(j)" << x.elem(j) << std::endl;
        //std::cout << "x.elem(j)+zErrorTol: " << x.elem(j)+zErrorTol << std::endl;
        //std::cout << "d.elem(j)" << d.elem(j) << std::endl;
        //std::cout << "(x.elem(j)+zErrorTol)/d.elem(j)" << (x.elem(j)+zErrorTol)/d.elem(j) << std::endl;
        //std::cout << "min av ovan = theta: " << theta << std::endl;


        arma::vec divVec = x.elem(j)/d.elem(j);
        //std::cout << "divVex: " << divVec << std::endl;
        arma::uvec tmpInd = arma::find(divVec <=theta);
        //std::cout << "tmpInd=arma::find(divVec <=theta): " << tmpInd << std::endl;
        j=j.elem(tmpInd);
        //std::cout << "j=j.elem(tmpInd): " << j << std::endl;

        arma::uvec sbo = arma::find(bas.elem(j)==t);
        if(!sbo.empty())
            lvindex = sbo(0);
        //std::cout << "sbo.size(): " << sbo.size() << std::endl;
        //std::cout << "!sbo.empty(): " << !sbo.empty() << std::endl;
        if(sbo.size()>1)
        {
            LOG("SBO.size() > 1");
            abort();
        }

              if(!sbo.empty())
                  lvindex=j(lvindex);
                else
                {
                    theta = arma::max(d.elem(j));
                    //std::cout << "d.elem(j): " << d.elem(j) << std::endl;
                    arma::uvec index = arma::find(d.elem(j)==theta);
                    double ed = index.size();
                    arma::vec edd = arma::randu<arma::vec>(1);
                    double random = edd[0];
                    double bajs = ed*random;

                    //std::cout << "ed = index.size(): " << ed << std::endl;
                    //std::cout << "random: " << random << std::endl;
                    //std::cout << "bajs= ed*random: " << bajs << std::endl;
                    arma::uword tmp = (arma::uword)ceil(bajs)-1;
                    //std::cout << "ciel(bajs): " << tmp << std::endl;
                    lvindex = j(tmp);
                    //std::cout << "lvindex = j(tmp): " << lvindex << std::endl;
                }
                leaving = bas(lvindex);
                //std::cout << "bas: " << bas << std::endl;
                //std::cout << "leaving = bas(lvindex): " << leaving << std::endl;
                //Perform pivot
                double ratio = x(lvindex)/d(lvindex);
                //std::cout << "x(lvindex): " << x(lvindex) << std::endl;
                //std::cout << "d(lvindex): " << d(lvindex) << std::endl;
                //std::cout << "ratio: " << ratio << std::endl;
                //std::cout << "d: " << d << std::endl;
                //std::cout << "x before: " << x << std::endl;
                x = x - (ratio*d);
                //std::cout << "x after: " << x << std::endl;


                x(lvindex) = ratio;
                //std::cout << "x(lvindex) = ratio, visar x: " << x << std::endl;
                //std::cout << "B: " << B << std::endl;
                //std::cout << "Be: " << Be << std::endl;
                B.col(lvindex) = Be;
               // std::cout << "B after B.col(lvindex) =Be: " << B << std::endl;
                bas(lvindex) = entering;
                //std::cout << "bas after bas(lvindex) = entering: " << bas << std::endl;
    }
    /*
    std::cout << "After loop" << std::endl;
    std::cout << "Bas: "<< bas << std::endl;
    std::cout << "leaving: " << leaving << std::endl;
    std::cout << "t: " << t << std::endl;
    std::cout << "iter: " << iter << std::endl;
    */
    if(iter>=maxIter && leaving !=t)
        err=1;

   // z.rows(n+1,2*n) = x;
   // LOG("x: " << x);
  //  LOG("z before last assignment: " << z);
    z.elem(bas) = x;
   // LOG("z in middle of last assignment: " << z);
    z = z.rows(0,n-1);
   // z.rows(n,2*n-1) = x;
   // z = z.rows(0,n-1);
  //  LOG("z after last assignment: " << z);
   // z = z.rows(1,n); //wierd? z(bas)=x; z=z(1:n);

    if(err!=0)
    {
        if(err==1)
        {
            LOG("Exceeded max iterations!");
        }
        else if(err==2)
        {
            LOG("Unbounded ray! Should not happen in our case!");
        }
    }
}

void
minimize(arma::mat& A,arma::vec& ddot,arma::vec& f)
{
    //|Af +b|^2
    //Skriv om till LCP: Mz + q = w;
    //Skapa M skapa q.
    // M och q beror av A,b,c och vi har A så skapa b och c
    int n = ddot.n_rows;

    arma::vec b = arma::zeros<arma::vec>(n,1);
    arma::vec c = arma::zeros<arma::vec>(n,1);
   for(int i = 0;i < n;++i)
   {
      //Här skapas b, se boksidan 493.
       if(ddot(i) < 0.0)
           b(i) = 2.0*ddot(i); //else behövs inte ty alla b element är 0 från början
       //Här skapas c, se boksidan 496.
       c(i) = fabs(ddot(i));
   }

   //Yey nu har vi A, b och c så nu kan vi skapa M och q.
   arma::mat M = arma::zeros<arma::mat>(3*n,3*n);
   arma::vec q = arma::zeros<arma::vec>(3*n,1);
   //.submat( first_row, first_col, last_row, last_col )

  //Här skapas M, se boksidan 844.
   M.submat(0,0  ,n-1, n-1) = 2.0*arma::trans(A)*A;
   M.submat(0,n  ,n-1, 2*n-1) = -A;
   M.submat(0,2*n,n-1, 3*n-1) = A;

   M.submat(n,0,2*n-1,  n-1) = A;
   M.submat(2*n,0,3*n-1,n-1) = -A;

   //Här skapas q, se boksidan 844.
   q.subvec(0,n-1) = 2.0*arma::trans(A)*b;
   q.subvec(n,2*n-1) = b;
   q.subvec(2*n,3*n-1) = c-b;
/*
   std::cout << "b: " << b << std::endl;
   std::cout << "c: " << c << std::endl;
   std::cout << "M: " << M << std::endl;
   std::cout << "q: " << q << std::endl;
*/
   //Se Er1ks mobilbild också!

   arma::vec z = arma::zeros<arma::vec>(3*n,1);
   arma::vec w = arma::zeros<arma::vec>(3*n,1);
   //LCP STUFF GOES HERE
   //---------------------------
   lemke(M,q,z);

   //w = M*z+q;

   f = z.rows(0,n-1);
   /*
   std::cout << "w: " << w << std::endl;
   std::cout << "z: " << z << std::endl;
   std::cout << "w%z: " << w%z << std::endl;
   std::cout << "f: " << f << std::endl;
   */
   //---------------------------

}

void
applyImpulse (std::vector<Contact>& contactArray,
              arma::vec& f)
{
    Contact ci;
    arma::vec impulse;
    impulse.set_size(3,1);
    for (unsigned int i = 0; i < contactArray.size(); ++i)
    {
        ci = contactArray.at(i);
        impulse = f(i) * ci.N;
        //Optimization: Add up all impulses to P and L in separate pass
        //and then assign the new velocities?
        //Parallelization: Find a way to skip +=
        ci.A->P += impulse;
        ci.B->P -= impulse;

        ci.A->L += arma::cross(ci.P - ci.A->X,impulse);
        ci.B->L -= arma::cross(ci.P - ci.B->X,impulse);

        ci.A->V = ci.A->inv_mass * ci.A->P;
        ci.B->V = ci.B->inv_mass * ci.B->P;

        ci.A->W = ci.A->inv_inertia * ci.A->L;
        ci.B->W = ci.B->inv_inertia * ci.B->L;
    }
}

void
computeRestingVector(std::vector<Contact>& contactArray,arma::vec& b)
{
//LOG("Debug in computeRestingVector.");
    for (unsigned int i = 0; i < contactArray.size(); ++i)
    {

        Contact ci = contactArray.at(i);
        // body A terms
        arma::vec rAi = ci.P - ci.A->X;
        arma::vec wAxrAi = arma::cross(ci.A->W,rAi);
        arma::vec At1 = ci.A->inv_mass * ci.A->m_externalForce;
        arma::vec At2 = arma::cross(ci.A->inv_inertia * (ci.A->m_externalTorque + arma::cross(ci.A->L,ci.A->W)),rAi);
        arma::vec At3 = arma::cross(ci.A->W,wAxrAi);
        arma::vec At4 = ci.A->V + wAxrAi;

        // body B terms
        arma::vec rBi = ci.P - ci.B->X;
        arma::vec wBxrBi = arma::cross(ci.B->W,rBi);

        arma::vec Bt1 = ci.B->inv_mass * ci.B->m_externalForce;
        arma::vec Bt2 = arma::cross(ci.B->inv_inertia * (ci.B->m_externalTorque + arma::cross(ci.B->L,ci.B->W)),rBi);
        arma::vec Bt3 = arma::cross(ci.B->W,wBxrBi);
        arma::vec Bt4 = ci.B->V + wBxrBi;

        // compute the derivative of the contact normal

        arma::vec3 Ndot;
        if (ci.isVFContact)
        {
            Ndot = arma::cross(ci.B->W,ci.N);
        }
        else
        {
            arma::vec3 EAdot = arma::cross(ci.A->W,ci.EA);
            arma::vec3 EBdot = arma::cross(ci.B->W,ci.EB);
            arma::vec3 U = arma::cross(ci.EA,EBdot) + arma::cross(EAdot,ci.EB);
            Ndot = (U - arma::dot(U,ci.N) * ci.N) / arma::norm(ci.N,2);
        }
        b(i) = arma::dot(ci.N,At1 + At2 + At3 - Bt1 - Bt2 - Bt3) + 2.0*arma::dot(Ndot,At4 - Bt4);
        /*
        LOG("rAi = ci.P - ci.A->X: " << rAi);
        LOG("wAxrAi = arma::cross(ci.A->W,rAi): " << wAxrAi);
        LOG("At1 = ci.A->inv_mass * ci.A->m_externalForce: " << At1);
        LOG("ci.A->inv_inertia" << ci.A->inv_inertia);
        LOG("ci.A->inv_inertia * (ci.A->m_externalTorque + arma::cross(ci.A->L,ci.A->W)): " << ci.A->inv_inertia * (ci.A->m_externalTorque + arma::cross(ci.A->L,ci.A->W)));
        LOG("At2 = arma::cross(ci.A->inv_inertia * (ci.A->m_externalTorque + arma::cross(ci.A->L,ci.A->W)),rAi): " << At2);
        LOG("At3 = arma::cross(ci.A->W,wAxrAi): " << At3);
        LOG("At4 = ci.A->V + wAxrAi: " << At4);

        LOG("norm At1: " << arma::norm(At1,2));
        LOG("norm At2: " << arma::norm(At2,2));
        LOG("norm At3: " << arma::norm(At3,2));
        LOG("norm At4: " << arma::norm(At4,2));

        LOG("norm Bt1: " << arma::norm(Bt1,2));
        LOG("norm Bt2: " << arma::norm(Bt2,2));
        LOG("norm Bt3: " << arma::norm(Bt3,2));
        LOG("norm Bt4: " << arma::norm(Bt4,2));

        LOG("arma::dot(ci.N,At1 + At2 + At3 - Bt1 - Bt2 - Bt3): " << arma::dot(ci.N,At1 + At2 + At3 - Bt1 - Bt2 - Bt3));
        LOG("2.0*arma::dot(Ndot,At4 - Bt4): " << 2.0*arma::dot(Ndot,At4 - Bt4));
        */

    }

             // LOG("Leaving computeRestingVector.");


}

void updateInternalValues(std::vector<Contact>& contactArray,arma::vec& g)
{
    // update internal force/torque
    for (unsigned int i = 0; i < contactArray.size(); ++i)
    {
        Contact ci = contactArray.at(i);
        arma::vec resting = g(i) * ci.N;
        ci.A->AppendInternalForce(resting);
        arma::vec temp = arma::cross(ci.P - ci.A->X,resting);
        ci.A->AppendInternalTorque(temp);
        temp = -resting;
        ci.B->AppendInternalForce(temp);
        temp = -arma::cross(ci.P - ci.B->X,resting);
        ci.B->AppendInternalTorque(temp);
    }

}

void
collisionResponse(std::vector<Contact>& contactArray)
{
    arma::mat A    = arma::zeros<arma::mat>(contactArray.size(),contactArray.size());
    arma::vec ddot = arma::zeros<arma::vec>(contactArray.size(),1);
    arma::vec f    = arma::zeros<arma::vec>(contactArray.size(),1);
    arma::vec b    = arma::zeros<arma::vec>(contactArray.size(),1); //RestingContactVector
    arma::vec g    = arma::zeros<arma::vec>(contactArray.size(),1); //RestingContactMagnitude
    arma::vec relA = arma::zeros<arma::vec>(contactArray.size(),1); //Relative acceleration
  //  A.set_size(contactArray.size(),contactArray.size());
  //  A.zeros(contactArray.size(),contactArray.size());
    createA(contactArray,A);
    LOG("A: " << A);
    computePreImpulseVelocity(contactArray,ddot);
    LOG("ddot: " << ddot);
    minimize(A,ddot,f);
    LOG("f: " << f);
    applyImpulse(contactArray,f);

    computeRestingVector(contactArray,b);
    LOG("b: " << b);
    //LCP STUFF GOES HERE
    //---------------------------
    //LCPSolver (int numEquations, double** M, double* Q, double* Z,double* W, int& status, int maxRetries = 100,double zeroTolerance = 0.0, double ratioError = 0.0);
    //LCPsolver(contactArray.size(),A,b,relAcc,g);


    lemke(A,b,relA);
    //g = A*b + relA;
    g = relA.rows(0,contactArray.size()-1);
    LOG("g: " << g);
    //---------------------------

    updateInternalValues(contactArray,g);

}


#endif
