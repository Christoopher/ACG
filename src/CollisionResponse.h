#ifndef COLLISION_RESPONSE_H
#define COLLISION_RESPONSE_H
//#define checkResting
#include <cstdlib>
#include <cmath>
#include "UTLog.h"
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

		arma::vec riAlphaCrossN = arma::cross(ci.P-ci.A->X,ci.N); //A->P vector
		arma::vec riBetaCrossN = arma::cross(ci.P-ci.B->X,ci.N); //B->P vector
		/*
		LOG( "ci.N: " << ci.N );
		LOG( "ci.P: " << ci.P );
		LOG( "ci.A->X: " << ci.A->X );
		LOG( "ci.B->X: " << ci.B->X );
		LOG( "ci.P-ci.A->X: " << ci.P-ci.A->X );
		LOG( "ci.P-ci.B->X: " << ci.P-ci.B->X );
		/*
		LOG( "riAlphaCrossN: " << riAlphaCrossN );
		LOG( "riBetaCrossN: " << riBetaCrossN );
		*/
		for(unsigned int j=0; j < contactArray.size(); ++j)
		{
			Contact cj = contactArray.at(j);

			arma::vec rjAlphaCrossN = arma::cross(cj.P-cj.A->X,cj.N);
			arma::vec rjBetaCrossN = arma::cross(cj.P-cj.B->X,cj.N);


			A(i,j)=0.0;

			if (ci.A == cj.A)
			{
				A(i,j) += ci.A->inv_mass * arma::dot(ci.N,cj.N);
				A(i,j) += arma::dot(riAlphaCrossN,ci.A->inv_inertia * rjAlphaCrossN);
			}
			else if (ci.A == cj.B)
			{
				A(i,j) -= ci.A->inv_mass * arma::dot(ci.N,cj.N);
				A(i,j) -= arma::dot(riAlphaCrossN,ci.A->inv_inertia * rjAlphaCrossN);
			}
			if (ci.B == cj.A)
			{
				A(i,j) -= ci.B->inv_mass * arma::dot(ci.N,cj.N);
				A(i,j) -= arma::dot(riBetaCrossN,ci.B->inv_inertia * rjBetaCrossN);
			}
			else if (ci.B == cj.B)
			{
				A(i,j) += ci.B->inv_mass * arma::dot(ci.N,cj.N);
				A(i,j) += arma::dot(riBetaCrossN,ci.B->inv_inertia * rjBetaCrossN);
			}

			/*
			float dotCommon = arma::dot(ci.N,cj.N);
			float a = ci.A->inv_mass * dotCommon;

			a += arma::dot(riAlphaCrossN,ci.A->inv_inertia * rjAlphaCrossN);



			A(i,j) +=(ci.A==cj.A)*a;
			A(i,j) -=(ci.A==cj.B)*a;

			//För de sista två ifsen
			float b = ci.B->inv_mass * dotCommon + arma::dot(riBetaCrossN,ci.B->inv_inertia * rjBetaCrossN);
			A(i,j) -= (ci.B==cj.A)*b;
			A(i,j) += (ci.B==cj.B)*b;
			*/
			/*
			LOG( "cj.N: " << cj.N );
			LOG( "cj.P: " << cj.P );
			LOG( "cj.A->X: " << cj.A->X );
			LOG( "cj.B->X: " << cj.B->X );
			LOG( "cj.P-cj.A->X: " << cj.P-ci.A->X );
			*/
			/*
			LOG( "rjAlphaCrossN: " << rjAlphaCrossN );
			LOG( "rjBetaCrossN: " << rjBetaCrossN );
			LOG( "dotCommon: " << dotCommon );
			LOG( "ci.A->inv_mass*dotCommon: " << ci.A->inv_mass * dotCommon );
			LOG( "+=dot(riAlphaCrossN,ci.A->inv_inertia * rjAlphaCrossN): " << a );
			LOG( "i: " << i << "j: "<< j );
			LOG( "(ci.A==cj.A)*a: " << (ci.A==cj.A)*a );
			LOG( "(ci.A==cj.B)*a: " << (ci.A==cj.B)*a );
			LOG( "ci.B->inv_mass: " << ci.B->inv_mass );
			LOG( "dot(riBetaCrossN,ci.B->inv_inertia * rjBetaCrossN): " << dot(riBetaCrossN,ci.B->inv_inertia * rjBetaCrossN) );
			LOG( "b: " << b );
			LOG( "(ci.B==cj.A)*b: " << (ci.B==cj.A)*b );
			LOG( "(ci.B==cj.B)*b: " << (ci.B==cj.B)*b );
			LOG( "ci.A->inv_inertia" << ci.A->inv_inertia );
			LOG( "ci.B->inv_inertia" << ci.B->inv_inertia );
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


	//LOG( "ci.N: " << ci.N );
	//LOG( "ci.A->V: " << ci.A->V );
	//LOG( "ci.A->W: " << ci.A->W );
	//LOG( "ci.B->V: " << ci.B->V );
	//LOG( "ci.B->W: " << ci.B->W );

	LOG( "ci.A->V + + arma::cross(ci.A->W,rAi): " << ci.A->V + arma::cross(ci.A->W,rAi) );
	//LOG( "(ci.B->V + arma::cross(ci.B->W,rBi)): " << (ci.B->V + arma::cross(ci.B->W,rBi)) );

}

void
	lemke(arma::mat& M,arma::vec& q, arma::vec& z)
{
	unsigned int n = q.n_rows;
	float zErrorTol = 1e-5;
	float pivTol    = 1e-8;
	int maxIter = std::min(1000,25*(int)n);
	unsigned int err = 0;

	//Check trivial solution
	bool trivial = true;
	for(unsigned int i = 0;i<n;++i)
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
	for(unsigned int i=0;i<n;++i)
		bas.at(i)= n + i;

	//Determine initial values

	arma::vec x; x.set_size(n,1);

	//   x = -arma::solve(B,q); Ändrat här själv
	x = q;
	/*
	LOG( "q: " << q );
	LOG( "B: " << B );
	LOG( "-arma::solve(B,q): " << x );
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
	//LOG( "x: " << x );
	//LOG( "tval: " << tval );
	//LOG( "lvindex: " << lvindex );
	//LOG( "bas: " << bas );
	//LOG( "t: " << t );
	tval = -tval;
	arma::uword leaving = bas(lvindex);
	bas(lvindex)=t;
	//LOG( "leaving: " << leaving );
	//LOG( "bas(lvindex)=t: " << bas );
	// LOG( "x before addition of scalar: " << x );
	x=x+tval; //vector + double samma som i matlab? JA
	//LOG( "x after addition of scalar: " << x );
	x(lvindex)=tval;
	//LOG( "x(lvindex)=tval: " << x );
	//LOG( "B before changing B.col(lvindex): " << B );
	//B.col(lvindex)=-B*arma::ones<arma::vec>(n,1);
	B.col(lvindex) = arma::ones<arma::vec>(n,1);
	//LOG( "B after changing B.col(lvindex): " << B );
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
			// LOG( "bas i loop: " << bas );
			//  LOG( "leaving: " << leaving );
			// LOG( "n: " << n );
			entering = n+leaving;
			Be.zeros();
			Be(leaving) = -1.0;//sparse(leaving,1,-1.0,n,1);
			//  LOG( "Be: " << Be );
		}
		else
		{
			entering = leaving-n;
			if(entering >= n)
				LOG( "entering index out of bounds: " << entering << "n: " << n );// LOG("entering index out of bounds: " << entering << "n: " << n);
			Be = M.col(entering);

		}
		arma::vec d;

		d = arma::solve(B,Be);

		/*
		LOG( "M: " << M );
		LOG( "entering = leaving -n: " << entering );
		LOG( "Be = M.col(entering): " << Be );
		LOG( "B: " << B );
		LOG( "d=solve(B,Be): " << d );
		*/


		//Find new leaving variable
		//  LOG( "d: " << d );
		//  LOG( "tol: " << pivTol );
		j = arma::find(d>pivTol);
		// LOG( "j = arma::find(d>pivTol): " << j );
		if(j.empty())
		{
			LOG("Inne i unbounded ray!");
			//  LOG("j = arma::find(d>pivTol): " << j);
			// LOG("d: " << d);
			err=2;
			break;
		}

#if _WIN32
		arma::vec thetaVec = arma::min((x.elem(j)+zErrorTol)/d.elem(j)); //double + vector samma som matlab?
		//LOG("arma::vec thetaVec = arma::min((x.elem(j)+zErrorTol)/d.elem(j)): " << thetaVec);
		if(thetaVec.n_rows > 1)
		{
			LOG("size of arma::min((x.elem(j)+zErrorTol)/d.elem(j)) was greater than 1");
			abort();
		}
		double theta = thetaVec(0);
#else

		double theta = arma::min((x.elem(j)+zErrorTol)/d.elem(j)); 
#endif

		//LOG( "x.elem(j)" << x.elem(j) );
		//LOG( "x.elem(j)+zErrorTol: " << x.elem(j)+zErrorTol );
		//LOG( "d.elem(j)" << d.elem(j) );
		//LOG( "(x.elem(j)+zErrorTol)/d.elem(j)" << (x.elem(j)+zErrorTol)/d.elem(j) );
		//LOG( "min av ovan = theta: " << theta );


		arma::vec divVec = x.elem(j)/d.elem(j);
		//LOG( "divVex: " << divVec );
		arma::uvec tmpInd = arma::find(divVec <=theta);
		//LOG( "tmpInd=arma::find(divVec <=theta): " << tmpInd );
		j=j.elem(tmpInd);
		//LOG( "j=j.elem(tmpInd): " << j );

		arma::uvec sbo = arma::find(bas.elem(j)==t);
		if(!sbo.empty())
			lvindex = sbo(0);
		//LOG( "sbo.size(): " << sbo.size() );
		//LOG( "!sbo.empty(): " << !sbo.empty() );
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
			//LOG( "d.elem(j): " << d.elem(j) );
			arma::uvec index = arma::find(d.elem(j)==theta);
			double ed = index.size();
			arma::vec edd = arma::randu<arma::vec>(1);
			double random = edd[0];
			double random2 = ed*random;

			//LOG( "ed = index.size(): " << ed );
			//LOG( "random: " << random );
			//LOG( "random2= ed*random: " << random2 );
			arma::uword tmp = (arma::uword)random2;
			//LOG( "ciel(random2): " << tmp );
			lvindex = j(tmp);
			//LOG( "lvindex = j(tmp): " << lvindex );
		}
		leaving = bas(lvindex);
		//LOG( "bas: " << bas );
		//LOG( "leaving = bas(lvindex): " << leaving );
		//Perform pivot
		double ratio = x(lvindex)/d(lvindex);
		//LOG( "x(lvindex): " << x(lvindex) );
		//LOG( "d(lvindex): " << d(lvindex) );
		//LOG( "ratio: " << ratio );
		//LOG( "d: " << d );
		//LOG( "x before: " << x );
		x = x - (ratio*d);
		//LOG( "x after: " << x );


		x(lvindex) = ratio;
		//LOG( "x(lvindex) = ratio, visar x: " << x );
		//LOG( "B: " << B );
		//LOG( "Be: " << Be );
		B.col(lvindex) = Be;
		// LOG( "B after B.col(lvindex) =Be: " << B );
		bas(lvindex) = entering;
		//LOG( "bas after bas(lvindex) = entering: " << bas );
	}
	/*
	LOG( "After loop" );
	LOG( "Bas: "<< bas );
	LOG( "leaving: " << leaving );
	LOG( "t: " << t );
	LOG( "iter: " << iter );
	*/
	if(iter>=maxIter && leaving !=t)
		err=1;

	if(err!=0)
	{
		if(err==1)
		{
			LOG("Exceeded max iterations!");
		}
		else if(err==2)
		{
			LOG("Unbounded ray! Should not happen in our case!");
			return;
		}
	}
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
	double epsilon = 0.4;
	for(int i = 0;i < n;++i)
	{
		//Här skapas b, se boksidan 493.
		if(ddot(i) < 0.0)
			b(i) = (1.0+epsilon)*ddot(i);//  b(i) = ddot(i) + epsilon*ddot(i);
		// b(i) = 1.6*ddot(i); //else behövs inte ty alla b element är 0 från början
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
	LOG( "b: " << b );
	LOG( "c: " << c );
	LOG( "M: " << M );
	LOG( "q: " << q );
	*/
	//Se Er1ks mobilbild också!

	arma::vec z = arma::zeros<arma::vec>(3*n,1);
	arma::vec w = arma::zeros<arma::vec>(3*n,1);
	//LCP STUFF GOES HERE
	//---------------------------
	lemke(M,q,z);

	//w = M*z+q;

	f = z.rows(0,n-1);

	bool print = false;
	for(int i =0; i < f.n_rows;++i)
	{
		if(f(i) > 0)
		{
			print = true;
			break;
		}
	}
	if(print)
	{
		LOG("f: " << f);
		LOG("b: " << b);
		LOG("Af+b: " << A*f+b);
		LOG("Af+ddot: " << A*f+ddot);
		LOG("c: " << c);
	}
	//DEBUG
	//-------------------------------------------------------------------------
	arma::vec tmp = A*f + b;
	for(unsigned int i = 0; i < f.n_rows; ++i)
	{
		if(tmp(i) > c(i))
		{
			LOG( "Af+b > c. Något fel i lemke? Af+b: " << tmp << "c: " << c );
			//abort();
		}
		else if(tmp(i) < 0.0)
		{
			LOG( "Af+b < 0. Något fel i lemke? Af+b: " << tmp );

		}
		else if(f(i) < 0.0)
		{
			LOG( "f < 0.0. Något fel i lemke? f: " << f );
			abort();
		}
	}
	//-------------------------------------------------------------------------
	/*
	LOG( "w: " << w );
	LOG( "z: " << z );
	LOG( "w%z: " << w%z );
	LOG( "f: " << f );
	*/
	//---------------------------

}

void
	applyImpulse (std::vector<Contact>& contactArray,arma::vec& f)
{
	Contact ci;
	arma::vec impulse;
	impulse.set_size(3,1);
	for (unsigned int i = 0; i < contactArray.size(); ++i)
	{
		ci = contactArray.at(i);
		impulse = 0.7*f(i) * ci.N;
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

		/*
		if(f(i) > 0.0)
		{
		double scale = 0.90;
		ci.A->P *= scale;
		ci.B->P *= scale;
		ci.A->L *= scale;
		ci.B->L *= scale;
		}
		*/
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

		if(b(i)>1e3)
		{
			LOG("rAi = ci.P - ci.A->X: " << rAi);
			LOG("wAxrAi = arma::cross(ci.A->W,rAi): " << wAxrAi);
			LOG("At1 = ci.A->inv_mass * ci.A->m_externalForce: " << At1);
			LOG("ci.A->inv_inertia" << ci.A->inv_inertia);
			LOG("ci.A->inv_inertia * (ci.A->m_externalTorque + arma::cross(ci.A->L,ci.A->W)): " << ci.A->inv_inertia * (ci.A->m_externalTorque + arma::cross(ci.A->L,ci.A->W)));
			LOG("At2 = arma::cross(ci.A->inv_inertia * (ci.A->m_externalTorque + arma::cross(ci.A->L,ci.A->W)),rAi): " << At2);
			LOG("At3 = arma::cross(ci.A->W,wAxrAi): " << At3);
			LOG("At4 = ci.A->V + wAxrAi: " << At4);

			LOG("Ndot: " << Ndot);

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
		}

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
		//N3L
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
	//LOG("f: " << f);
	applyImpulse(contactArray,f);


	//LCP STUFF GOES HERE
	//---------------------------
	//LCPSolver (int numEquations, double** M, double* Q, double* Z,double* W, int& status, int maxRetries = 100,double zeroTolerance = 0.0, double ratioError = 0.0);
	//LCPsolver(contactArray.size(),A,b,relAcc,g);
#ifdef checkResting
	bool resting = false;
	for(unsigned int i =0;i<ddot.n_rows;++i)
	{
		if(fabs(ddot(i)) < 1e-3)
		{
			resting = true;
			break;
		}

	}

	if(resting)
	{
#endif

		computeRestingVector(contactArray,b);

		LOG("b: " << b);
		lemke(A,b,relA);
		//LOG("relA: " << relA);
		//LOG("b: " << b);
		//LOG("ArelA+b: " << A*relA+b);

		//LOG("ArelA+ddot: " << A*f+ddot);
		//LOG("c: " << c);
		//g = A*b + relA;
		g = relA.rows(0,contactArray.size()-1);
		if(g(0)!=0)
		{
			arma::vec tmp = A*g+b;
			arma::vec tmp2 = arma::zeros<arma::vec>(g.n_rows,1);
			for(int i = 0; i < g.n_rows;i++)
				tmp2(i) = tmp(i)*g(i);

			LOG("Resting contact!");
			LOG("A: " << A);
			LOG("b: " << b);
			LOG("g: " << g << "<-- Ska vara > 0");
			LOG("ddotdot = Ag+b: " << tmp << " <-- Ska vara > 0");
			LOG("ddotdot boll g: " << tmp2 << "<-- Ska vara = 0");

		}

		//---------------------------

		updateInternalValues(contactArray,g);
#ifdef checkResting
	}
#endif
}


#endif