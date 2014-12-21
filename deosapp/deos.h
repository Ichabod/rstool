#ifndef DEOSAPPDEOS_H
#define DEOSAPPDEOS_H

#include <IpIpoptApplication.hpp>
#include "../algorithm/initialguess.h"
#include "../algorithm/application.h"

namespace DEOSApp
{

	class DEOS : public Algorithm::Application
	{
	public:
		DEOS(int N, double hstep, double c_rx = 7071000.0, double c_G = 398e12, double c_Ms = 200.0,
			double c_js1 = 2000.0, double c_js2 = 5000.0, double c_js3 = 2000.0,
			double c_dsx = 0.0, double c_dsy = 1.0, double c_dsz = 0.0,	double c_vmax = 0.15, double c_mmax = 1.0, double c_safety = 1.9);

		virtual ~DEOS();

		bool generateTarget(double * x0, double c_Mt = 200.0, double c_jt1 = 1000.0, double c_jt2 = 2000.0, double c_jt3 = 1000.0, double c_dtx = 0.0, double c_dty = -1.0, double c_dtz = 0.0);
		virtual int generateInitialGuess(Algorithm::InitialGuess * & initialguess, double * poshint);

		int getInitialGuessSize() const
		{
			return N*6 + (N+1)*13 + N*13*stages;
		}

		virtual int horizonSize() const
		{
			return N;
		}

		virtual int stepSize() const
		{
			return hstep;
		}

		virtual double startTime() const
		{
			return 0.0;
		}

		virtual double * getParameters()
		{
			return par;
		}

		double * getTargetPos()
		{
			return target;
		}

		double * getTargetDockpoint()
		{
			return Rtd;
		}


	protected:
		static void ODEwrapper( int * n, double * t, double * y, double * dy, double * rpar, int * ipar );

		bool rkstepWrapper(double hstep, double * x, double * x_to, double * k, double * par);
		bool rkstepWrapper(double hstep, double * x, double * par);

		bool IpoptWrapper(Ipopt::TNLP * nlp);

		static void printVec(int step, double * x, int n);

		int N;
		double hstep;

		double c_rx;
		double c_G;
		double c_Ms;
		double c_js1;
		double c_js2;
		double c_js3;
		double c_dsx;
		double c_dsy;
		double c_dsz;
		double c_safety;
		double c_vmax;
		double c_mmax;

		Ipopt::SmartPtr<Ipopt::IpoptApplication> ipoptapp; 


		/* Targetberechnung */
		double * par;
		double * target;
		double * vtargetN;
		double * Rtd;
		double * dRtd;

		/* Dockingpoint berechnung */
		int * nqx;
		int * nrx;
		double * dqxNtemp;
		double * drxNtemp;

		/* Trajectory berechnung */
		double * mpcsolution;
		int mpcN;
		double * mpcstart;
		double * rkA;
		double * rkb;
		double * rkc;
		double * k;
		double * rkwork;
		int * rkiwork;
		int stages;


	};

}


#endif