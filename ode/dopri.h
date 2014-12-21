
#ifndef ODEDOPRI_H
#define ODEDOPRI_H

#include "solver.h"

namespace ODE
{

	class DoPri : public Solver
	{
		public:

			/**
			* Integriert nach dem Initialisieren bis zum Zeitpunkt t
			* @param t Endzeitpunkt der Integration
			* @param rpar Benutzerdefinierte Doubleparameter zur freien Verwendung
			* @param ipar Benutzerdefinierte Integerparameter zur freien Verwendung
			* @throw wrongSolverInputException Wird geworfen, wenn die Funktionsparameter ungültig sind (sollte eigentlich nicht vorkommen)
			* @throw tooManyStepsException Wird geworfen, wenn während der Integration zu viele Schritte benötigt wurden
			* @throw tooSmallStepSizeException Wird geworfen, wenn während der Integration die Schrittweite zu klein wird
			* @throw stiffException Wird geworfen, wenn eine steife Funktion integriert werden soll
			*/
			virtual void calc ( double t, double * rpar=0, int* ipar=0 );

			virtual void setMaxSteps ( int steps = 10000000 );
			virtual void setMaxStepsize ( double stepsize = 0.0 );
			virtual void setInitialStepsize ( double init = 0.0 );
			virtual void setStepsizeSelectionParameters ( double fac1=0.2, double fac2=10.0 );

			/**
			* Setzt Sicherheitsfaktor/Saftyfactor
			* @param safetyfactor Sicherheitsfaktor ;)
			*/
			virtual void setSafetyFactor ( double safetyfactor = 0.9 );

			/**
			 * Gibt an, über wie viele Schritte hinweg das Problem auf Steifheit getestet werden soll
			 * @param teststeps Schrittanzahl des Steifheitstests
			 * @throw wrongSolverInputException Es wurde eine ungültige Schrittanzahl eingegeben (negativ?)
			 */
			void setStiffnessTest ( int teststeps = 10 );

			/**
			 * Setzt den Stabilisierungsparameter beta
			 * @param beta Neuer Wert für beta (0 <= beta <= 0.1)
			 * @throw wrongSolverInputException Es wurde eine ungültiger Wert für Beta angegeben
			 */
			void setBeta ( double beta );

			virtual int getFunctionEvalCount();
			virtual int getComputedSteps();
			virtual int getAcceptedSteps();
			virtual int getRejectedSteps();

			/**
				Funktion zur Festlegung von Standardwerten für notwendige Parameter
			*/
			void setDefaults();

			/**
				Destruktor
			*/
			virtual ~DoPri();

		protected:

			/**
				Konstruktor
			*/
			DoPri ( ODEFUNC func, int ndgl );

			/**
				Funktionsprototyp zum Aufruf des Differentialgleichungslösers
				\param n Dimension des Problem
				\param FCN Zeiger der Differentialgleichung
				\param x aktuelle Zeit
				\param y aktueller Vektor der Trajektorienwerte
				\param xend Endwert der Trajektorie
				\param rtol relative Toleranz
				\param atol absolute Toleranz
				\param itol Festlegungsvariable für Toleranzen
					- \f$ itol = 0 \f$ dann sind atol und rtol skalar
					- \f$ itol \geq 0 \f$ dann sind atol unr rtol Vektoren
				\param solout Ausgabefunktion
				\param iout Ausgabeparameter
					- \f$ iout \geq 1 \f$ Ausgabe des Vektor der Trajektorienwerte für jeden akzeptierten Schritt
					- \f$ iout = 0 \f$ keine Ausgabe, Dummy-Funktion ist notwendig
				\param work reservierter Arbeitsspeicher für double-Werte
				\param lwork Größe des reservierten Arbeitsspeichers für double-Werte
				\param iwork reservierter Arbeitsspeicher für integer-Werte
				\param liwork Größe des reservierten Arbeitsspeichers für integer-Werte
				\param rpar Array der Double-Parameter des Problems
				\param ipar Array der Integer-Parameter des Problems
				\param did Rückgabewert der Routine:
					- \f$ did = 1 \f$ Erfolgreiche Berechnung
					- \f$ did = 2 \f$ Erfolgreiche Berechnung, unterbrochen durch Ausgabe mit SOLOUT
					- \f$ did = -1 \f$ Eingabe ist nicht konsistent
					- \f$ did = -2 \f$ Maximalzahl der durchzuführenden Schritte ist zu gering
					- \f$ did = -3 \f$ Schrittweite wird zu klein
					- \f$ did = -4 \f$ Problemstellung ist steif
			*/
			void ( *DOPRIFUNC ) ( int * n, ODE::ODEFUNC FCN, double * x, double * y, double * xend, double * rtol, double * atol, int * itol, ODE::SOLOUTFUNC solout, int * iout, double * work, int* lwork, int * iwork, int * liwork, double * rpar, int * ipar, int * did );

	};

}

#endif
