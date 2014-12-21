
#ifndef ODE_RADAU_H
#define ODE_RADAU_H

#include "solver.h"

namespace ODE
{

	typedef void ( *JACFUNC ) ( int * n, double * t, double * y, double * dfy, int * ldfy, double * rpar, double* ipar );
	typedef void ( *MASSFUNC ) ( int * n, double ** mas, int * lmas, double * rpar, double* ipar );

	/** \class btmb::OdeSol::Radau
		\brief Klasse der Radau-Löser für gewöhnliche Differentialgleichungen
	*/
	class Radau : public Solver
	{
		public:

			/*! Enumeration zur Auswahl der Schrittweitensteuerung */
			typedef enum
			{
				/*! Liefert bei einfachen Problemen im Vergleich zur klassischen Methode sicherere Resultate*/
				Gustaffson,

				/*! Die Klassische Methode ist etwas schneller */
				Classic
			} STEPSIZESTRATEGY;

			/**
			* Integriert nach dem Initialisieren bis zum Zeitpunkt t
			* @param t Endzeitpunkt der Integration
			* @param rpar Benutzerdefinierte Doubleparameter zur freien Verwendung
			* @param ipar Benutzerdefinierte Integerparameter zur freien Verwendung
			* @throw wrongSolverInputException Wird geworfen, wenn die Funktionsparameter ungültig sind (sollte eigentlich nicht vorkommen)
			* @throw tooManyStepsException Wird geworfen, wenn während der Integration zu viele Schritte benötigt wurden
			* @throw tooSmallStepSizeException Wird geworfen, wenn während der Integration die Schrittweite zu klein wird
			* @throw singularyMatrixException Wird geworfen, wenn während der Integration das Newtonverfahren keine eindeutige Lösung mehr finden kann
			 */
			virtual void calc ( double t, double * rpar=0, int* ipar=0 );

			virtual void setMaxSteps ( int steps = 10000000 );
			virtual void setMaxStepsize ( double stepsize = 0.0 );
			virtual void setInitialStepsize ( double init = 0.0 );
			virtual void setStepsizeSelectionParameters ( double fac1=0.2, double fac2=10.0 );

			virtual int getFunctionEvalCount();
			virtual int getComputedSteps();
			virtual int getAcceptedSteps();
			virtual int getRejectedSteps();

			/**
			* Ändert Sicherheitsfaktor/Saftyfactor
			* @param safetyfactor Sicherheitsfaktor ;)
			*/
			virtual void setSafetyFactor ( double safetyfactor = 0.9 );

			/**
			 * Wenn gesetzt, wird die Jacobimatrix standardmäßig erst in Hessenbergform transformiert, was bei großen System von Vorteil sein kann
			 * @param enable Flag
			 */
			void enableHessenbergTransformation ( bool enable );

			/**
			 * Gibt an, wie viele Newtonschritte pro Integrationsschritt maximal erlaubt sind
			 * @param steps Anzahl der erlaubten Newtonschritte
			 * @throw wrongSolverInputException Wird geworfen, wenn eine ungültige Schrittanzahl eingegeben wurde
			 */
			void setMaxNewtonSteps ( int steps = 7 );

			/**
			 * Aktiviert die Verwendung von Startwerten bei dem Newtonverfahren. Diese Funktion sollte genutzt werden,
			 * wenn das Newtonverfahren Konvergenzprobleme hat.
			 * @sa badConvergence()
			 * @param enable Flag
			 */
			void enableStartingValues ( bool enable = false );

			/**
			 * Gibt nach einer Integration an, ob das Newtonverfahren Probleme mit der Konvergenz hatte
			 * @return
			 */
			bool badConvergence();

			/**
			 * Setzt gewisse Parameter, die für das intrene Algebraische System zur Bestimmung der Jacobi- bzw. Hessematrix
			 * @param dim_index1 Parameter 1, standardmäßig die Dimension des Differentialgleichungssystems
			 * @param dim_index2 Parameter 2, standardmäßig 0
			 * @param dim_index3 Parameter 3, standardmäßig 0
			 * @param m1 Parameter 4, standardmäßig 0
			 * @param m2 Parameter 5, standardmäßig 0
			 * @sa http://www.unige.ch/~hairer/software.html
			 */
			void setAlgebraicSystemParameters ( int dim_index1, int dim_index2, int dim_index3, int m1, int m2 );

			/**
			 * Definiert, wie oft die Jacobi-Matrix neu berechnet werden soll. Für große Systeme mit berechnungsintensiver Jacobimatrix
			 * ist 0.1 ein vernünftiger Wert. Negative Werte erzwingen eine Neiberechnung nach jedem Schritt.
			 * @param step Schrittweite der Jacobineuberechnung
			 */
			void setJacobiRecomputation ( double step=0.001 );

			/**
			 * Gibt die Strategie an, mit der die Schrittweiten bestimmt werden
			 * @sa STEPSIZESTRATEGY
			 * @param strategy Schrittweitenstrategie
			 */
			void setStepsizeStrategy ( STEPSIZESTRATEGY strategy = Gustaffson );

			/**
			 * Ist die Bedingung fac1 <= hnew/hold <= fac2 erfüllt, dann wird die Schrittweite nicht geändert.
			 * Für Große Systeme liefern die Parameter 0.99 und 2.0 vernünftigere Resultate
			 * @param fac1 Untergrenze
			 * @param fac2 Obergrenze
			 */
			void setStepsizeFixingParameters ( double fac1 = 1.0, double fac2 = 1.2 );

			virtual ~Radau();

		protected:

			Radau ( ODEFUNC func, int ndgl, JACFUNC jacfunc = 0, int ljac = -1, int ujac = -1, MASSFUNC massfunc = 0, int lmas = -1, int umas = -1 );

			int ijac;		/**< Entscheidungsvariable zur Jacobi-Matrix
							- \f$ ijac = 0 \f$ Berechne Jacobi-Matrix numerisch
							- \f$ ijac = 1 \f$ Berechne Jacobi-Matrix durch vom Nutzer implementierte Funktion "jac" */

			int mljac;		/**< Einstellungsmöglichkeit bei besonderer Struktur der Jacobi-Matrix
							- \f$ mljac = ndgl \f$ Jacobi-Matrix ist voll, Lösung mit Gauss
							- \f$ 0 \leq mljac \leq ndgl \f$ Jacobi-Matrix hat untere Bandbreite mljac */

			int mujac;		/**< Einstellungsmöglichkeit bei besonderer Struktur der Jacobi-Matrix
							- \f$ mujac = ndgl \f$ Jacobi-Matrix ist voll, Lösung mit Gauss
							- \f$ 0 \leq mujac \leq ndgl \f$ Jacobi-Matrix hat obere Bandbreite mljac */

			int imas;		/**< Entscheidungsvariable zur Massen-Matrix
							- \f$ imas = 0 \f$ Massen-Matrix ist Einheitsmatrix
							- \f$ imas = 1 \f$ Massen-Matrix ist durch Nutzer in der Funktion "mass" gegeben */

			int mlmas;		/**< Einstellungsmöglichkeit bei besonderer Struktur der Massen-Matrix
							- \f$ mlmas = ndgl \f$ Massen-Matrix ist voll, Lösung mit Gauss
							- \f$ 0 \leq mlmas \leq ndgl \f$ Massen-Matrix hat untere Bandbreite mljac */

			int mumas;		/**< Einstellungsmöglichkeit bei besonderer Struktur der Massen-Matrix
							- \f$ mumas = ndgl \f$ Massen-Matrix ist voll, Lösung mit Gauss
							- \f$ 0 \leq mumas \leq ndgl \f$ Massen-Matrix hat untere Bandbreite mljac */

			double hinit;		/**< Anfangsschrittweite */

			JACFUNC jac;		/**< Funktion zur Bereitstellung der Jacobi-Matrix durch Nutzer */
			MASSFUNC mass;		/**< Funktion zur Bereitstellung der Massen-Matrix durch Nutzer */

			/**
				Funktionsprototyp zum Aufruf des Differentialgleichungslösers
				\param n Dimension des Problem
				\param FCN Zeiger der Differentialgleichung
				\param x aktuelle Zeit
				\param y aktueller Vektor der Trajektorienwerte
				\param xend Endwert der Trajektorie
				\param hinit Anfangsschrittweite
				\param rtol relative Toleranz
				\param atol absolute Toleranz
				\param itol Festlegungsvariable für Toleranzen
					- \f$ itol = 0 \f$ dann sind atol und rtol skalar
					- \f$ itol \geq 0 \f$ dann sind atol unr rtol Vektoren
				\param jac Funktion zur Bereitstellung der Jacobi-Matrix durch Nutzer
				\param ijac Entscheidungsvariable zur Jacobi-Matrix
							- \f$ ijac = 0 \f$ Berechne Jacobi-Matrix numerisch
							- \f$ ijac = 1 \f$ Berechne Jacobi-Matrix durch vom Nutzer implementierte Funktion "jac"
				\param mljac Einstellungsmöglichkeit bei besonderer Struktur der Jacobi-Matrix
							- \f$ mljac = ndgl \f$ Jacobi-Matrix ist voll, Lösung mit Gauss
							- \f$ 0 \leq mljac \leq ndgl \f$ Jacobi-Matrix hat untere Bandbreite mljac
				\param mujac Einstellungsmöglichkeit bei besonderer Struktur der Jacobi-Matrix
							- \f$ mujac = ndgl \f$ Jacobi-Matrix ist voll, Lösung mit Gauss
							- \f$ 0 \leq mujac \leq ndgl \f$ Jacobi-Matrix hat obere Bandbreite mljac
				\param mas Funktion zur Bereitstellung der Massen-Matrix durch Nutzer
				\param imas Entscheidungsvariable zur Massen-Matrix
							- \f$ imas = 0 \f$ Massen-Matrix ist Einheitsmatrix
							- \f$ imas = 1 \f$ Massen-Matrix ist durch Nutzer in der Funktion "mas" gegeben
				\param mlmas Einstellungsmöglichkeit bei besonderer Struktur der Massen-Matrix
							- \f$ mlmas = ndgl \f$ Massen-Matrix ist voll, Lösung mit Gauss
							- \f$ 0 \leq mlmas \leq ndgl \f$ Massen-Matrix hat untere Bandbreite mljac
				\param mumas Einstellungsmöglichkeit bei besonderer Struktur der Massen-Matrix
							- \f$ mumas = ndgl \f$ Massen-Matrix ist voll, Lösung mit Gauss
							- \f$ 0 \leq mumas \leq ndgl \f$ Massen-Matrix hat untere Bandbreite mljac
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
			void ( *RADAUFUNC ) ( int * n, ODEFUNC FCN, double * x, double * y, double * xend,
			                      double *hinit, double * rtol, double * atol, int * itol,
			                      JACFUNC jac ,int * ijac, int * mljac, int * mujac,
			                      MASSFUNC mass , int * imas, int * mlmas, int * mumas,
			                      SOLOUTFUNC solout, int * iout,
			                      double * work, int* lwork, int * iwork, int * liwork, double * rpar, int * ipar, int * did );

	};

}

#endif
