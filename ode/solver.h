#ifndef ODESOLVER_H
#define ODESOLVER_H

namespace ODE
{

	typedef void ( *ODEFUNC ) ( int * n, double * t, double * y, double * dy, double * rpar, int * ipar );
	typedef void ( *SOLOUTFUNC ) ( int * nr, double * t_old, double * t, double * y, int * n, double * cont, int * icont, int * nrd, double * rpar, int* ipar, int* ireturn );

	class Solver
	{
		public:
			/**
				Destruktor
			*/
			virtual ~Solver();

			/**
			 * Initialisiert den Integrator.
			 * @param t0 Startzeitpunkt
			 * @param y0 Startwert per Referenz, Beinhaltet nach einer Integration die Lösung zum Zeitpunkt t
			 * @param dy0 Startwert per Referenz, Beinhaltet hier nicht, ist aber nötig, damit init falls nötig überschrieben werden kann.
			 */
			virtual void init ( double t0, double * y0);
	
			/**
			* Setzt Fehlertoleranzen für die Schrittweitensteuerung.
			* @param rtol Relative Fehlertoleranz
			* @param atol Absolute Fehlertoleranz
			*/
			void setTol(double rtol, double atol);

			/**
			 * Integriert nach dem Initialisieren bis zum Zeitpunkt t
			 * @param t Endzeitpunkt der Integration
			 * @param rpar Benutzerdefinierte Doubleparameter zur freien Verwendung
			 * @param ipar Benutzerdefinierte Integerparameter zur freien Verwendung
			 * @throw wrongSolverInputException Wird geworfen, wenn die Funktionsparameter ungültig sind (sollte eigentlich nicht vorkommen)
			 * @throw tooManyStepsException Wird geworfen, wenn während der Integration zu viele Schritte benötigt wurden
			 * @throw tooSmallStepSizeException Wird geworfen, wenn während der Integration die Schrittweite zu klein wird
			 */
			virtual void calc ( double t, double * rpar=0, int* ipar=0 ) = 0;

			/**
			 * Setzt die Höchstgrenze an Schritten, die pro Aufruf von calc() durchgeführt werden dürfen.
			 * @param steps Höchstgrenze an Schritten
			 * @throw wrongSolverInputException Es wurde eine ungültige Höchstgrenze eingegeben (negativ?)
			 */
			virtual void setMaxSteps ( int steps = 10000000 ) = 0;

			/**
			 * Ändert Parameter, die für die Schrittweitenerkennung erforderlich sind
			 * @param fac1 Gibt eine Begrenzung für die neue Schrittweite an: fac1 <= h_new/h_old <= fac2
			 * @param fac2 siehe Parameterbeschreibung für fac1
			 */
			virtual void setStepsizeSelectionParameters ( double fac1=0.2, double fac2=10.0 ) = 0;

			/**
			 * Setzt eine Obergrenze für die Schrittweite.
			 * @param stepsize Maximale Schrittweite, wird der Parameter gleich 0 gesetzt, dann ist das Integrationsintervall die maximale Schrittweite.
			 * @throw wrongSolverInputException Es wurde eine ungültige Schrittweite
			 */
			virtual void setMaxStepsize ( double stepsize = 0.0 ) = 0;

			/**
			 * Setzt die Startschrittweite
			 * @param hinit Startschrittweite, wird hinit = 0.0 gesetzt, dann wird eine Schätzung für die Startschrittweite jeweils numerisch bestimmt.
			 * @throw wrongSolverInputException Es wurde eine ungültige Schrittweite
			 */
			virtual void setInitialStepsize ( double hinit = 0.0 ) = 0;

			/**
			 * Gibt nach einer Integration die Anzahl der Funktionsauswertungen zurück
			 * @return Anzahl der Funktionsauswertungen
			 */
			virtual int getFunctionEvalCount() = 0;

			/**
			 * Gibt nach einer Integration die Anzahl der berechneten Schritte zurück
			 * @return Anzahl der Schritte
			 */
			virtual int getComputedSteps() = 0;

			/**
			 * Gibt nach einer Integration die Anzahl der angenommenen Schrittweiten zurück
			 * @return Anzahl der angenommenen Schrittweiten
			 */
			virtual int getAcceptedSteps() = 0;

			/**
			 * Gibt nach einer Integration die Anzahl der abgelehnten Schrittweiten zurück
			 * @return Anzahl der abgelehnten Schrittweiten
			 */
			virtual int getRejectedSteps() = 0;

		protected:
			/**
				Konstruktor
			*/
			Solver ( ODEFUNC func, int ndgl );
			Solver() {}

			ODEFUNC f;		/**< Funktion der rechten Seite */
			SOLOUTFUNC solout;	/**< Funktion für kontinuierliche Ausgabe */

			bool initialized;	/**< Bool-Variable zur Überprüfung der Initialisierung eines Problems */
			int ndgl;		/**< Dimension des Problems */
			double t;		/**< aktueller Zeitpunkt */
			double * y;		/**< Vektor der Trajektorienwerte zum aktuellen Zeitpunkt */

			int iout;		/**< Ausgabeparameter
							- \f$ iout \geq 1 \f$ Ausgabe des Vektor der Trajektorienwerte für jeden akzeptierten Schritt
							- \f$ iout = 0 \f$ keine Ausgabe, Dummy-Funktion ist notwendig */
			int did;		/**< Rückgabewert der Routine:
							- \f$ did = 1 \f$ Erfolgreiche Berechnung
							- \f$ did = 2 \f$ Erfolgreiche Berechnung, unterbrochen durch Ausgabe mit SOLOUT
							- \f$ did = -1 \f$ Eingabe ist nicht konsistent
							- \f$ did = -2 \f$ Maximalzahl der durchzuführenden Schritte ist zu gering
							- \f$ did = -3 \f$ Schrittweite wird zu klein
							- \f$ did = -4 \f$ Problemstellung ist steif */
			int itol;		/**< Festlegungsvariable für Toleranzen
						- \f$ itol = 0 \f$ dann sind atol und rtol skalar
						- \f$ itol \geq 0 \f$ dann sind atol unr rtol Vektoren */
			double rtol;		/**< relative Toleranz */
			double atol;		/**< absolute Toleranz */

			double *work;		/**< reservierter Arbeitsspeicher für double-Werte */
			int lwork;		/**< Größe des reservierten Arbeitsspeichers für double-Werte */
			int * iwork;		/**< reservierter Arbeitsspeicher für integer-Werte */
			int liwork;		/**< Größe des reservierten Arbeitsspeichers für integer-Werte */

	};

}

#endif
