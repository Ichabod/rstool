#ifndef ODERADAU5913_H
#define ODERADAU5913_H

#include "radau.h"

namespace ODE
{

	class Radau5913 : public Radau
	{
		public:
			/**
			 * Initialisiert den Differentialgleichungslöser
			 * @param func Zeiger auf die "rechte Seite" der Differentialgleichung
			 * @param ndgl Dimension des Problems
			 * @param jacfunc Zeiger auf die Funktion der Jacobimatrix, wenn leer wird die Matrix algebraisch berechnet.
			 * @param ljac Bandbreite der unteren Dreiecksmatrix der Jacobimatrix (ngdl oder -1 für volle Matrix)
			 * @param ujac Bandbreite der oberen Dreiecksmatrix der Jacobimatrix (ngdl oder -1 für volle Matrix)
			 * @param massfunc Zeiger auf die Funktion der Mass-Matrix, wenn leer wird die Einheitsmatrix angenommen
			 * @param lmas Bandbreite der unteren Dreiecksmatrix der Mass-Matrix (ngdl oder -1 für volle Matrix)
			 * @param umas Bandbreite der oberen Dreiecksmatrix der Mass-Matrix (ngdl oder -1 für volle Matrix)
			 * @throw memoryException Wird geworfen, wenn kein Speicher für das Problem reserviert werden konnte.
			 */
			Radau5913 ( ODEFUNC func, int ndgl, JACFUNC jacfunc = 0, int ljac = -1, int ujac = -1, MASSFUNC massfunc = 0, int lmas = -1, int umas = -1 );

			/**
			 * Setzt Parameter für die Ordnungsbestimmung des Verfahrens
			 * @param inc_factor Die Ordnung wird erhöht, wenn die Kontraktion kleiner als dieser Parameter ist
			 * @param dec_factor Die Ordnung wird Verringert, wenn die Kontraktion größer als dieser Parameter ist
			 * @param hfac1 Die Ordnung wird nur dann verringert, wenn gilt: hfac1 <= hnew/hold <= hfac2
			 * @param hfac2 siehe hfac2
			 */
			void setOrderSelectionParameters ( double inc_factor = 0.002, double dec_factor = 0.8, double hfac1 = 1.2, double hfac2 = 0.8 );

			/**
			 * Setzt die Grenzen der Stages (gültige Werte sind jeweils  nur 1, 3, 5 und 7)
			 * @param minstage Kleinste Anzahl der Stages
			 * @param maxstage Größte Anzahl der Stages
			 * @param firststage Stageanzahl für den ersten Schritt
			 * @throw wrongSolverInputException Wird geworfen, wenn die Stageangaben keine gültigen Werte sind
			 */
			void setStages ( int minstage = 3, int maxstage= 7, int firststage = 3 );

			/**
				Funktion zur Festlegung von Standardwerten für notwendige Parameter
			*/
			void setDefaults();

		private:

			static bool isValidStage ( int stage );

			int maxstages;

	};

}

#endif
