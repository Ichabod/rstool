#ifndef ODERADAU5_H
#define ODERADAU5_H

#include "radau.h"

namespace ODE
{

	class Radau5 : public Radau
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
			Radau5 ( ODEFUNC func, int ndgl, JACFUNC jacfunc = 0, int ljac = -1, int ujac = -1, MASSFUNC massfunc = 0, int lmas = -1, int umas = -1 );

			/**
				Funktion zur Festlegung von Standardwerten für notwendige Parameter
			*/
			void setDefaults();

		private:
	};

}

#endif
