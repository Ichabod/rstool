#include <iostream>
#include <fstream>
#include <iomanip>

#include "discr/grid.h"
#include "tools/rtclock.h"

using namespace std;
using namespace Discr;

int main()
{

	Grid grid("deos_result.grid");
	vector<Vertex> vertices;
	vector<Drawable> triangles;
	if (grid.generate3DTriangulation(vertices, triangles, true, 2))
	{

		ofstream vertfile("deos_result_points.txt");
		for (unsigned int i=0; i<vertices.size(); i++)
		{
			double point[3];
			vertices[i].readCoordinates(point);
			vertfile << point[0] << " "<< point[1] << " "<< point[2] << endl;
		}
		vertfile.close();

		ofstream facefile("deos_result_faces.txt");
		for (unsigned int i=0; i<triangles.size(); i++)
		{
			facefile << triangles[i][0] << " "<< triangles[i][1] << " "<< triangles[i][2] << endl;
		}
		facefile.close();

	} else
	{
		cout << "generation failed" << endl;
	}

}
