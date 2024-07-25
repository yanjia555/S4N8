#include "ForceVector.h"

mat ForceVector(Mesh mesh, Quadrature quadrature, string shapeOption, double p0) {


	mat F = zeros(mesh.nn * 5, 1);
	F(720 * 5 - 2 - 1) = 1000;




	return F;

}