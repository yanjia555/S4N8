#include <iostream>
#include <armadillo>
#include "Constraints.h"
#include "GaussQuadrature.h"
#include "MeshGenerator.h"
#include "Material.h"
#include "StiffnessMatrix.h"
#include "ForceVector.h"
#include "StaticSolver.h"
using namespace arma;
using namespace std;





int main() {
	double a = 100; // length in x
	double b = 100; // width in y

	int ex = 15; // number of elements in x
	int ey = 15; // number of elements in y

	double tmax = 1; // maximum thickness

 Mesh mesh = MeshGenerator(a, b, ex, ey, tmax);
 Quadrature quadrature = GaussQuadrature("gauss2");


 
 Material material{};//材料
	material.E = 72000; //AL
	material.v = 0.33;
	material.rho = 2.7e-9;
	material.G = material.E / (2 + 2 * material.v);

	string ShapeOption = "Q8";
	mat K = StiffnessMatrix(material, mesh, quadrature, ShapeOption);

	double p0=1e-3;
	mat F=ForceVector(mesh,quadrature,ShapeOption,p0);
	ivec nc = mesh.lato1;


   ConstraintsRes csr = Constraints(nc, K, F);



	 mat w = StaticSolver(csr.K_c, csr.F_c, mesh ,csr.nctot);//solve

	 int k = 0;
	imat Elements = zeros<imat>(mesh.ne, 8);
	for (int i = 0; i < mesh.ne; i+=15) 
	{
		Elements(span(i, i + ey - 1), 0) = regspace<imat>(1 + k * (3 * ey + 2), 2, 1 + k * (3 * ey + 2) + 2 * (ey - 1));
		Elements(span(i, i + ey - 1), 1) = Elements(span(i, i + ey - 1), 0) + 2;

		Elements(span(i, i + ey - 1), 2) = Elements(span(i, i + ey - 1), 1) + 3 * ey + 2;
		Elements(span(i, i + ey - 1), 3) = Elements(span(i, i + ey - 1), 2) - 2;
		Elements(span(i, i + ey - 1), 4) = regspace<imat>(2 + k * (3 * ey + 2), 2, 2 + k * (3 * ey + 2) + 2 * (ey - 1));

		Elements(span(i, i + ey - 1), 5) = regspace<imat>(2 * ey + 3 + k * (3 * ey + 2), 1, 2 * ey + 3 + k * (3 * ey + 2) + ey - 1);
		Elements(span(i, i + ey - 1), 6) = Elements(span(i, i + ey - 1), 2) - 1;
		Elements(span(i, i + ey - 1), 7) = Elements(span(i, i + ey - 1), 5) - 1;
		k++;



   }


	for (int j = 0; j < mesh.ne; ++j) {

		mesh.Eid(j, 0) = j;
		mesh.Eid(j, span(1, 8)) = Elements(j, span(0, 7));

       }

	mat Nodes = mesh.Nid;
	mat gNTu = zeros<mat>(Nodes.n_rows, 3);

	for (int id = 0; id < Nodes.n_rows; ++id) {
		   gNTu(id, span(0, 2)) = w(span(5 * id, 5 * id +2),0).t();
	   }

	double max_disp = max(gNTu(span(0, gNTu.n_rows - 1), 2));
	double min_disp = min(gNTu(span(0, gNTu.n_rows - 1), 2));



	cout << "Max displacement: " << max_disp << endl;  //输出最大最小位移
	cout << "Min displacement: " << min_disp << endl;
	







	return 0;




  }




