#ifndef MESHGENERATOR_H
#define MESHGENERATOR_H
#include <string>
#include <armadillo>





using namespace arma;
using namespace std;


struct Mesh {


	int ne; // number of elements
	int nn; // number of nodes (including mid node)
	mat coordinates; // node coordinates (x, y, z)



	mat H1i; // normal vector of each node point
	mat H2i; // tangent vector of each node point
	mat H3i; // binormal vector of each node point
	mat t; // thickness of each node point

	mat up_coordinates; // top surface coordinates (x, y, z)
	mat bot_coordinates; // bottom surface coordinates (x, y, z)
	ivec lato1; // left side node indices
	ivec lato2; // down side node indices

	ivec lato3; // right side node indices
	mat Nid; // node ID matrix (node_id, x, y, z)
	imat Eid; // element ID matrix (ele_id, node1, node2, node3, node4)










};

Mesh MeshGenerator(double a, double b, int ex, int ey, double tmax);

#endif //MESHGENERATOR_H
