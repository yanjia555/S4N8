#include "MeshGenerator.h"

Mesh MeshGenerator(double a, double b, int ex, int ey, double tmax) {
	Mesh mesh;
	
	int nx = 2 * ex + 1; //x方向node
	int ny = 2 * ey + 1;//y方向node

    mesh.ne = ex * ey;
	mesh.nn = nx * (ey + 1) + (ex + 1) * ey;

	double dx = a / ex; //x方向单元长度
	double dy = b / ey;
	

	//得到xyz的坐标
	vec x = linspace(-a / 2, b / 2, 2 * ex + 1); // node coordinates from left to right
	ivec rt = reshape(repmat(ivec({2 * ey + 1, ey + 1}), 1, ex), ex * 2, 1);//每列节点数量 7 4 7 4 7 4……

	rt.insert_rows(rt.n_rows, 1); //
	rt(rt.n_rows - 1) = 2 * ey + 1;
	
	vec x_expanded;
	for (int i = 0; i < x.n_rows; ++i) {
	x_expanded.insert_rows(x_expanded.n_rows, rt(i));
		x_expanded.subvec(x_expanded.n_rows - rt(i), x_expanded.n_rows - 1).fill(x(i));
	}




	x = x_expanded;
	 vec y;
	for (int i = 0; i < ex; ++i) {
		y = join_cols(y, linspace(-b / 2, b / 2, 2 * ey + 1));	
		y = join_cols(y, linspace(-b / 2, b / 2, ey + 1));


	}
	y = join_cols(y, linspace(-b / 2, b / 2, 2 * ey + 1));
	
	vec z = 0.01 * square(y);
	vec t = tmax - abs(y) * 0.0;//厚度先采用均匀厚度
	
	mat n_vector(mesh.nn, 3);
	for (int i = 0; i < mesh.nn; ++i) {
		vec temp = vec({0, -0.02 * y(i), 1});
		n_vector.row(i) = trans(temp / norm(temp));
}
	
	
	mesh.coordinates.set_size(mesh.nn, 3);
	mesh.coordinates.col(0) = x; // x coordinates
	mesh.coordinates.col(1) = y; // y coordinates
	mesh.coordinates.col(2) = z; // z coordinates

	mesh.H3i = n_vector; // normal vector of each node point //论文公式14
	mesh.t = t; // thickness of each node point
	vec ii = vec({1, 0, 0}); // x-axis
	mat H2i(mesh.nn, 3, fill::zeros);
	
	for (int i = 0; i < mesh.nn; ++i) {

		vec temp = trans(cross(n_vector.row(i), ii));
		H2i.row(i) = temp.t() / norm(temp);
}
	
	mat H1i(mesh.nn, 3, fill::zeros);
	for (int i = 0; i < mesh.nn; ++i) {

		vec temp = trans(cross(n_vector.row(i), H2i.row(i)));
		 H1i.row(i) = temp.t() / norm(temp);

	}
	
	mesh.H2i = H2i; //
	mesh.H1i = H1i;
	
	//壳的上下表面
	mat x_up = x + n_vector.col(0) % (0.5 * t);
	mat y_up = y + n_vector.col(1) % (0.5 * t);
	mat z_up = z + n_vector.col(2) % (0.5 * t);
	
	mat x_bot = x - n_vector.col(0) % (0.5 * t);
	mat y_bot = y - n_vector.col(1) % (0.5 * t);
	mat z_bot = z - n_vector.col(2) % (0.5 * t);
	

	////上下表面
	 mesh.up_coordinates.set_size(mesh.nn, 3);
	 mesh.up_coordinates.col(0) = x_up;
	 mesh.up_coordinates.col(1) = y_up;
	 mesh.up_coordinates.col(2) = z_up;
	
	mesh.bot_coordinates.set_size(mesh.nn, 3);
	mesh.bot_coordinates.col(0) = x_bot;
	mesh.bot_coordinates.col(1) = y_bot;
	
	mesh.bot_coordinates.col(2) = z_bot;



	mesh.lato1 = linspace<ivec>(1, nx, nx); //边
	ivec num_temp;
	for (int i = 0; i < ex; ++i) {

		num_temp.insert_rows(num_temp.n_rows, 2);
		num_temp(num_temp.n_rows - 2) = ny + 1 + i * (ny + ey + 1);
		num_temp(num_temp.n_rows - 1) = ny + 4 + i * (ny + ey + 1);
}  


	mesh.lato2 = num_temp;
	mesh.lato3 = linspace<ivec>(mesh.nn - ny + 1, mesh.nn, ny);
	
	mesh.Nid.set_size(mesh.nn, 4);
	mesh.Nid.zeros();
	mesh.Eid.set_size(mesh.ne, 9);
	mesh.Eid.zeros();
	
	for (int j = 0; j < mesh.nn; ++j) {
		mesh.Nid(j, 0) = j + 1;
		mesh.Nid(j, 1) = x(j);
		mesh.Nid(j, 2) = y(j);
		mesh.Nid(j, 3) = z(j);
	      }
	
	int k = 0;
	imat Elements(mesh.ne, 8, fill::zeros);  //ne行，8列
	
	for (int i = 0; i < mesh.ne; i += ey) {
		for (int j = 0; j < ey; ++j) {
			int row = i + j;
		Elements(row, 0) = 1 + k * (3 * ey + 2) + 2 * j;
		Elements(row, 1) = Elements(row, 0) + 2;
		Elements(row, 2) = Elements(row, 1) + 3 * ey + 2;
		Elements(row, 3) =  Elements(row, 2) - 2;

		Elements(row, 4) = 2 + k * (3 * ey + 2) + 2 * j;
		Elements(row, 5) = 2 * ey + 3 + k * (3 * ey + 2) + j;
		Elements(row, 6) = Elements(row, 2) - 1;
	    Elements(row, 7) =  Elements(row, 5) - 1;

		}

		k++;

	}

	for (int j = 0; j < mesh.ne; ++j)
	{


		mesh.Eid(j, 0) = j + 1;
		mesh.Eid.row(j).cols(1, 8) = Elements.row(j);


	 }
	












	return mesh;
	
}