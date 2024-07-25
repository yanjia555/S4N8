#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H
#include <armadillo>
using namespace arma;



struct ConstraintsRes {

    imat nctot;
    mat K_c;
    mat F_c;

};


ConstraintsRes Constraints(ivec nc,mat K,mat F);





#endif //CONSTRAINTS_H
