#include "Constraints.h"


ConstraintsRes Constraints(ivec nc,mat K,mat F) {
    ConstraintsRes cr;
    mat t=zeros(1,numel(nc)*5);
      imat nctot = conv_to<imat>::from(t);


    for (int i = 0; i < numel(nc); ++i) {


        nctot(0,span(i*5,4+i*5)) = {nc(i)*5-4,nc(i)*5-3,nc(i)*5-2,
         nc(i)*5-1,nc(i)*5};


  }


    mat K_c = K;
    mat F_c = F;
    uvec nctot_uvec = conv_to<uvec>::from(nctot);
    K_c.shed_rows(nctot_uvec-1);
    K_c.shed_cols(nctot_uvec-1);

     F_c.shed_rows(nctot_uvec-1);
     cr.nctot = nctot;
     cr.K_c = K_c;
     cr.F_c = F_c;













    return cr;



}
