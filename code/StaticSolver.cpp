#include "StaticSolver.h"





mat StaticSolver(mat K_c,mat F_c,Mesh mesh,imat nctot) 
{
    // s=K_c\F_c;
    mat s = solve(K_c,F_c);

    // index=1:mesh.nn*5;
    // q=zeros(mesh.nn*5,1);
    // index(nctot)=[];%删除被约束的自由度
    //
    // for i=1:numel(index)
    //     q(index(i))=s(i);
    // end
    //
    // w=q;

    mat q = zeros<mat>(mesh.nn*5,1);
    uvec index = regspace<uvec>(0,mesh.nn*5-1);
    uvec nctot_uvec = conv_to<uvec>::from(nctot);



    index.shed_rows(nctot_uvec-1);
    for (int i = 0; i < index.n_elem; ++i) {


        q(index(i)) = s(i);


    }

   
    
    
    
    
    
    return q;
}
