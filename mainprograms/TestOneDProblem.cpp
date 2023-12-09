//
//  TestOneDProblem.cpp
//  FemSC
//
//  Created by Eduardo Ferri on 08/17/15.
//
//
//TestOneDProblem cpp
/*
        Os testes foram preparados com um proposito educacional,
        recomenda-se que o aluno entenda a funcionalidade de cada
        teste e posteriormente use com seu código caso a caso
*/
//      Obs: O xmax e xmin estao tomados como 4 e 0, respectivamente,
//      caso estes valores sejam alterados, editar o teste TestNodes.
//
//
#include <iostream>
#include <math.h>
#include "GeoMesh.h"
#include "ReadGmsh.h"
#include "CompMesh.h"
#include "GeoElement.h"
#include "GeoElementTemplate.h"
#include "MathStatement.h"
#include "L2Projection.h"
#include "Analysis.h"
#include "IntRule.h"
#include "PostProcessTemplate.h"
#include "Poisson.h"

using std::cout;
using std::endl;
using std::cin;

void exact(const VecDouble &point,VecDouble &val, MatrixDouble &deriv);

int main ()
{
    // Criacao da malha geometrica e preenchimento das propriedades da mesma
    // baseado na malha do Gmsh chamada oneD.msh
    GeoMesh gmesh;
    ReadGmsh read;
    std::string filename("oneD.msh");
#ifdef MACOSX
    filename = "../"+filename;
#endif
    read.Read(gmesh,filename);

    // Criacao de malha computacional e de mathstatement de Poisson com matid 1
    CompMesh cmesh(&gmesh);
    MatrixDouble perm(3,3);
    perm.setZero();
    perm(0,0) = 1.;
    perm(1,1) = 1.;
    perm(2,2) = 1.;
    Poisson *mat1 = new Poisson(1/*matid*/,perm);
    mat1->SetDimension(1);
    
    // Função conhecida da nossa equação diferencial que irá para o vetor de carga
    auto force = [](const VecDouble &x, VecDouble &res)
    {
        res[0] = 1.;
    };
    mat1->SetForceFunction(force);
    
    // Condicoes de contorno usando o conceito de L2 projection
    MatrixDouble proj(1,1),val1(1,1),val2(1,1);
    proj.setZero();
    val1.setZero();
    val2.setZero();
    // 0 bctype significa dirichlet, ou imposição em u
    L2Projection *bc_linha = new L2Projection(0/*bctype*/,2/*matid*/,proj,val1,val2);
    L2Projection *bc_point = new L2Projection(0/*bctype*/,3/*matid*/,proj,val1,val2);
        
    // Setando todos os mathstatements do problema na malha computacional, incluindo,
    // as condicoes de contorno
    std::vector<MathStatement *> mathvec = {nullptr,mat1,bc_linha,bc_point};
    cmesh.SetMathVec(mathvec);
    cmesh.SetDefaultOrder(2);
    cmesh.AutoBuild();
    cmesh.Resequence();

    // Criação do Analysis e rodando a simulação, que envolve Assemble() e Solve()
    Analysis AnalysisLoc(&cmesh);
    AnalysisLoc.RunSimulation();
    
    std::cout << cmesh.Solution() << std::endl;
    
    // Posprocessamento do erro baseado na solução exata (exact)
    PostProcessTemplate<Poisson> postprocess;
    postprocess.SetExact(exact);
    VecDouble errvec;
    errvec = AnalysisLoc.PostProcessError(std::cout, postprocess);

    return 0;
}
void exact(const VecDouble &point,VecDouble &val, MatrixDouble &deriv){

    deriv(0,0) = 4-point[0];
    val[0]=point[0]*(8.-point[0])/2.;
    return;
}


