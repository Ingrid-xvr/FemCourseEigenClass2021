#include <iostream>
#include <math.h>
#include "VTKGeoMesh.h"
#include "Geom1d.h"
#include "CompElementTemplate.h"
#include "Shape1d.h"
#include "CompMesh.h"
#include "GeoElementTemplate.h"
#include "GeoMesh.h"
#include "Poisson.h"

using namespace std;


void CallStiffFunc();

int main() {
    CallStiffFunc();
    return 0;
}

void CallStiffFunc(){
    
    
    // Criar nós com coordenadas e adicionando na malha geometrica
    const int dim = 1;
    const int numnodes = 2;
    GeoMesh gmesh;
    gmesh.SetDimension(dim);
    gmesh.SetNumNodes(numnodes);
    GeoNode gnod0, gnod1;
    VecDouble co0(1), co1(1);
    co0 << 0.;
    co1 << 4.;
    gnod0.SetCo(co0);
    gnod1.SetCo(co1);
    gmesh.Node(0) = gnod0;
    gmesh.Node(1) = gnod1;
    
    // Criando elemento geometrico e se auto adicionando na malha geometrica no construtor
    int materialid = 0;
    VecInt nodeindices(2);// vetor com indices dos nos na malha
    nodeindices << 0,1;
    int index = 0;
    GeoElementTemplate<Geom1d> geo(nodeindices,materialid,&gmesh,index);
        
    // Criando malha computacional e setando a ordem polinomial padrap
    CompMesh cmesh(&gmesh);
    cmesh.SetDefaultOrder(1);
    
    // Criando Mathstatement de equacao de Poisson e adicionando na malha computacional
    MatrixDouble perm(3,3);
    perm.setIdentity();
    Poisson poi(materialid, perm);
    cmesh.SetMathStatement(materialid, &poi);
    
    // Criando elemento computacional baseado em elemento geometrico, se adicionando na malha computacional
    // e calculando a matriz de rigidez e vetor de carga.
    CompElementTemplate<Shape1d> cel(index,&cmesh,&geo);
    MatrixDouble ek(2,2),ef(2,1);
    cel.CalcStiff(ek, ef);
    cout << ek << endl;
}
