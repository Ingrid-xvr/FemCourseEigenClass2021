//
//  ShapeQuad.cpp
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "Shape1d.h"
#include "ShapeQuad.h"

/// computes the shape functions in function of the coordinate in parameter space and orders of the shape functions (size of orders is number of sides of the element topology)
void ShapeQuad::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, MatrixDouble &dphi){
    
    for (int i = 0; i < orders.size(); i++)
    {
        if (orders[i] < 0) {
            std::cout << "ShapeQuad::Shape: Invalid dimension for arguments: order\n";
            DebugStop();
        }
    }
    if (orders[0] > 1 || orders[1] > 1 || orders[2] > 1 || orders[3] > 1) {
        std::cout << "ShapeQuad::Shape: Invalid dimension for arguments: order\n";
        DebugStop();
    }

    auto nf = NShapeFunctions(orders);

    if (orders[nf-1] > 2) {
        std::cout << "ShapeQuad::Shape, only implemented until order = 2" << std::endl;
        DebugStop();
    }

    phi.resize(nf);
    dphi.resize(Dimension, nf);

    double qsi = xi[0];
    double eta = xi[1];

    phi[0] = (1.-qsi)*(1.-eta)/4.;
    phi[1] = (1.+qsi)*(1.-eta)/4.;
    phi[2] = (1.+qsi)*(1.+eta)/4.;
    phi[3] = (1.-qsi)*(1.+eta)/4.;

    dphi(0,0) = (eta-1.)/4.;
    dphi(1,0) = (qsi-1.)/4.;
    dphi(0,1) = (1.-eta)/4.;
    dphi(1,1) = (-1.-qsi)/4.;
    dphi(0,2) = (1.+eta)/4.;
    dphi(1,2) = (1.+qsi)/4.;
    dphi(0,3) = (-1.-eta)/4.;
    dphi(1,3) = (1.-qsi)/4.;

//p quadratico
    double phiBolha = 4.*phi[0]*phi[2]; 
    
    VecDouble dphibolha(2);
    dphibolha[0] = 4.*(dphi(0,0)*phi[2]+phi[0]*dphi(0,2));
    dphibolha[1] = 4.*(dphi(1,0)*phi[2]+phi[0]*dphi(1,2));

    int count = 4;

    for(int i=4; i<8; i++){ // arestas
        if(orders[i]==2){
            int aux1 = SideNodeLocIndex(i, 0);
            int aux2 = SideNodeLocIndex(i, 1);

            phi[count] = 4.*phi[aux1]*phi[aux2]+phiBolha;

            dphi(0,count) = 4.*(dphi(0,aux1)*phi[aux2]+phi[aux1]*dphi(0,aux2))+dphibolha[0];
            dphi(1,count) = 4.*(dphi(1,aux1)*phi[aux2]+phi[aux1]*dphi(1,aux2))+dphibolha[1];

            count++; 
        }
    }

    if(orders[8]==2){ // face

        phi[count] = 4.*phiBolha;
        dphi(0,count) = 4.*dphibolha[0];
        dphi(1,count) = 4.*dphibolha[1];
    }
}

/// returns the number of shape functions associated with a side
int ShapeQuad::NShapeFunctions(int side, int order){
    if(order < 1 || order >2) DebugStop();
    if(side<4)
        return 1;//0 a 4
    else if(side<8)
        return (order-1);//6 a 14
    else if(side==8)
        return ((order-1)*(order-1));
    
    std::cout << "ShapeQuad::NShapeFunctions, bad parameter side " << side << std::endl;
    DebugStop();
    
    return 0;
}

/// returns the total number of shape functions
int ShapeQuad::NShapeFunctions(VecInt &orders){
    
    int res=4;
    for(int in=4; in<orders.size(); in++) {
        res += NShapeFunctions(in, orders[in]);
    }
    
    return res;
}
