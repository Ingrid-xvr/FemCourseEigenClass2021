/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "GeomQuad.h"

GeomQuad::GeomQuad() {
}

GeomQuad::~GeomQuad() {
}

GeomQuad::GeomQuad(const GeomQuad &copy) {
    fNodeIndices = copy.fNodeIndices;
}

GeomQuad& GeomQuad::operator=(const GeomQuad& copy) {
    fNodeIndices = copy.fNodeIndices;
    return *this;
}

void GeomQuad::Shape(const VecDouble &xi, VecDouble &phi, MatrixDouble &dphi) {
    if(xi.size() != Dimension || phi.size() != nCorners || dphi.rows() != Dimension || dphi.cols() != nCorners) DebugStop();
    
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
}

void GeomQuad::X(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x) {
    if(xi.size() != Dimension) DebugStop();
    if(x.size() < NodeCo.rows()) DebugStop();
    if(NodeCo.cols() != nCorners) DebugStop();
    
    VecDouble phi(nCorners);
    MatrixDouble dphi(Dimension, nCorners); 

    Shape(xi, phi, dphi);

    x.setZero();

    for(int i=0; i<nCorners; i++){
        x[0] += phi[i]*NodeCo(0, i); 
        x[1] += phi[i]*NodeCo(1, i);
    }

//    std::cout << "Xi: "<< std::endl << xi << std::endl;
//     std::cout << "Phi: "<< std::endl << phi << std::endl;
//     std::cout << "x: "<< std::endl << x << std::endl;
//     std::cout << "NodeCo: "<< std::endl << NodeCo << std::endl;
}

void GeomQuad::GradX(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x, MatrixDouble &gradx) {
    if(xi.size() != Dimension) DebugStop();
    if(x.size() != NodeCo.rows()) DebugStop();
    if(NodeCo.cols() != nCorners) DebugStop();
    
    gradx.resize(2,2);
    gradx.setZero();
    
    VecDouble phi(nCorners);
    MatrixDouble dphi(Dimension,nCorners);

    Shape(xi,phi,dphi);

    for(int i = 0; i<Dimension;i++){
        for(int j = 0; j<nCorners;j++){
            gradx(i,0) += NodeCo(i,j)*dphi(0,j); // dxi/dqsi
            gradx(i,1) += NodeCo(i,j)*dphi(1,j); // dxi/deta
        } // xi := x (i=0), y (i=1)
    }
}

void GeomQuad::SetNodes(const VecInt &nodes) {
    if(nodes.size() != nCorners) DebugStop();
    fNodeIndices = nodes;
}

void GeomQuad::GetNodes(VecInt &nodes) const{
    nodes = fNodeIndices;
}

int GeomQuad::NodeIndex(int node) const {
    return fNodeIndices[node];
}

int GeomQuad::NumNodes() {
    return nCorners;
}

GeoElementSide GeomQuad::Neighbour(int side) const {
    return fNeighbours[side];
}

void GeomQuad::SetNeighbour(int side, const GeoElementSide &neighbour) {
    fNeighbours[side] = neighbour;
}
