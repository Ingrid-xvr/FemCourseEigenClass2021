/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "GeomTriangle.h"

GeomTriangle::GeomTriangle() {
}

GeomTriangle::~GeomTriangle() {
}

GeomTriangle::GeomTriangle(const GeomTriangle &copy) {
    fNodeIndices = copy.fNodeIndices;

}

GeomTriangle& GeomTriangle::operator=(const GeomTriangle& copy) {
    fNodeIndices = copy.fNodeIndices;

    return *this;
}

void GeomTriangle::Shape(const VecDouble& xi, VecDouble& phi, MatrixDouble& dphi) {
    if(xi.size() != Dimension || phi.size() != nCorners || dphi.rows() != Dimension || dphi.cols() != nCorners) DebugStop();
    
    double qsi = xi[0];
    double eta = xi[1];

    phi[0] = 1.-qsi-eta;
    phi[1] = qsi;
    phi[2] = eta;
    
    dphi(0,0) = -1.;
    dphi(1,0) = -1.;
    dphi(0,1) = 1.;
    dphi(1,1) = 0.;
    dphi(0,2) = 0.;
    dphi(1,2) = 1.;

}

void GeomTriangle::X(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x) {
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
}

void GeomTriangle::GradX(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x, MatrixDouble &gradx) {
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

void GeomTriangle::SetNodes(const VecInt &nodes) {
    if(nodes.size() != nCorners) DebugStop();
    fNodeIndices = nodes;
}

void GeomTriangle::GetNodes(VecInt &nodes) const  {
    nodes = fNodeIndices;
}

int GeomTriangle::NodeIndex(int node) const {
    return fNodeIndices[node];
}

int GeomTriangle::NumNodes() {
    return nCorners;
}

GeoElementSide GeomTriangle::Neighbour(int side)  const {
    return fNeighbours[side];
}

void GeomTriangle::SetNeighbour(int side, const GeoElementSide &neighbour) {
    fNeighbours[side] = neighbour;
}
