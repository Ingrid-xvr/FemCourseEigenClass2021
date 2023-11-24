/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

///\cond
#include <iostream> 
///\endcond
#include "IntRuleTriangle.h"

IntRuleTriangle::IntRuleTriangle(){
    SetOrder(0);
}

IntRuleTriangle::IntRuleTriangle(int order) {
    SetOrder(order);
}

void IntRuleTriangle::SetOrder(int order) {
    
   fOrder = order; 

   switch (order)
   {
   case 0:
   case 1:
        fPoints.resize(1,2);
        fWeights.resize(1);

        fPoints(0,0) = 1./3.;
        fPoints(0,1) = 1./3.;

        fWeights[0] = 1./2.; //o peso é 1, mas divide por 2 por ser triangulo

    break;
   case 2:
        fPoints.resize(3,2);
        fWeights.resize(3);

        fPoints(0,0) = 1./6.;
        fPoints(0,1) = 2./3.;
        fPoints(1,0) = 1./6.;
        fPoints(1,1) = 1./6.;
        fPoints(2,0) = 2./3.;
        fPoints(2,1) = 1./6.;
       
        fWeights[0] = 1./6.;
        fWeights[1] = 1./6.;
        fWeights[2] = 1./6.;
        
        break;
   case 3:
        fPoints.resize(4,2);
        fWeights.resize(4);

        fPoints(0,0) = 1./3.;
        fPoints(0,1) = 1./3.;
        fPoints(1,0) = 1./5.;
        fPoints(1,1) = 3./5.;
        fPoints(2,0) = 1./5.;
        fPoints(2,1) = 1./5.;
        fPoints(3,0) = 3./5.;
        fPoints(3,1) = 1./5.;
       
        fWeights[0] = -0.5625/2.;
        fWeights[1] = 0.520833333333333/2.;
        fWeights[2] = 0.520833333333333/2.;
        fWeights[3] = 0.520833333333333/2.;
        break;

    case 4:
        fPoints.resize(6,2);
        fWeights.resize(6);

        fPoints(0,0) = 0.445948490915965;
        fPoints(0,1) = 0.108103018168070;
        fPoints(1,0) = 0.445948490915965;
        fPoints(1,1) = 0.445948490915965;
        fPoints(2,0) = 0.108103018168070;
        fPoints(2,1) = 0.445948490915965;
        fPoints(3,0) = 0.091576213509771;
        fPoints(3,1) = 0.816847572980459;
        fPoints(4,0) = 0.091576213509771;
        fPoints(4,1) = 0.091576213509771;
        fPoints(5,0) = 0.816847572980459;
        fPoints(5,1) = 0.091576213509771;

        fWeights[0] = 0.223381589678011/2.;
        fWeights[1] = 0.223381589678011/2.;
        fWeights[2] = 0.223381589678011/2.;
        fWeights[3] = 0.109951743655322/2.;
        fWeights[4] = 0.109951743655322/2.;
        fWeights[5] = 0.109951743655322/2.;
        break;
    case 5:
        fPoints.resize(7,2);
        fWeights.resize(7);

        fPoints(0,0) = 0.333333333333333;
        fPoints(0,1) = 0.333333333333333;
        fPoints(1,0) = 0.470142064105115;
        fPoints(1,1) = 0.059715871789770;
        fPoints(2,0) = 0.470142064105115;
        fPoints(2,1) = 0.470142064105115;
        fPoints(3,0) = 0.059715871789770;
        fPoints(3,1) = 0.470142064105115;
        fPoints(4,0) = 0.101286507323456;
        fPoints(4,1) = 0.797426985353087;
        fPoints(5,0) = 0.101286507323456;
        fPoints(5,1) = 0.101286507323456;
        fPoints(6,0) = 0.797426985353087;
        fPoints(6,1) = 0.101286507323456;

        fWeights[0] = 0.225/2.;
        fWeights[1] = 0.132394152788506/2.;
        fWeights[2] = 0.132394152788506/2.;
        fWeights[3] = 0.132394152788506/2.;
        fWeights[4] = 0.125939180544827/2.;
        fWeights[5] = 0.125939180544827/2.;
        fWeights[6] = 0.125939180544827/2.;

        break;
        
   default:
    DebugStop();
    break;
   }
}
