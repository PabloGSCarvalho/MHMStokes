

#include <cmath>
#include <set>

#include <iostream>
#include <fstream>
#include <string>
#include "DarcyPTest.h"
#include "StokesTest.h"
#include "BrinkmanTest.h"
#include "CoupledTest.h"

#include "TPZCouplingDSMaterial.h"
#include "TPZStokesMaterial.h"
#include "TPZDarcyPMaterial.h"

#include <pzlog.h>

#define TEST_DOMAINS
//#define APP_CURVE

//HDivPiola = 1;
const int SpaceHDiv = 1; //Velocidade em subespaço de H(div)
const int SpaceContinuous = 2; //Velocidade em subespaço de [H1]ˆ2
const int SpaceDiscontinuous = 3; //Velociadade em subespaço de H(Ph) - Ph: partição
const REAL Pi=M_PI;

const REAL visco=1., permeability=1., theta=-1.; //Coeficientes: viscosidade, permeabilidade, fator simetria

bool DarcyDomain = false, StokesDomain = false , BrinkmanDomain = true, CoupledDomain = false;

int main(int argc, char *argv[])
{
//    gRefDBase.InitializeAllUniformRefPatterns();

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    //Dados do problema:

    REAL hx=2.,hy=2.; //Dimensões em x e y do domínio
    //double hx=Pi,hy=2.;
    int h_level = 2;
    int nx=h_level+1 ,ny=h_level+1; //Número de nos em x  y
    int pOrder = 2; //Ordem polinomial de aproximação
    
    TPZVec<REAL> h_s(3,0);
    h_s[0]=2.,h_s[1]=2.,h_s[2]=2.; //Dimensões em x e y do domínio

    if (DarcyDomain) {
        DarcyPTest * Test1 = new DarcyPTest();
        Test1->Run(SpaceHDiv, pOrder, nx, ny, hx, hy,visco,permeability,theta);
    }
    else if (StokesDomain)
    {
        pOrder = 3;

        TPZVec<STATE> S0(13,0.);
        S0[0]=0.0000001,S0[1]=1.,S0[2]=3.,S0[3]=5.,S0[4]=10.,S0[5]=15.,S0[6]=20.,S0[7]=25.,S0[8]=30.,S0[9]=35.,S0[10]=40.,S0[11]=45.,S0[12]=50.;
        
        hx=2., hy=2.;

       for (int it=0; it<=12; it++) {
           //h_level = pow(2., 2+it);
           h_level = 2;
           //Coeficiente estabilização (Stokes)
           STATE hE=hx/h_level;
           STATE s0=S0[it];
           STATE sigma=s0*(pOrder*pOrder)/hE;


           nx=h_level+1 ,ny=h_level+1;
           hE=hx/h_level;
           sigma=s0*(pOrder*pOrder)/hE;
           StokesTest  * Test1 = new StokesTest();
           Test1->Run(SpaceHDiv, pOrder, nx, ny, hx, hy,visco,theta,sigma);
           //h_level = h_level*2;
       }
        
    }
    else if (BrinkmanDomain)
    {
        pOrder = 2;
        hx=2.,hy=2.;
        
        TPZVec<STATE> S0(13,0.);
        S0[0]=0.0000001,S0[1]=1.,S0[2]=3.,S0[3]=5.,S0[4]=10.,S0[5]=15.,S0[6]=20.,S0[7]=25.,S0[8]=30.,S0[9]=35.,S0[10]=40.,S0[11]=45.,S0[12]=50.;
        
        
        for (int it=0; it<=0; it++) {
            //h_level = pow(2., 2+it);
            h_level = 4;
            
            //Coeficiente estabilização (Stokes)
            STATE hE=hx/h_level;
            STATE s0=4.;
            STATE sigma=s0*(pOrder*pOrder)/hE;
            
            
            nx=h_level+1 ,ny=h_level+1;
            hE=hx/h_level;
            sigma=s0*(pOrder*pOrder)/hE;
            
            BrinkmanTest  * Test2 = new BrinkmanTest();
            //Test2->SetTriangularMesh();
            //Test2->SetHdivPlus();
            Test2->SetBrinkmanCoef(1.);
            Test2->SetViscosity(0.);
            Test2->Run(SpaceHDiv, pOrder, nx, ny, hx, hy, theta, sigma);

            //            BrinkmanTest  * Test1 = new BrinkmanTest();
            //            Test1->SetTriangularMesh();
            //            Test1->SetFullHdiv();
            //            Test1->Run(SpaceHDiv, pOrder, nx, ny, hx, hy,visc,theta,sigma);
            //
            //            BrinkmanTest  * Test3 = new BrinkmanTest();
            //            Test3->SetTriangularMesh();
            //            Test3->Run(SpaceDiscontinuous, pOrder, nx, ny, hx, hy,visc,theta,sigma);

   
            //h_level = h_level*2;
        }
        
    }
    else  if(CoupledDomain)
    {
        int h_level = 64;
        
        //double hx=1.,hy=1.; //Dimensões em x e y do domínio
        double hx=Pi,hy=2.; //Dimensões em x e y do domínio (acoplamento)
        int nelx=h_level, nely=h_level; //Número de elementos em x e y
        int nx=nelx+1 ,ny=nely+1; //Número de nos em x  y
        int pOrder = 2; //Ordem polinomial de aproximação
        STATE hE=hx/h_level;
        STATE s0=12.;
        STATE sigma=s0*(pOrder*pOrder)/hE;
        
        CoupledTest  * Test3 = new CoupledTest();
        Test3->Run(SpaceHDiv, pOrder, nx, ny, hx, hy,visco,permeability,theta,sigma);
    }
    
    return 0;
}