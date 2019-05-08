/*
 *  MHMBrinkmanTest.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 28/07/2017.
 *  Copyright 2017 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __PZ__MHMBrinkmanTest__
#define __PZ__MHMBrinkmanTest__

#include <cmath>
#include <set>

#include <iostream>
#include <fstream>
#include <string>
#include "pzgmesh.h"
#include "pzstack.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"
#include "TPZMHMBrinkmanMaterial.h"
#include "TPZMHMBrinkmanBC.h"
#include "TPZMultiphysicsCompMesh.h"

#include <pzgeoel.h>
#include "pzgeoelbc.h"
#include "pzfmatrix.h"
#include "pzbstrmatrix.h"
#include <TPZGeoElement.h>
#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzmat1dlin.h"
#include "pzmat2dlin.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZParSkylineStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZGeoLinear.h"
#include "tpzgeoelrefpattern.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzanalysis.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSpStructMatrix.h"

using namespace std;
using namespace pzshape;

class MHMBrinkmanTest{
private:
    
    int fdim; //Dimensão do problema
    int fmatID; //Materia do elemento volumétrico
    
    //Materiais das condições de contorno
    int fmatBCbott;
    int fmatBCtop;
    int fmatBCleft;
    int fmatBCright;
    
    //Material do elemento de interface
    int fmatLambda; //Multiplier material
    int fmatLambdaBC;
    
    int fmatLambdaBC_bott;
    int fmatLambdaBC_top;
    int fmatLambdaBC_left;
    int fmatLambdaBC_right;
    
    
    int fmatInterfaceLeft;
    int fmatInterfaceRight;
    int fmatWrap;
    
    //Materiais das condições de contorno (elementos de interface)
    int fmatIntBCbott;
    int fmatIntBCtop;
    int fmatIntBCleft;
    int fmatIntBCright;
    
    //Materia de um ponto
    int fmatPoint;
    
    //Condições de contorno do problema
    int fdirichlet;
    int fneumann;
    int fpenetration;
    int fpointtype;
    int fdirichletvar;
    
    
    int fquadmat1; //Parte inferior do quadrado
    int fquadmat2; //Parte superior do quadrado
    int fquadmat3; //Material de interface
    
    STATE fviscosity;
    STATE fpermeability;
    STATE ftheta;
    
    int fSpaceV;
    
    REAL fphi_r;
    
    bool f_is_hdivFull;
    
    bool f_hdivPlus;
    
    bool fTriang;
    
    TPZVec<TPZCompMesh * > f_mesh_vector;
    
public:

    MHMBrinkmanTest();
    
    ~MHMBrinkmanTest();
    
    void Run(int Space, int pOrder, int nx, int ny, double hx, double hy, STATE visco, STATE theta, STATE sigma);
    
    /*  Malhas geometricas */
    TPZGeoMesh *CreateGMesh(int nx, int ny, double hx, double hy);
    
    //   TPZGeoMesh *GMeshDeformed(int dim, bool ftriang, int ndiv);
    
    void ChangeExternalOrderConnects(TPZCompMesh *mesh, int addToOrder);
    /* Malhas computacionais */
    
    TPZCompEl *CreateInterfaceEl(TPZGeoEl *gel,TPZCompMesh &mesh, int64_t &index);
    
    TPZCompMesh *CMesh_v(TPZGeoMesh *gmesh, int Space, int pOrder);
    TPZCompMesh *CMesh_p(TPZGeoMesh *gmesh, int Space, int pOrder);
    //TPZCompMesh *CMesh_St(TPZGeoMesh *gmesh, int Space, int pOrder);
    TPZMultiphysicsCompMesh *CMesh_m(TPZGeoMesh *gmesh, int Space, int pOrder, STATE visco, STATE theta, STATE sigma);
    
    void SetHdivPlus(){
        f_hdivPlus = true;
    };
    
    
    void SetFullHdiv(){
        f_is_hdivFull = true;
    };

    void SetTriangularMesh(){
        fTriang = true;
    };
    
    //solucao exata
    static void Sol_exact(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);
    
    //lado direito da equacao
    static void F_source(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE>& gradu);
    
    // static void AddMultiphysicsInterfaces(TPZCompMesh &cmesh, int matfrom, int mattarget);
    void AddMultiphysicsInterfaces(TPZMultiphysicsCompMesh &cmesh);
    
    void AddMultiphysicsInterfaces(TPZMultiphysicsCompMesh &cmesh, int matfrom, int mattarget);
    
    void AddMultiphysicsInterfacesLeftNRight(TPZMultiphysicsCompMesh &cmesh, int matfrom);
    
    // Rotate function
    void Rotate(TPZVec<REAL> &co, TPZVec<REAL> &co_r, bool rotate);
    
    // Insere interfaces na malha multifísica
    void InsertInterfaces(TPZMultiphysicsCompMesh *cmesh);
};


#endif 
