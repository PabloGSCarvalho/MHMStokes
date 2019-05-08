/*
 *  MHMBrinkmanTest.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 28/07/2017.
 *  Copyright 2017 __MyCompanyName__. All rights reserved.
 *
 */

#include "MHMBrinkmanTest.h"
#include "pzcheckgeom.h"
#include "pzstack.h"
#include "TPZParSkylineStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "TPZGmshReader.h"
#include "TPZInterfaceInsertion.h"
#include "TPZVecL2.h"

using namespace std;

const REAL Pi=M_PI;

const REAL phi_r = 0.;

MHMBrinkmanTest::MHMBrinkmanTest()
{
    
    fdim=2; //Dimensão do problema
    fmatID=1; //Materia do elemento volumétrico
    
    //Materiais das condições de contorno
    fmatBCbott=-1;
    fmatBCtop=-2;
    fmatBCleft=-3;
    fmatBCright=-4;
    
    //Material do elemento de interface
    fmatLambda=4; // Multiplier material
    fmatLambdaBC=3;
    
    fmatLambdaBC_bott=11;
    fmatLambdaBC_top=12;
    fmatLambdaBC_left=13;
    fmatLambdaBC_right=14;

    fmatInterfaceLeft=5;
    fmatInterfaceRight=6;
    fmatWrap = 7;
    
    //Materiais das condições de contorno (elementos de interface)
    fmatIntBCbott=-11;
    fmatIntBCtop=-12;
    fmatIntBCleft=-13;
    fmatIntBCright=-14;
    
    //Materia de um ponto
    fmatPoint=-5;
    
    //Condições de contorno do problema
    fdirichlet=0;
    fneumann=1;
    fpenetration=2;
    fpointtype=5;
    fdirichletvar=4;
    
    
    fquadmat1=1; //Parte inferior do quadrado
    fquadmat2=2; //Parte superior do quadrado
    fquadmat3=3; //Material de interface
    
    fviscosity=1.;
    fpermeability=1.;
    ftheta=-1.;
    
    fSpaceV=0;
    
    fphi_r=0;
    
    f_is_hdivFull = false;
    
    f_hdivPlus = false;
    
    fTriang = false;
    
    f_mesh_vector.resize(2);
}

MHMBrinkmanTest::~MHMBrinkmanTest()
{
    
}

void MHMBrinkmanTest::Run(int Space, int pOrder, int nx, int ny, double hx, double hy, STATE visco, STATE theta, STATE sigma)
{
    
    
    //Gerando malha geométrica:
    fSpaceV = Space;
    TPZGeoMesh *gmesh = CreateGMesh(nx, ny, hx, hy); //Função para criar a malha geometrica
    
#ifdef PZDEBUG
    std::ofstream fileg("MalhaGeo.txt"); //Impressão da malha geométrica (formato txt)
    std::ofstream filegvtk("MalhaGeo.vtk"); //Impressão da malha geométrica (formato vtk)
    gmesh->Print(fileg);
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk,true);
#endif
    
    //Gerando malha computacional:
    
    int n_mais = 0;
    if (f_hdivPlus) {
        n_mais = 1;
    }
    
    
    TPZCompMesh *cmesh_v = this->CMesh_v(gmesh, Space, pOrder); //Função para criar a malha computacional da velocidade
    TPZCompMesh *cmesh_p = this->CMesh_p(gmesh, Space, pOrder+n_mais); //Função para criar a malha computacional da pressão
    //TPZCompMesh *cmesh_St = this->CMesh_St(gmesh, Space, pOrder); //Função para criar a malha computacional da pressão
    
    ChangeExternalOrderConnects(cmesh_v,n_mais);
    
    f_mesh_vector[0]=cmesh_v;
    f_mesh_vector[1]=cmesh_p;
    
    TPZMultiphysicsCompMesh *cmesh_m = this->CMesh_m(gmesh, Space, pOrder, visco, theta, sigma); //Função para criar a malha computacional multifísica
    
#ifdef PZDEBUG
    {
        std::ofstream filecv("MalhaC_v.txt"); //Impressão da malha computacional da velocidade (formato txt)
        std::ofstream filecp("MalhaC_p.txt"); //Impressão da malha computacional da pressão (formato txt)
        std::ofstream filecSt("MalhaC_St.txt"); //Impressão da malha computacional da pressão (formato txt)
        cmesh_v->Print(filecv);
        cmesh_p->Print(filecp);
        //cmesh_St->Print(filecSt);
        
        std::ofstream filecm("MalhaC_m.txt"); //Impressão da malha computacional multifísica (formato txt)
        cmesh_m->Print(filecm);
    }
#endif

    
    cmesh_m->LoadReferences();
    InsertInterfaces(cmesh_m);
    

//    AddMultiphysicsInterfacesLeftNRight(*cmesh_m,fmatLambda); // Rever isto aqui
//    AddMultiphysicsInterfaces(*cmesh_m,fmatIntBCbott,fmatBCbott);
//    AddMultiphysicsInterfaces(*cmesh_m,fmatIntBCtop,fmatBCtop);
//    AddMultiphysicsInterfaces(*cmesh_m,fmatIntBCleft,fmatBCleft);
//    AddMultiphysicsInterfaces(*cmesh_m,fmatIntBCright,fmatBCright);
    
//    AddMultiphysicsInterfaces(*cmesh_m);
    
#ifdef PZDEBUG
    std::ofstream fileg1("MalhaGeo.txt"); //Impressão da malha geométrica (formato txt)
    std::ofstream filegvtk1("MalhaGeo.vtk"); //Impressão da malha geométrica (formato vtk)
    gmesh->Print(fileg1);
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk1,true);
    
    std::ofstream filecm("MalhaC_m.txt"); //Impressão da malha computacional multifísica (formato txt)
    cmesh_m->Print(filecm);
#endif
    
    //Resolvendo o Sistema:
    int numthreads = 0;
    
    bool optimizeBandwidth = true; //Impede a renumeração das equacoes do problema (para obter o mesmo resultado do Oden)
    TPZAnalysis an(cmesh_m, optimizeBandwidth); //Cria objeto de análise que gerenciará a analise do problema
    
//            TPZSpStructMatrix struct_mat(cmesh_m);
//            struct_mat.SetNumThreads(numthreads);
//            an.SetStructuralMatrix(struct_mat);
    
    
    //TPZParSkylineStructMatrix matskl(cmesh_m, numthreads);

    TPZSkylineNSymStructMatrix matskl(cmesh_m); //OK para Hdiv
    matskl.SetNumThreads(numthreads);
    an.SetStructuralMatrix(matskl);
//    
//    if (Space==1) {
//        TPZFStructMatrix matsklD(cmesh_m); //caso nao simetrico *** //OK para discont.
//        matsklD.SetNumThreads(numthreads);
//        an.SetStructuralMatrix(matsklD);
//    }


    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    an.SetSolver(step);
    
    
    
    std::cout << "Assemble matrix with NDoF = " << cmesh_m->NEquations() << std::endl;
    
    an.Assemble();//Assembla a matriz de rigidez (e o vetor de carga) global
     {
         std::ofstream filestiff("stiffness_before.txt");
         an.Solver().Matrix()->Print("K1 = ",filestiff,EMathematicaInput);
         
     }
    std::cout << "Solving Matrix " << std::endl;
    an.Solve();
    
#ifdef PZDEBUG
    //Imprimir Matriz de rigidez Global:
    {
        std::ofstream filestiff("stiffness.txt");
        an.Solver().Matrix()->Print("K1 = ",filestiff,EMathematicaInput);
        
        std::ofstream filerhs("rhs.txt");
        an.Rhs().Print("R = ",filerhs,EMathematicaInput);
    }
#endif
    
#ifdef PZDEBUG
    //Imprimindo vetor solução:
    {
        TPZFMatrix<STATE> solucao=cmesh_m->Solution();//Pegando o vetor de solução, alphaj
        std::ofstream solout("sol.txt");
        solucao.Print("Sol",solout,EMathematicaInput);//Imprime na formatação do Mathematica
        
        std::ofstream fileAlpha("alpha.txt");
        an.Solution().Print("Alpha = ",fileAlpha,EMathematicaInput);
    }
#endif
    
    //Calculo do erro
    std::cout << "Comuting Error " << std::endl;
    TPZManVector<REAL,6> Errors;
    ofstream ErroOut("Error_Brinkman.txt", std::ofstream::app);
    an.SetExact(Sol_exact);
    an.PostProcessError(Errors,false);
    
    ErroOut <<"Sigma = "<< sigma/(pOrder*pOrder*(nx-1)) << "  //  Ordem = "<< pOrder << "  //  Tamanho da malha = "<< nx-1 <<" x "<< ny-1 << std::endl;
    ErroOut <<" " << std::endl;
    //ErroOut <<"Norma H1/HDiv - V = "<< Errors[0] << std::endl;
    ErroOut <<"Norma L2 - V = "<< Errors[1] << std::endl;
    ErroOut <<"Semi-norma H1/Hdiv - V = "<< Errors[2] << std::endl;
    ErroOut <<"Norma L2 - P = "<< Errors[4] << std::endl;
    ErroOut <<"-------------" << std::endl;
    ErroOut.flush();
    
    //Pós-processamento (paraview):
    std::cout << "Post Processing " << std::endl;
    std::string plotfile("Brinkman.vtk");
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("P");
    vecnames.Push("V");
    vecnames.Push("f");
    vecnames.Push("V_exact");
    scalnames.Push("P_exact");
    scalnames.Push("Div");


    int postProcessResolution = 3; //  keep low as possible

    int dim = gmesh->Dimension();
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(postProcessResolution,dim);
    
    std::cout << "FINISHED!" << std::endl;
    
}


void MHMBrinkmanTest::Rotate(TPZVec<REAL> &co, TPZVec<REAL> &co_r, bool rotate){
    
    if (rotate==true) {
        //rotação +
        co_r[0] = co[0]*cos(phi_r) - co[1]*sin(phi_r);
        co_r[1] = co[0]*sin(phi_r) + co[1]*cos(phi_r);
        
    }else{
        
        co_r[0] = co[0]*cos(phi_r) + co[1]*sin(phi_r);
        co_r[1] = - co[0]*sin(phi_r) + co[1]*cos(phi_r);
        
    }
    
    
}

TPZGeoMesh *MHMBrinkmanTest::CreateGMesh(int nx, int ny, double hx, double hy)
{
    
    
    if(fTriang){
        int i,j;
        int64_t id, index;
        
        
        //Criando malha geométrica, nós e elementos.
        //Inserindo nós e elementos no objeto malha:
        
        TPZGeoMesh *gmesh = new TPZGeoMesh();
        gmesh->SetDimension(2);
        
        //Vetor auxiliar para armazenar coordenadas:
        
        TPZVec<REAL> coord (3,0.), coord_r(3,0.);
        
        //Inicialização dos nós:
        
        for(i = 0; i < ny; i++){
            for(j = 0; j < nx; j++){
                id = i*nx + j;
                coord[0] = (j)*hx/(nx - 1);
                coord[1] = -1.+(i)*hy/(ny - 1);
                //using the same coordinate x for z
                coord[2] = 0.;
                //cout << coord << endl;
                
                //rottação phi
                Rotate(coord, coord_r, true);
                
                //Get the index in the mesh nodes vector for the new node
                index = gmesh->NodeVec().AllocateNewElement();
                //Set the value of the node in the mesh nodes vector
                gmesh->NodeVec()[index] = TPZGeoNode(id,coord_r,*gmesh);
            }
        }
        
        
        //Ponto 1
        TPZVec<int64_t> pointtopology(1);
        pointtopology[0] = nx-1;
        
        gmesh->CreateGeoElement(EPoint,pointtopology,fmatPoint,id);
        
        
        //Vetor auxiliar para armazenar as conecções entre elementos:
        
        TPZVec <int64_t> connectD(3,0);
        TPZVec <int64_t> connectU(3,0);
        
        
        //Conectividade dos elementos:
        
        for(i = 0; i < (ny - 1); i++){
            for(j = 0; j < (nx - 1); j++){
                index = (i)*(nx - 1)+ (j);
                connectD[0] = (i)*ny + (j);
                connectD[1] = connectD[0]+1;
                connectD[2] = connectD[0]+nx+1;
                gmesh->CreateGeoElement(ETriangle,connectD,fmatID,id);
                
                connectU[0] = connectD[2];
                connectU[1] = connectD[2]-1;
                connectU[2] = connectD[0];
                gmesh->CreateGeoElement(ETriangle,connectU,fmatID,id);
                
                //   std::cout<<connectD<<std::endl;
                //   std::cout<<connectU<<std::endl;
                
                id++;
            }
        }
        
        
        //Gerando informação da vizinhança:
        
        gmesh->BuildConnectivity();
        
        {
            TPZCheckGeom check(gmesh);
            check.CheckUniqueId();
        }
        int64_t el, numelements = gmesh->NElements();
        
        TPZManVector <int64_t> TopolPlate(4);
        
        for (el=0; el<numelements; el++)
        {
            int64_t totalnodes = gmesh->ElementVec()[el]->NNodes();
            TPZGeoEl *plate = gmesh->ElementVec()[el];
            for (int i=0; i<4; i++){
                TopolPlate[i] = plate->NodeIndex(i);
            }
            
            //Colocando as condicoes de contorno:
            TPZManVector <TPZGeoNode> Nodefinder(totalnodes);
            TPZManVector <REAL,3> nodecoord(3,0.),nodecoord_r(3,0.);
            
            //Na face x = 1
            TPZVec<int64_t> ncoordzbottVec(0); int64_t sizeOfbottVec = 0;
            TPZVec<int64_t> ncoordztopVec(0); int64_t sizeOftopVec = 0;
            TPZVec<int64_t> ncoordzleftVec(0); int64_t sizeOfleftVec = 0;
            TPZVec<int64_t> ncoordzrightVec(0); int64_t sizeOfrightVec = 0;
            
            
            for (int64_t i = 0; i < totalnodes; i++)
            {
                Nodefinder[i] = gmesh->NodeVec()[TopolPlate[i]];
                Nodefinder[i].GetCoordinates(nodecoord_r);
                int id_node = Nodefinder[i].Id();
                
                for (int64_t j = 0; j < ny; j++){
                    
                    
                    if (id_node==j)
                    {
                        sizeOfbottVec++;
                        ncoordzbottVec.Resize(sizeOfbottVec);
                        ncoordzbottVec[sizeOfbottVec-1] = TopolPlate[i];
                    }
                    if (id_node==j+nx*(nx-1))
                    {
                        sizeOftopVec++;
                        ncoordztopVec.Resize(sizeOftopVec);
                        ncoordztopVec[sizeOftopVec-1] = TopolPlate[i];
                    }
                    
                    
                    if (id_node==j*nx)
                    {
                        sizeOfleftVec++;
                        ncoordzleftVec.Resize(sizeOfleftVec);
                        ncoordzleftVec[sizeOfleftVec-1] = TopolPlate[i];
                    }
                    if (id_node==(j+1)*nx-1)
                    {
                        sizeOfrightVec++;
                        ncoordzrightVec.Resize(sizeOfrightVec);
                        ncoordzrightVec[sizeOfrightVec-1] = TopolPlate[i];
                    }
                    
                }
            }
            
            
            if (sizeOfbottVec == 2) {
                int sidesbott = plate->WhichSide(ncoordzbottVec);
                TPZGeoElSide platesidebott(plate, sidesbott);
                TPZGeoElBC(platesidebott,fmatBCbott);
                if(fSpaceV!=2){
                    TPZGeoElBC(platesidebott,fmatIntBCbott);
                }
            }
            
            if (sizeOftopVec == 2) {
                int sidestop = plate->WhichSide(ncoordztopVec);
                TPZGeoElSide platesidetop(plate, sidestop);
                TPZGeoElBC(platesidetop,fmatBCtop);
                if(fSpaceV!=2){
                    TPZGeoElBC(platesidetop,fmatIntBCtop);
                }
            }
            
            if (sizeOfleftVec == 2) {
                int sidesleft = plate->WhichSide(ncoordzleftVec);
                TPZGeoElSide platesideleft(plate, sidesleft);
                TPZGeoElBC(platesideleft,fmatBCleft);
                if(fSpaceV!=2){
                    TPZGeoElBC(platesideleft,fmatIntBCleft);
                }
            }
            
            if (sizeOfrightVec == 2) {
                int sidesright = plate->WhichSide(ncoordzrightVec);
                TPZGeoElSide platesideright(plate, sidesright);
                TPZGeoElBC(platesideright,fmatBCright);
                if(fSpaceV!=2){
                    TPZGeoElBC(platesideright,fmatIntBCright);
                }
            }
            
            
            ncoordzbottVec.Resize(0);
            sizeOfbottVec = 0;
            ncoordztopVec.Resize(0);
            sizeOftopVec = 0;
            ncoordzleftVec.Resize(0);
            sizeOfleftVec = 0;
            ncoordzrightVec.Resize(0);
            sizeOfrightVec = 0;
            
        }
        
        
        //Criando 1D material (Lambda multiplier):
        
        if(fSpaceV!=2){
            TPZVec<int64_t> nodint(2);
            for(i = 0; i < (ny - 1); i++){
                for(j = 0; j <= (nx - 1); j++){
                    if(j>0&&j<(nx-1)){
                        nodint[0]=j+nx*i;
                        nodint[1]=j+nx*(i+1);
                        gmesh->CreateGeoElement(EOned, nodint, fmatLambda, index); //Criando elemento de interface (GeoElement)
                        
                    }
                    
                    
                    if(i>0&&j<(ny-1)){
                        nodint[0]=j+ny*i;
                        nodint[1]=j+ny*i+1;
                        gmesh->CreateGeoElement(EOned, nodint, fmatLambda, index); //Criando elemento de interface (GeoElement)
                        
                    }
                    
                    if(j<(nx-1)&&i<(ny-1)){
                        nodint[0]=j+nx*i;
                        nodint[1]=j+nx*i+nx+1;
                        gmesh->CreateGeoElement(EOned, nodint, fmatLambda, index); //Criando elemento de interface (GeoElement)
                        
                    }
                    
                }
            }
        }
        
        //new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (nodind3,matInterface,*gmesh); //Criando elemento de interface (RefPattern)
        id++;
        
        gmesh->AddInterfaceMaterial(fquadmat1, fquadmat2, fquadmat3);
        gmesh->AddInterfaceMaterial(fquadmat2, fquadmat1, fquadmat3);
        
        TPZCheckGeom check(gmesh);
        check.CheckUniqueId();
        
        gmesh->BuildConnectivity();
        
        //Impressão da malha geométrica:
        
        ofstream bf("before.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
        return gmesh;
        
    }else{
        
        int i,j;
        int64_t id, index;
        
        
        //Criando malha geométrica, nós e elementos.
        //Inserindo nós e elementos no objeto malha:
        
        TPZGeoMesh *gmesh = new TPZGeoMesh();
        gmesh->SetDimension(2);
        
        //Vetor auxiliar para armazenar coordenadas:
        
        TPZVec <REAL> coord (3,0.);
        
        
        //Inicialização dos nós:
        
        for(i = 0; i < ny; i++){
            for(j = 0; j < nx; j++){
                id = i*nx + j;
                coord[0] = (j)*hx/(nx - 1);
                coord[1] = (i)*hy/(ny - 1);
                //using the same coordinate x for z
                coord[2] = 0.;
                //cout << coord << endl;
                //Get the index in the mesh nodes vector for the new node
                index = gmesh->NodeVec().AllocateNewElement();
                //Set the value of the node in the mesh nodes vector
                gmesh->NodeVec()[index] = TPZGeoNode(id,coord,*gmesh);
            }
        }
        
        //Ponto 1
//        TPZVec<int64_t> pointtopology(1);
//        pointtopology[0] = 0;
//
//        gmesh->CreateGeoElement(EPoint,pointtopology,fmatPoint,id);

        //Ponto 2
//        TPZVec<int64_t> pointtopology2(1);
//        pointtopology2[0] = nx-1;
//
//        gmesh->CreateGeoElement(EPoint,pointtopology2,fmatPoint,id);

//
//        //Ponto 3
//        TPZVec<int64_t> pointtopology3(1);
//        pointtopology3[0] = 1;
//
//        gmesh->CreateGeoElement(EPoint,pointtopology3,fmatPoint,id);
//
        
        //Vetor auxiliar para armazenar as conecções entre elementos:
        
        TPZVec <int64_t> connect(4,0);
        
        
        //Conectividade dos elementos:
        
        for(i = 0; i < (ny - 1) ; i++){
            for(j = 0; j < (nx - 1); j++){
                index = (i)*(nx - 1)+ (j);
                connect[0] = (i)*ny + (j);
                connect[1] = connect[0]+1;
                connect[2] = connect[1]+(nx);
                connect[3] = connect[0]+(nx);
                gmesh->CreateGeoElement(EQuadrilateral,connect,fmatID,id);
            }
        }
        
        
        //Gerando informação da vizinhança:
        
        gmesh->BuildConnectivity();
        
        {
            TPZCheckGeom check(gmesh);
            check.CheckUniqueId();
        }
        int64_t el, numelements = gmesh->NElements();
        
        TPZManVector <int64_t> TopolPlate(4);
        
        for (el=0; el<numelements; el++)
        {
            int64_t totalnodes = gmesh->ElementVec()[el]->NNodes();
            TPZGeoEl *plate = gmesh->ElementVec()[el];
            for (int i=0; i<4; i++){
                TopolPlate[i] = plate->NodeIndex(i);
            }
            
            //Colocando as condicoes de contorno:
            TPZManVector <TPZGeoNode> Nodefinder(totalnodes);
            TPZManVector <REAL,3> nodecoord(3);
            
            //Na face x = 1
            TPZVec<int64_t> ncoordzbottVec(0); int64_t sizeOfbottVec = 0;
            TPZVec<int64_t> ncoordztopVec(0); int64_t sizeOftopVec = 0;
            TPZVec<int64_t> ncoordzleftVec(0); int64_t sizeOfleftVec = 0;
            TPZVec<int64_t> ncoordzrightVec(0); int64_t sizeOfrightVec = 0;
            
            for (int64_t i = 0; i < totalnodes; i++)
            {
                Nodefinder[i] = gmesh->NodeVec()[TopolPlate[i]];
                Nodefinder[i].GetCoordinates(nodecoord);
                if (nodecoord[1] == 0.)
                {
                    sizeOfbottVec++;
                    ncoordzbottVec.Resize(sizeOfbottVec);
                    ncoordzbottVec[sizeOfbottVec-1] = TopolPlate[i];
                }
                if (nodecoord[1] == 0.+hy)
                {
                    sizeOftopVec++;
                    ncoordztopVec.Resize(sizeOftopVec);
                    ncoordztopVec[sizeOftopVec-1] = TopolPlate[i];
                }
                if (nodecoord[0] == 0.)
                {
                    sizeOfleftVec++;
                    ncoordzleftVec.Resize(sizeOfleftVec);
                    ncoordzleftVec[sizeOfleftVec-1] = TopolPlate[i];
                }
                if (nodecoord[0] == hx)
                {
                    sizeOfrightVec++;
                    ncoordzrightVec.Resize(sizeOfrightVec);
                    ncoordzrightVec[sizeOfrightVec-1] = TopolPlate[i];
                }
            }
            
            if (sizeOfbottVec == 2) {
                int sidesbott = plate->WhichSide(ncoordzbottVec);
                TPZGeoElSide platesidebott(plate, sidesbott);
                TPZGeoElBC(platesidebott,fmatBCbott);
//                if(fSpaceV!=2){
//                    TPZGeoElBC(platesidebott,fmatIntBCbott);
//                }
            }
            
            if (sizeOftopVec == 2) {
                int sidestop = plate->WhichSide(ncoordztopVec);
                TPZGeoElSide platesidetop(plate, sidestop);
                TPZGeoElBC(platesidetop,fmatBCtop);
//                if(fSpaceV!=2){
//                    TPZGeoElBC(platesidetop,fmatIntBCtop);
//                }
            }
            
            if (sizeOfleftVec == 2) {
                int sidesleft = plate->WhichSide(ncoordzleftVec);
                TPZGeoElSide platesideleft(plate, sidesleft);
                TPZGeoElBC(platesideleft,fmatBCleft);
//                if(fSpaceV!=2){
//                    TPZGeoElBC(platesideleft,fmatIntBCleft);
//                }
            }
            
            if (sizeOfrightVec == 2) {
                int sidesright = plate->WhichSide(ncoordzrightVec);
                TPZGeoElSide platesideright(plate, sidesright);
                TPZGeoElBC(platesideright,fmatBCright);
//                if(fSpaceV!=2){
//                    TPZGeoElBC(platesideright,fmatIntBCright);
//                }
            }
            
            
            ncoordzbottVec.Resize(0);
            sizeOfbottVec = 0;
            ncoordztopVec.Resize(0);
            sizeOftopVec = 0;
            ncoordzleftVec.Resize(0);
            sizeOfleftVec = 0;
            ncoordzrightVec.Resize(0);
            sizeOfrightVec = 0;
            
        }
        
        //Criando 1D material
        TPZVec<int64_t> nodint(2);
        
        // Criando elementos 1D internos, para tração tangencial (Lambda multiplier)
        if(fSpaceV!=2){
            for(i = 0; i < (ny - 1); i++){
                for(j = 0; j < (nx - 1); j++){
                    if(j>0&&j<(nx-1)){
                        nodint[0]=j+nx*i;
                        nodint[1]=j+nx*(i+1);
                        gmesh->CreateGeoElement(EOned, nodint, fmatLambda, index); //Criando elemento de interface (GeoElement)
                        
                    }
                    if(i>0&&j<(ny-1)){
                        nodint[0]=j+ny*i;
                        nodint[1]=j+ny*i+1;
                        gmesh->CreateGeoElement(EOned, nodint, fmatLambda, index); //Criando elemento de interface (GeoElement)
                        
                    }
                }
            }
            
        }
        
        // Criando elementos 1D externos, para tração tangencial (Lambda multiplier BC)
        
        for(i = 0; i < ny; i++){
            for(j = 0; j < nx; j++){
                if ((i==0)&&j<(nx-1)) {
                    nodint[0]=j+nx*i;
                    nodint[1]=j+nx*i+1;
                    gmesh->CreateGeoElement(EOned, nodint, fmatLambdaBC_bott, index); //Criando elemento de interface (GeoElement)
                }
                if ((i==(ny-1))&&j<(nx-1)) {
                    nodint[0]=j+nx*i;
                    nodint[1]=j+nx*i+1;
                    gmesh->CreateGeoElement(EOned, nodint, fmatLambdaBC_top, index); //Criando elemento de interface (GeoElement)
                }


                if((j==0)&&i<(ny-1)) {
                    nodint[0]=j+nx*i;
                    nodint[1]=j+nx*(i+1);
                    gmesh->CreateGeoElement(EOned, nodint, fmatLambdaBC_left, index); //Criando elemento de interface (GeoElement)

                }
                if((j==(nx-1))&&i<(ny-1)) {
                    nodint[0]=j+nx*i;
                    nodint[1]=j+nx*(i+1);
                    gmesh->CreateGeoElement(EOned, nodint, fmatLambdaBC_right, index); //Criando elemento de interface (GeoElement)
                }
            }
        }
        
        
        //new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (nodind3,matInterface,*gmesh); //Criando elemento de interface (RefPattern)
        id++;
        
        //   gmesh->AddInterfaceMaterial(fquadmat1, fquadmat2, fquadmat3);
        //   gmesh->AddInterfaceMaterial(fquadmat2, fquadmat1, fquadmat3);
        
        TPZCheckGeom check(gmesh);
        check.CheckUniqueId();
        
        gmesh->BuildConnectivity();
        
        //Impressão da malha geométrica:
        
        ofstream bf("before.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
        return gmesh;
        
        
    }
    

}

TPZCompEl *MHMBrinkmanTest::CreateInterfaceEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    if(!gel->Reference() && gel->NumInterfaces() == 0)
        return new TPZInterfaceElement(mesh,gel,index);
    
    return NULL;
}

void MHMBrinkmanTest::Sol_exact(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol){

    dsol.Resize(3,2);
    sol.Resize(3);
    
    REAL x1 = x[0];
    REAL x2 = x[1];
    
    STATE v_1 = x2;
    STATE v_2 = 0.;
    STATE pressure= 1.;
    
    sol[0]=v_1;
    sol[1]=v_2;
    sol[2]=pressure;
    
    // vx direction
    dsol(0,0)= 0.;
    dsol(0,1)= 1.;
    
    // vy direction
    dsol(1,0)= 0.;
    dsol(1,1)= 0.;
    
    // Gradiente pressão
    dsol(2,0)= 0.;
    dsol(2,1)= 0.;
    
    // General form : : Artigo Botti, Di Pietro, Droniou

//    dsol.Resize(3,2);
//    sol.Resize(3);
//
//    REAL x1 = x[0];
//    REAL x2 = x[1];
//
//    REAL m_v= 1., m_u= 1.0;
//
//    REAL Cf=m_v/m_u;
//
//    STATE v_1 = -exp(-Cf)*sin(x1)*sin(x2)+(1./m_v)*(1.-exp(-Cf))*sin(x1)*sin(x2);
//    STATE v_2 = -exp(-Cf)*cos(x1)*cos(x2)-(1./m_v)*(1.-exp(-Cf))*cos(x1)*cos(x2);
//    STATE pressure= cos(x1)*sin(x2);
//
//    sol[0]=v_1;
//    sol[1]=v_2;
//    sol[2]=pressure;
//
//    // vx direction
//    dsol(0,0)= -exp(-Cf)*cos(x1)*sin(x2)+(1./m_v)*(1.-exp(-Cf))*cos(x1)*sin(x2);
//    dsol(0,1)= exp(-Cf)*cos(x2)*sin(x1)+(1./m_v)*(1.-exp(-Cf))*cos(x2)*sin(x1);
//
//    // vy direction
//    dsol(1,0)= -exp(-Cf)*cos(x2)*sin(x1)+(1./m_v)*(1.-exp(-Cf))*cos(x2)*sin(x1);
//    dsol(1,1)= exp(-Cf)*cos(x1)*sin(x2)+(1./m_v)*(1.-exp(-Cf))*cos(x1)*sin(x2);
//
//    // Gradiente pressão
//    dsol(2,0)= -sin(x1)*sin(x2);
//    dsol(2,1)= cos(x1)*cos(x2);
    
    
    // Brinkman : : Artigo Botti, Di Pietro, Droniou
    
//    dsol.Resize(3,2);
//    sol.Resize(3);
//
//    REAL x1 = x[0];
//    REAL x2 = x[1];
//
//    REAL e = exp(1.);
//
//    STATE v_1 = (1.-2./e)*sin(x1)*sin(x2);
//    STATE v_2 = -1.*cos(x1)*cos(x2);
//    STATE pressure= cos(x1)*sin(x2);
//
//    sol[0]=v_1;
//    sol[1]=v_2;
//    sol[2]=pressure;
//
//    // vx direction
//    dsol(0,0)= (1.-2./e)*cos(x1)*sin(x2);
//    dsol(0,1)= cos(x2)*sin(x1);
//
//    // vy direction
//    dsol(1,0)= (1.-2./e)*cos(x2)*sin(x1);
//    dsol(1,1)= cos(x1)*sin(x2);
//
//    // Gradiente pressão
//    dsol(2,0)= -sin(x1)*sin(x2);
//    dsol(2,1)= cos(x1)*cos(x2);

    // Stokes : : Artigo Botti, Di Pietro, Droniou
    
//    dsol.Resize(3,2);
//    sol.Resize(3);
//
//    REAL x1 = x[0];
//    REAL x2 = x[1];
//
//    REAL e = exp(1.);
//
//    STATE v_1 = -1.*sin(x1)*sin(x2);
//    STATE v_2 = -1.*cos(x1)*cos(x2);
//    STATE pressure= cos(x1)*sin(x2);
//
//    sol[0]=v_1;
//    sol[1]=v_2;
//    sol[2]=pressure;
//
//    // vx direction
//    dsol(0,0)= -1.*cos(x1)*sin(x2);
//    dsol(0,1)= cos(x2)*sin(x1);
//
//    // vy direction
//    dsol(1,0)= -1.*cos(x2)*sin(x1);
//    dsol(1,1)= cos(x1)*sin(x2);
//
//    // Gradiente pressão
//    dsol(2,0)= -sin(x1)*sin(x2);
//    dsol(2,1)= cos(x1)*cos(x2);
    

    // Darcy : : Artigo Botti, Di Pietro, Droniou
    
//        dsol.Resize(3,2);
//        sol.Resize(3);
//
//        REAL x1 = x[0];
//        REAL x2 = x[1];
//
//        STATE v_1 = sin(x1)*sin(x2);
//        STATE v_2 = -1.*cos(x1)*cos(x2);
//        STATE pressure= cos(x1)*sin(x2);
//
//        sol[0]=v_1;
//        sol[1]=v_2;
//        sol[2]=pressure;
//
//        // vx direction
//        dsol(0,0)= cos(x1)*sin(x2);
//        dsol(0,1)= cos(x2)*sin(x1);
//
//        // vy direction
//        dsol(1,0)= cos(x2)*sin(x1);
//        dsol(1,1)= cos(x1)*sin(x2);
//
//        // Gradiente pressão
//        dsol(2,0)= -sin(x1)*sin(x2);
//        dsol(2,1)= cos(x1)*cos(x2);

    
}

void MHMBrinkmanTest::F_source(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE>& gradu){
    
    f.resize(3);
    REAL x1 = x[0];
    REAL x2 = x[1];
    STATE f_1 =0., f_2=0.;

    
    // General form : : Artigo Botti, Di Pietro, Droniou
    
//    REAL m_v= 1., m_u= 1.0;
//
//    REAL Cf=m_v/m_u;
//
//        f_1 = -sin(x1)*sin(x2)-exp(-Cf)*sin(x1)*sin(x2)+(1./m_v)*(1.-exp(-Cf))*sin(x1)*sin(x2)-m_u*(2.*exp(-Cf)*sin(x1)*sin(x2)-(1./m_v)*4.*(1-exp(-Cf))*sin(x1)*sin(x2));
//
//        f_2 = cos(x1)*cos(x2)-exp(-Cf)*cos(x1)*cos(x2)-(1./m_v)*(1.-exp(-Cf))*cos(x1)*cos(x2)-m_u*(2.*exp(-Cf)*cos(x1)*cos(x2)+(1./m_v)*4.*(1-exp(-Cf))*cos(x1)*cos(x2));
//        STATE g_1 = (2./m_v)*(1.-exp(-Cf))*cos(x1)*sin(x2);
//
//        f[0] = f_1; // x direction
//        f[1] = f_2; // y direction
//
//        f[2] = g_1; // g source
    
    
    // Brinkman : : Artigo Botti, Di Pietro, Droniou

//    REAL e = exp(1.);
//
//    f_1 = (-8./e+ 4.)*sin(x1)*sin(x2);
//    f_2 = (2./e- 4.)*cos(x1)*cos(x2);
//    STATE g_1 = 2.*(1.-1./e)*cos(x1)*sin(x2);
//
//    f[0] = f_1; // x direction
//    f[1] = f_2; // y direction
//
//    f[2] = g_1; // g source

    // Stokes : : Artigo Botti, Di Pietro, Droniou
    

//    f_1 = -3.*sin(x1)*sin(x2);
//    f_2 = -1.*cos(x1)*cos(x2);

    f[0] = f_1; // x direction
    f[1] = f_2; // y direction
    f[2] = 0.;

    
    // Darcy : : Artigo Botti, Di Pietro, Droniou
    
//        f_1 = 0.;
//        f_2 = 0.;
//
//        f[0] = f_1; // x direction
//        f[1] = f_2; // y direction
//        f[2] = 2.*cos(x1)*sin(x2);



}

void MHMBrinkmanTest::ChangeExternalOrderConnects(TPZCompMesh *mesh, int addToOrder){
    
    int nEl= mesh-> NElements();
    int dim = mesh->Dimension();
    
    for (int iel=0; iel<nEl; iel++) {
        TPZCompEl *cel = mesh->ElementVec()[iel];
        if (!cel) continue;
        int ncon = cel->NConnects();
        int corder = 0;
        int nshape = 0;
        int nshape2 = 0;
        
        if(cel->Dimension()== dim)
        {
            TPZConnect &conel = cel->Connect(ncon-1);
            corder = conel.Order();
            nshape = conel.NShape();
            
            int neworder = corder + addToOrder;//Aqui = +1
            int64_t cindex = cel->ConnectIndex(ncon-1);
            conel.SetOrder(neworder,cindex);
            
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
            intel->SetPreferredOrder(neworder);
            nshape = intel->NConnectShapeF(ncon-1,neworder);
            
            if(dim==2 && addToOrder==1)
            {
                if(fTriang){
                    nshape2 = (corder + 2)*(corder + 2)-1;
                }else{//Quadrilateral
                    nshape2 = 2*(corder + 1)*(corder + 2);
                }
                if(nshape2!=nshape)
                {
                    DebugStop();
                }
            }
            
            conel.SetNShape(nshape);
            mesh->Block().Set(conel.SequenceNumber(),nshape);
        }
    }
    mesh->CleanUpUnconnectedNodes();
    mesh->ExpandSolution();
}


TPZCompMesh *MHMBrinkmanTest::CMesh_v(TPZGeoMesh *gmesh, int Space, int pOrder)
{
 
    //Criando malha computacional:
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);//Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim);//Insere dimensão do modelo
    
    
    // 1 - Material volumétrico 2D 
    TPZVecL2 *material = new TPZVecL2(fmatID);
    cmesh->InsertMaterialObject(material);
    
    if (Space==1) {
        cmesh->SetAllCreateFunctionsHDiv(); //Criando funções HDIV:
        //cmesh->ApproxSpace().CreateDisconnectedElements(true); //HDIV-Full:
        //Dimensões do material (para HDiv):
        //TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xfin(1,1,0.);
        //material->SetMaterial(xkin, xcin, xfin);
        
    }else{
        DebugStop();
    }
    
    // 1 - Condições de contorno
    TPZFMatrix<STATE> val1(1,1,0.), val2(2,1,0.);
    {
        TPZMaterial * BCond0 = material->CreateBC(material, fmatBCbott, fdirichlet, val1, val2);
        cmesh->InsertMaterialObject(BCond0);
        
        TPZMaterial * BCond1 = material->CreateBC(material, fmatBCtop, fdirichlet, val1, val2);
        cmesh->InsertMaterialObject(BCond1);
        
        TPZMaterial * BCond2 = material->CreateBC(material, fmatBCleft, fdirichlet, val1, val2);
        cmesh->InsertMaterialObject(BCond2);
        
        TPZMaterial * BCond3 = material->CreateBC(material, fmatBCright, fdirichlet, val1, val2);
        cmesh->InsertMaterialObject(BCond3);
    }
    {
//        TPZMaterial * BCond0 = material->CreateBC(material, fmatIntBCbott, fdirichlet, val1, val2);
//        cmesh->InsertMaterialObject(BCond0);
//
//        TPZMaterial * BCond1 = material->CreateBC(material, fmatIntBCtop, fdirichlet, val1, val2);
//        cmesh->InsertMaterialObject(BCond1);
//
//        TPZMaterial * BCond2 = material->CreateBC(material, fmatIntBCleft, fdirichlet, val1, val2);
//        cmesh->InsertMaterialObject(BCond2);
//
//        TPZMaterial * BCond3 = material->CreateBC(material, fmatIntBCright, fdirichlet, val1, val2);
//        cmesh->InsertMaterialObject(BCond3);
        
//        TPZMaterial * BCond4 = material->CreateBC(material, fmatInterLambda, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
//        cmesh->InsertMaterialObject(BCond4); //Insere material na malha
    }
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximacao da malha:
    
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    
        
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    
    return cmesh;
   
    
}


TPZCompMesh *MHMBrinkmanTest::CMesh_p(TPZGeoMesh *gmesh, int Space, int pOrder)
{
    
    //Criando malha computacional:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);

    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim); //Insere dimensão do modelo
    
    cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    
    // 1 - Material volumétrico 2D
    TPZVecL2 *material = new TPZVecL2(fmatID);
    cmesh->InsertMaterialObject(material);

    //Dimensões do material (para H1 e descontínuo):
    TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xbin(1,1,0.), xfin(1,1,0.);
//    material->SetMaterial(xkin, xcin, xfin);
    
    // 2 - Material para tração tangencial 1D

    TPZMat1dLin *matLambda = new TPZMat1dLin(fmatLambda);
    matLambda->SetMaterial(xkin, xcin, xbin, xfin);
    cmesh->InsertMaterialObject(matLambda);
    
    // 3 - Material para tração tangencial nos contornos
    
    TPZMat1dLin *matLambdaBC_bott = new TPZMat1dLin(fmatLambdaBC_bott);
    matLambdaBC_bott->SetMaterial(xkin, xcin, xbin, xfin);
    cmesh->InsertMaterialObject(matLambdaBC_bott);

    TPZMat1dLin *matLambdaBC_top = new TPZMat1dLin(fmatLambdaBC_top);
    matLambdaBC_top->SetMaterial(xkin, xcin, xbin, xfin);
    cmesh->InsertMaterialObject(matLambdaBC_top);

    TPZMat1dLin *matLambdaBC_left = new TPZMat1dLin(fmatLambdaBC_left);
    matLambdaBC_left->SetMaterial(xkin, xcin, xbin, xfin);
    cmesh->InsertMaterialObject(matLambdaBC_left);

    TPZMat1dLin *matLambdaBC_right = new TPZMat1dLin(fmatLambdaBC_right);
    matLambdaBC_right->SetMaterial(xkin, xcin, xbin, xfin);
    cmesh->InsertMaterialObject(matLambdaBC_right);
    
//    TPZMaterial * BCond4 = material->CreateBC(material, fmatInterLambda, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
//    cmesh->InsertMaterialObject(BCond4); //Insere material na malha
    
    
    //    Ponto de pressao:
    //
    TPZFMatrix<STATE> val3(1,1,0.), val4(1,1,0.);
    ////
    TPZMaterial * BCPoint = material->CreateBC(material, fmatPoint, fpointtype, val3, val4); //Cria material que implementa um ponto para a pressao
    cmesh->InsertMaterialObject(BCPoint); //Insere material na malha
    
//    TPZMaterial * BCPoint2 = material->CreateBC(material, fmatPoint, fpointtype, val3, val4); //Cria material que implementa um ponto para a pressao
//    cmesh->InsertMaterialObject(BCPoint2); //Insere material na malha
    
//
    //    //    Ponto de pressao2:
    //    //
    //    TPZFMatrix<STATE> val5(1,1,0.), val6(1,1,0.);
    //    ////
    //    TPZMaterial * BCPoint2 = material->CreateBC(material, matPoint2, pointtype, val5, val6); //Cria material que implementa um ponto para a pressao
    //    cmesh->InsertMaterialObject(BCPoint2); //Insere material na malha
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha
    
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    std::set<int> materialids;
    materialids.insert(fmatID);
    materialids.insert(fpointtype);
    cmesh->AutoBuild(materialids);
    
    gmesh->ResetReference();
  //  cmesh->LoadReferences();
    
    materialids.clear();
    materialids.insert(fmatLambda);
    materialids.insert(fmatLambdaBC_bott);
    materialids.insert(fmatLambdaBC_top);
    materialids.insert(fmatLambdaBC_left);
    materialids.insert(fmatLambdaBC_right);

    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->SetDefaultOrder(pOrder-1);
    cmesh->SetDimModel(1);
    cmesh->AutoBuild(materialids);

    cmesh->LoadReferences();
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    cmesh->AutoBuild();
    
    
    // @omar::
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
}

TPZMultiphysicsCompMesh *MHMBrinkmanTest::CMesh_m(TPZGeoMesh *gmesh, int Space, int pOrder, STATE visco, STATE theta, STATE sigma)
{
    
    //Criando malha computacional:
    int bc_inte_order = 10;
    TPZMultiphysicsCompMesh * cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim); //Insere dimensão do modelo
    cmesh->SetAllCreateFunctionsMultiphysicElem();

    // Criando material:
    
    // 1 - Material volumétrico 2D
    TPZMHMBrinkmanMaterial *material = new TPZMHMBrinkmanMaterial(fmatID,fdim,Space,visco,theta,sigma);

    TPZAutoPointer<TPZFunction<STATE> > fp = new TPZDummyFunction<STATE> (F_source, 5);
    TPZAutoPointer<TPZFunction<STATE> > solp = new TPZDummyFunction<STATE> (Sol_exact, 5);
    ((TPZDummyFunction<STATE>*)fp.operator->())->SetPolynomialOrder(5);
    ((TPZDummyFunction<STATE>*)solp.operator->())->SetPolynomialOrder(5);
//    material->SetForcingFunction(fp); //Caso simples sem termo fonte
    material->SetForcingFunctionExact(solp);
    
    cmesh->InsertMaterialObject(material);
    
    // 1 - Condições de contorno:
    
    TPZFMatrix<STATE> val1(3,3,0.), val2(3,1,0.);
    val2(0,0) = 0.0; // vx -> 0
    val2(1,0) = 0.0; // vy -> 0
    

    val2(1,0) = 1.0;
    TPZBndCond * BC_bott = material->CreateBC(material, fmatBCbott, fneumann, val1, val2);
    BC_bott->SetBCForcingFunction(0, solp);
    cmesh->InsertMaterialObject(BC_bott);
    
    val2(1,0) = 1.0; // vx -> 0
    TPZBndCond * BC_top = material->CreateBC(material, fmatBCtop, fdirichlet, val1, val2);
    BC_top->SetBCForcingFunction(0, solp);
    cmesh->InsertMaterialObject(BC_top);
        
    val2(0,0) = 0.0;
    TPZBndCond * BC_left = material->CreateBC(material, fmatBCleft, fneumann, val1, val2);
    BC_left->SetBCForcingFunction(0, solp);
    cmesh->InsertMaterialObject(BC_left);
        
    val2(0,0) = 0.0;
    TPZBndCond * BC_right = material->CreateBC(material, fmatBCright, fneumann, val1, val2);
    BC_right->SetBCForcingFunction(0, solp);
    cmesh->InsertMaterialObject(BC_right);

    {
//        TPZMaterial * BCond0 = material->CreateBC(material, fmatIntBCbott, fneumann, val1, val2);
//        cmesh->InsertMaterialObject(BCond0);
//
//        TPZMaterial * BCond1 = material->CreateBC(material, fmatIntBCtop, fneumann, val1, val2);
//        cmesh->InsertMaterialObject(BCond1);
//
//        TPZMaterial * BCond2 = material->CreateBC(material, fmatIntBCleft, fneumann, val1, val2);
//        cmesh->InsertMaterialObject(BCond2);
//
//        TPZMaterial * BCond3 = material->CreateBC(material, fmatIntBCright, fneumann, val1, val2);
//        cmesh->InsertMaterialObject(BCond3);
//
//        TPZMaterial * BCond4 = material->CreateBC(material, fmatInterLambda, fneumann, val1, val2);
//        cmesh->InsertMaterialObject(BCond4);
    }
    
    // 2.1 - Material para tração tangencial 1D (Interior)
    TPZMHMBrinkmanMaterial *matLambda = new TPZMHMBrinkmanMaterial(fmatLambda,fdim-1,Space,visco,theta,sigma);
    cmesh->InsertMaterialObject(matLambda);
    
    // 2.2 - Material for interfaces (Interior)
    TPZMHMBrinkmanMaterial *matInterfaceLeft = new TPZMHMBrinkmanMaterial(fmatInterfaceLeft,fdim-1,Space,visco,theta,sigma);
    matInterfaceLeft->SetMultiplier(1.);
    cmesh->InsertMaterialObject(matInterfaceLeft);

    TPZMHMBrinkmanMaterial *matInterfaceRight = new TPZMHMBrinkmanMaterial(fmatInterfaceRight,fdim-1,Space,visco,theta,sigma);
    matInterfaceRight->SetMultiplier(-1.);
    cmesh->InsertMaterialObject(matInterfaceRight);
    
    
    // 3.1 - Material para tração tangencial 1D nos contornos
    TPZMHMBrinkmanBC *matLambdaBC_bott = new TPZMHMBrinkmanBC(fmatLambdaBC_bott,fdim-1,Space,visco,theta,sigma);
    matLambdaBC_bott->SetBC(BC_bott);
    cmesh->InsertMaterialObject(matLambdaBC_bott);

    TPZMHMBrinkmanBC *matLambdaBC_top = new TPZMHMBrinkmanBC(fmatLambdaBC_top,fdim-1,Space,visco,theta,sigma);
    matLambdaBC_top->SetBC(BC_top);
    cmesh->InsertMaterialObject(matLambdaBC_top);

    TPZMHMBrinkmanBC *matLambdaBC_left = new TPZMHMBrinkmanBC(fmatLambdaBC_left,fdim-1,Space,visco,theta,sigma);
    matLambdaBC_left->SetBC(BC_left);
    cmesh->InsertMaterialObject(matLambdaBC_left);
    
    TPZMHMBrinkmanBC *matLambdaBC_right = new TPZMHMBrinkmanBC(fmatLambdaBC_right,fdim-1,Space,visco,theta,sigma);
    matLambdaBC_right->SetBC(BC_right);
    cmesh->InsertMaterialObject(matLambdaBC_right);

    // 3.2 - Interface material para tração tangencial 1D nos contornos
    TPZMHMBrinkmanBC *matInterfaceBC_bott = new TPZMHMBrinkmanBC(fmatIntBCbott,fdim-1,Space,visco,theta,sigma);
    matInterfaceBC_bott->SetBC(BC_bott);
    matInterfaceBC_bott->SetMultiplier(1.);
    cmesh->InsertMaterialObject(matInterfaceBC_bott);

    TPZMHMBrinkmanBC *matInterfaceBC_top = new TPZMHMBrinkmanBC(fmatIntBCtop,fdim-1,Space,visco,theta,sigma);
    matInterfaceBC_top->SetBC(BC_top);
    matInterfaceBC_top->SetMultiplier(1.);
    cmesh->InsertMaterialObject(matInterfaceBC_top);

    TPZMHMBrinkmanBC *matInterfaceBC_left = new TPZMHMBrinkmanBC(fmatIntBCleft,fdim-1,Space,visco,theta,sigma);
    matInterfaceBC_left->SetBC(BC_left);
    matInterfaceBC_left->SetMultiplier(1.);
    cmesh->InsertMaterialObject(matInterfaceBC_left);

    TPZMHMBrinkmanBC *matInterfaceBC_right = new TPZMHMBrinkmanBC(fmatIntBCright,fdim-1,Space,visco,theta,sigma);
    matInterfaceBC_right->SetBC(BC_right);
    matInterfaceBC_right->SetMultiplier(1.);
    cmesh->InsertMaterialObject(matInterfaceBC_right);
    

    //Ponto
    TPZFMatrix<STATE> val3(1,1,0.), val4(1,1,0.);
    val4(0,0)=1.;
    TPZMaterial * BCPoint = material->CreateBC(material, fmatPoint, fpointtype, val3, val4); //Cria material que implementa um ponto para a pressão
    cmesh->InsertMaterialObject(BCPoint); //Insere material na malha
    
    
    
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha:
    
    TPZManVector<int,5> active_approx_spaces(2,1);
    
    cmesh->BuildMultiphysicsSpace(active_approx_spaces,f_mesh_vector);
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    
    
    return cmesh;
    
    
}

void MHMBrinkmanTest::InsertInterfaces(TPZMultiphysicsCompMesh *cmesh_m){
    
    std::set<int> boundaries_ids;
    boundaries_ids.insert(fmatBCbott);
    boundaries_ids.insert(fmatBCleft);
    boundaries_ids.insert(fmatBCtop);
    boundaries_ids.insert(fmatBCright);
    
    TPZInterfaceInsertion InterfaceInsertion(cmesh_m, fmatLambda, boundaries_ids);
    TPZManVector<int64_t,3> Interfaces(2,0);
    Interfaces[0] = fmatInterfaceLeft;
    Interfaces[1] = fmatInterfaceRight;
    InterfaceInsertion.SetInterfaceVectorId(Interfaces);
    
    //InterfaceInsertion.InsertHdivBound(fmatWrap);
    InterfaceInsertion.AddMultiphysicsInterfacesLeftNRight(fmatLambda);
   // InterfaceInsertion.SetMultiplierBCMatId(fmatLambdaBC);
    InterfaceInsertion.AddMultiphysicsBCInterface(fmatLambdaBC_bott,fmatIntBCbott);
    InterfaceInsertion.AddMultiphysicsBCInterface(fmatLambdaBC_top,fmatIntBCtop);
    InterfaceInsertion.AddMultiphysicsBCInterface(fmatLambdaBC_left,fmatIntBCleft);
    InterfaceInsertion.AddMultiphysicsBCInterface(fmatLambdaBC_right,fmatIntBCright);

}
