//
//  TPZInterfaceInsertion.h
//  Benchmark0a
//  Class that stores Interface neighbor data in terms of geometric element indexes
//  Created by Pablo Carvalho on 02/08/18.
//

#ifndef TPZInterfaceInsertion_h
#define TPZInterfaceInsertion_h

#include <stdio.h>
#include <iostream>
#include <set>
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "TPZMultiphysicsCompMesh.h"

class TPZInterfaceInsertion {
    
private:

    /// Mesh geometry
    TPZCompMesh *m_cmesh;
    
    TPZGeoMesh * m_geometry;
        
    /// Wrap flux id
    int m_id_flux_wrap;
    
    /// Interface id
    int m_interface_id;

    /// Multiplier material id
    int m_multiplier_id;
    
    /// Set of boundary material ids
    std::set<int> m_boundaries_ids;

public:
    
    /// @TODO:: OD, Rename TPZInterfaceInsertion -> TPZInterfaceDescription
    /// @TODO:: OD, Refactor and rename the methods dependent on the approximation space
    
    /// Default constructor
    TPZInterfaceInsertion();
    
    /// Default desconstructor
    ~TPZInterfaceInsertion();
    
    /// Copy constructor
    TPZInterfaceInsertion(TPZInterfaceInsertion & other);
    
    /// Constructor based on a computational mesh and Interface material id
    TPZInterfaceInsertion(TPZCompMesh *m_cmesh, int Interface_id, std::set<int> & boundaries_ids);
    
    /// Set Interface Identifier
    void SetInterfaceIdentifier(int Interface_id);
    
    /// Get Interface material Identifier
    int & GetInterfaceMaterialId();

    /// Set wrap Identifier
    void SetWrapFluxIdentifier(int wrapFlux);
    
    /// Get wrap Identifier
    int & GetWrapFluxId();

    /// Set multiplier material id
    void SetMultiplierMatId(int wrapFlux);
    
    /// Get multiplier material id
    int & GetMultiplierMatId();
    
    /// Add multiphysics interfaces
    void AddMultiphysicsInterfaces();
    
    void AddMultiphysicsInterfaces(int matfrom, int mattarget);
    
    void AddMultiphysicsInterfacesLeftNRight(int matfrom);
    
    /// Open the connects of a Interface, create dim-1 Interface elements (Hdiv version)
    void InsertHdivBound(int mat_id_flux_wrap);
    
};


#endif /* TPZInterfaceInsertion_h */

