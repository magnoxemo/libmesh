#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/system.h"
#include "libmesh/mesh_generation.h"

#include <vector>
#include <cstdlib>
#include <ctime>

using namespace libMesh;

void CreateMesh(libMesh::Mesh & mesh, int nx, int ny){

    libMesh::MeshTools::Generation::build_square(mesh,
                                                 nx,
                                                 ny,
                                                 0.0, 10.0,
                                                 0.0, 10.0,
                                                 libMesh::QUAD4);
}

void CreateMesh(libMesh::Mesh & mesh, int nx, int ny, int nz ){
    libMesh::MeshTools::Generation::build_cube(mesh,
                                               nx,
                                               ny,
                                               nz,
                                               0.0, 10.0,
                                               0.0, 10.0,
                                               0.0, 10.0);
}

void PopulateRandomIntegers(libMesh::Mesh & mesh, LinearImplicitSystem& system, const unsigned int &index ) {

    DofMap& dof_map = system.get_dof_map();
    for (const auto& elem : mesh.element_ptr_range()) {

        std::vector<dof_id_type> dof_indices;
        dof_map.dof_indices(elem, dof_indices);
        int random_int = rand() %5;

        system.solution->set(dof_indices[1], random_int);
        elem->set_extra_integer(index, 5-random_int);

    }
}

void AddVaribalesToSystem(LinearImplicitSystem& system, libMesh::EquationSystems & equation_systems){

    system.add_variable("cluster_id", CONSTANT, MONOMIAL);
    system.add_variable("random_field", CONSTANT, MONOMIAL);

    equation_systems.init();
}

DofMap& CreateDegreesofFreedomMap(LinearImplicitSystem& system){
    return  system.get_dof_map();
}

void FindCluster(libMesh::Mesh & mesh, libMesh::DofMap & dof_map, LinearImplicitSystem& system, const unsigned int index ) {

    for (const auto &elem: mesh.element_ptr_range()) {

        std::vector <dof_id_type> dof_indices;
        dof_map.dof_indices(elem, dof_indices);

        bool belong_to_a_cluster = false;

        for (unsigned int side = 0; side < elem->n_sides(); ++side) {
            const Elem *neighbor = elem->neighbor_ptr(side);

            if (neighbor) {
                if (elem->get_extra_integer(index) == neighbor->get_extra_integer(index)) {
                    belong_to_a_cluster = true;
                }
            }
        }
        if (!belong_to_a_cluster) {
            system.solution->set(dof_indices[0], 0);
        } else if (belong_to_a_cluster) {
            const unsigned int cluster_id = elem->get_extra_integer(index);
            system.solution->set(dof_indices[0], cluster_id);
        }
    }
}
void CloseSystems(libMesh::Mesh & mesh, LinearImplicitSystem& system,libMesh::EquationSystems & equation_systems, bool print_system_info = true){

    system.solution->close();
    if (print_system_info){
        mesh.print_info();
        equation_systems.print_info();
    }
}

void WriteEquationSystemOnTheMesh(libMesh::Mesh & mesh, libMesh::EquationSystems & equation_systems ){

    ExodusII_IO(mesh).write_discontinuous_equation_systems("output_rand.e", equation_systems);
}

int main(int argc, char** argv) {

    LibMeshInit init(argc, argv);

    Mesh mesh(init.comm());
    unsigned int nx = 40;
    unsigned int ny = 40;
    srand(time(0));

    CreateMesh(mesh,nx,ny);
    const unsigned int index = mesh.add_elem_integer("cluster_id");

    EquationSystems equation_systems(mesh);
    LinearImplicitSystem& system = equation_systems.add_system<LinearImplicitSystem>("cluster");
    DofMap& dof_map = CreateDegreesofFreedomMap(system);
    AddVaribalesToSystem (system, equation_systems);
    PopulateRandomIntegers (mesh,system,index);
    FindCluster(mesh,dof_map,system,index);
    CloseSystems(mesh,system,equation_systems);

    WriteEquationSystemOnTheMesh(mesh,equation_systems);

    return 0;
}