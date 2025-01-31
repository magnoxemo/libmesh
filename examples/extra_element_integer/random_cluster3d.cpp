//
// Created by ebny_walid on 1/31/25.
//
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

int main(int argc, char** argv) {

    LibMeshInit init(argc, argv);

    Mesh mesh(init.comm());
    unsigned int nx = 100;
    unsigned int ny = 100;
    unsigned int nz = 100;
    srand(time(0));

    libMesh::MeshTools::Generation::build_cube(mesh,
                                               nx,
                                               ny,
                                               nz,
                                               0.0, 10.0,
                                               0.0, 10.0,
                                               0.0, 10.0);

    const unsigned int index = mesh.add_elem_integer("cluster_id");

    EquationSystems equation_systems(mesh);

    LinearImplicitSystem& system = equation_systems.add_system<LinearImplicitSystem>("cluster");
    system.add_variable("cluster_id", CONSTANT, MONOMIAL);
    system.add_variable("random_field", CONSTANT, MONOMIAL);

    equation_systems.init();

    std::cout << "Number of DoFs: " << system.n_dofs() << std::endl;

    DofMap& dof_map = system.get_dof_map();

    for (const auto& elem : mesh.element_ptr_range()) {

        std::vector<dof_id_type> dof_indices;
        dof_map.dof_indices(elem, dof_indices);


        /*I am gonna create a random number here from 1 to 4*/
        int random_int = rand() %5;

        for (auto dof_index : dof_indices) {
            if (dof_index %2 !=0){
                system.solution->set(dof_index, random_int);

                if (random_int ==4){
                    elem->set_extra_integer(index, 1);
                } else if (random_int ==3){
                    elem->set_extra_integer(index, 2);
                }else if (random_int ==2){
                    elem->set_extra_integer(index, 3);
                }else if (random_int ==1){
                    elem->set_extra_integer(index, 4);
                }else if (random_int ==0){
                    elem->set_extra_integer(index, 5);
                }
            }
        }
    }

    for (const auto& elem : mesh.element_ptr_range()) {

        std::vector<dof_id_type> dof_indices;
        dof_map.dof_indices(elem, dof_indices);

        for (auto dof_index : dof_indices) {
            if (dof_index %2 ==0){

                bool belong_to_a_cluster =false ;
                for (unsigned int side = 0; side < elem->n_sides(); ++side){
                    const Elem* neighbor = elem->neighbor_ptr(side);

                    if (neighbor){
                        if (elem->get_extra_integer(index)==neighbor->get_extra_integer(index)){
                            belong_to_a_cluster = true ;
                        }
                    }
                }
                if (!belong_to_a_cluster){
                    system.solution->set(dof_index,0);
                }else if (belong_to_a_cluster){
                    const unsigned int cluster_id =elem->get_extra_integer(index);
                    system.solution->set(dof_index, cluster_id);
                }

            }
        }
    }

    system.solution->close();

    mesh.print_info();
    equation_systems.print_info();

    ExodusII_IO(mesh).write_discontinuous_equation_systems("random_cluster_3D.e", equation_systems);

    return 0;
}
