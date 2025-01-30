//
// Created by ebny_walid on 1/30/25.
//
//
// Created by ebny_walid on 1/30/25.
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
    unsigned int nx = 40;
    unsigned int ny = 40;
    srand(time(0));

    libMesh::MeshTools::Generation::build_square(mesh,
                                                 nx,
                                                 ny,
                                                 0.0, 10.0,
                                                 0.0, 10.0,
                                                 libMesh::QUAD4);

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
                const unsigned int cluster_id =elem->get_extra_integer(index);
                system.solution->set(dof_index, cluster_id);
            }
        }
    }

    /* Now that I have some sort of random field here */



    system.solution->close();

    mesh.print_info();
    equation_systems.print_info();

    ExodusII_IO(mesh).write_discontinuous_equation_systems("output_rand.e", equation_systems);

    return 0;
}
