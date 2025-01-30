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

using namespace libMesh;

int main(int argc, char** argv) {

    LibMeshInit init(argc, argv);

    Mesh mesh(init.comm());
    unsigned int nx = 40;  
    unsigned int ny = 40;  
    libMesh::MeshTools::Generation::build_square(mesh,
                                                  nx,
                                                  ny,
                                                  0.0, 10.0,  
                                                  0.0, 10.0, 
                                                  libMesh::QUAD4);

    const unsigned int index = mesh.add_elem_integer("cluster_id");

    for (const auto& elem : mesh.element_ptr_range()) {
        
        int id = elem->id() * 2 + 1;

        if (elem->id() % 2 == 0) {
            elem->set_extra_integer(index, 2);
        } else if (elem->id() % 3 == 0) {
            elem->set_extra_integer(index,  3);
        } else {
            elem->set_extra_integer(index,  4);
        }
    }

    EquationSystems equation_systems(mesh);

    LinearImplicitSystem& system = equation_systems.add_system<LinearImplicitSystem>("cluster");
    system.add_variable("cluster_id", CONSTANT, MONOMIAL); 

    equation_systems.init();

    std::cout << "Number of DoFs: " << system.n_dofs() << std::endl;

    DofMap& dof_map = system.get_dof_map();

    for (const auto& elem : mesh.element_ptr_range()) {
        std::cout << "Element " << elem->id() << " type: " << elem->type() << std::endl;

        int cluster_id = elem->get_extra_integer(index);

        std::vector<dof_id_type> dof_indices;
        dof_map.dof_indices(elem, dof_indices);

        std::cout << "Element " << elem->id() << " DoF indices: ";
        for (auto dof_index : dof_indices) {
            std::cout << dof_index << " ";
        }
        std::cout << std::endl;

        for (auto dof_index : dof_indices) {
            system.solution->set(dof_index, cluster_id);
        }
    }

    system.solution->close();

    mesh.print_info();
    equation_systems.print_info();

    ExodusII_IO(mesh).write_equation_systems("output.exo", equation_systems);

    return 0;
}
