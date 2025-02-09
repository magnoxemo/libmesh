#include <cstdlib>
#include <ctime>
#include <random>
#include <vector>

#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/libmesh.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/system.h"

using namespace libMesh;

void CreateMesh(libMesh::Mesh &mesh, int nx, int ny) {
    libMesh::MeshTools::Generation::build_square(mesh, nx, ny, 0.0, 10.0, 0.0,
                                                 10.0);
}

void CreateMesh(libMesh::Mesh &mesh, int nx, int ny, int nz) {
    libMesh::MeshTools::Generation::build_cube(mesh, nx, ny, nz, 0.0, 10.0, 0.0,
                                               10.0, 0.0, 10.0);
}

void PopulateRandomIntegers(libMesh::Mesh &mesh, LinearImplicitSystem &system) {
    DofMap &dof_map = system.get_dof_map();
    std::vector<dof_id_type> dof_indices;
    for (const auto &elem : mesh.element_ptr_range()) {

        dof_map.dof_indices(elem, dof_indices);
        int random_int = rand() % 5;
        system.solution->set(dof_indices[0], 5 - random_int);
    }
    system.solution->close();
}

void AddVariablesToSystem(libMesh::EquationSystems &equation_systems,
                          LinearImplicitSystem &system, std::string var_name,
                          bool needs_re_init = false) {
    system.add_variable(var_name, CONSTANT, MONOMIAL);
    if (!needs_re_init) {
        equation_systems.init();
    } else {
        equation_systems.reinit();
    }
}

double GetElementDataFromMesh(LinearImplicitSystem &system,
                              const libMesh::Elem &elem,
                              const unsigned int variable_num) {

    DofMap &dof_map = system.get_dof_map();
    std::vector<dof_id_type> dof_indices;
    dof_map.dof_indices(&elem, dof_indices, variable_num);

    std::vector<double> solution_value(1);
    system.solution->get(dof_indices, solution_value);

    return static_cast<double>(solution_value[0]);
}

void SetElementIDasClusterID(libMesh::Mesh &mesh, unsigned int index) {
    for (const auto &elem : mesh.element_ptr_range()) {
        elem->set_extra_integer(index, elem->id());
    }
}

bool BelongToCluster(double elem_solution, double neighbor_solution) {
    if (static_cast<int>(elem_solution) == static_cast<int>(neighbor_solution)) {
        // horrible idea but for now OK
        return true;
    } else {
        return false;
    }
}

void ApplyRecursiveClustering(LinearImplicitSystem &system, unsigned int parent_element_id,
                              libMesh::Elem &current_elem,
                              const unsigned int variable_num,
                              const unsigned int index) {
    // Iterate through all sides of the current element
// rewrite everything
    for (unsigned int side = 0; side < current_elem.n_sides(); side++) {
        // Get the neighbor element
         libMesh::Elem *neighbor_elem = current_elem.neighbor_ptr(side);
        //need to ensure that curent_elem isn't messing with the main one
        if (neighbor_elem && (neighbor_elem->id() != parent_element_id)) {
            // goes on only if the test passes
            double element_solution = GetElementDataFromMesh(system, current_elem, variable_num);
            double neighbor_solution = GetElementDataFromMesh(system, *neighbor_elem, variable_num);
            if (BelongToCluster(element_solution, neighbor_solution)) {
                neighbor_elem->set_extra_integer(index, parent_element_id);
                std::cout<<"     new element to the cluster = "<<neighbor_elem->id()<< " derived from = "<<current_elem.id()<<std::endl;
                ApplyRecursiveClustering(system, current_elem.id(), *neighbor_elem, variable_num, index);
            }
        }
    }

}

const unsigned int FindCluster(libMesh::Mesh &mesh,
                               LinearImplicitSystem &system,
                               const std::string &variable_name) {
    const unsigned int variable_num = system.variable_number(variable_name);
    const unsigned int index = mesh.add_elem_integer(variable_name);

    SetElementIDasClusterID(mesh, index);

    for (const auto &elem : mesh.element_ptr_range()) {

            for (unsigned int side = 0; side < elem->n_sides(); side++) {
                if (elem->get_extra_integer(index) == elem->id()) {
                libMesh::Elem *neighbor_elem = elem->neighbor_ptr(side);
                if (neighbor_elem) {
                    double element_solution = GetElementDataFromMesh(system, *elem, variable_num);
                    double neighbor_solution = GetElementDataFromMesh(system, *neighbor_elem, variable_num);
                    if (BelongToCluster(element_solution, neighbor_solution)) {
                        /*I have to pass the parent element here cause there is a big chance when doing it recursively
                         * the neighbor element would end up setting the extra_integer to the first element
                         * I have to think about a way to prevent it.
                         */
                        neighbor_elem->set_extra_integer(index, elem->id());
                        std::cout<<"main_element_id = "<<elem->id()<<" cluster_elem_id = "<<neighbor_elem->id()<<std::endl;
                        ApplyRecursiveClustering(system, elem->id(), *neighbor_elem, variable_num, index);
                    }
                }
            }
        }
    }

    return index;
}

void CaptureClusterID(libMesh::Mesh &mesh,
                      libMesh::EquationSystems &equation_systems,
                      const unsigned int index) {

    LinearImplicitSystem &local_system =
            equation_systems.add_system<LinearImplicitSystem>("cluster_id");
    AddVariablesToSystem(equation_systems, local_system, "cluster_id_field",
                         true);

    DofMap &local_dof_map = local_system.get_dof_map();

    for (const auto &elem : mesh.element_ptr_range()) {
        std::vector<dof_id_type> local_dof_indices;
        local_dof_map.dof_indices(elem, local_dof_indices);
        int cluster_id = elem->get_extra_integer(index);
        local_system.solution->set(local_dof_indices[0], cluster_id);
    }

    local_system.solution->close();
}

void PrintSystemInformation(libMesh::Mesh &mesh,
                            libMesh::EquationSystems &equation_systems) {
    mesh.print_info();
    equation_systems.print_info();
}