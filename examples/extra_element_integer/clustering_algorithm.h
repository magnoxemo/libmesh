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

void SetElementIDasClusterID(libMesh::Mesh &mesh,unsigned int index){
    for (const auto &elem:mesh.element_ptr_range()){
        elem->set_extra_integer(index,elem->id());
    }
}

bool BelongToCluster(double elem_solution,double neighbor_solution){
    if (static_cast<int>(elem_solution)==static_cast<int>(neighbor_solution)){
        //horible idea but for now OK
        return true;
    }else{
        return false;
    }
}

void ApplyRecursiveClustering(LinearImplicitSystem& system, const libMesh::Elem& parent_elem, libMesh::Elem& elem, const unsigned int &variable_num, const unsigned int &index) {
    const libMesh::Elem& parent = parent_elem;
    if (elem.get_extra_integer(index) == elem.id()) /*failing here */{

        for (unsigned int side = 0; side < elem.n_sides(); side++) {
            const libMesh::Elem* neighbor = elem.neighbor_ptr(side);

            if (neighbor) {
                double neighbor_solution = GetElementDataFromMesh(system, *neighbor, variable_num);
                double elem_solution = GetElementDataFromMesh(system,elem, variable_num);

                if (BelongToCluster(elem_solution,neighbor_solution)) {
                    //const_cast<libMesh::Elem*>(neighbor)->set_extra_integer(index, elem.get_extra_integer(index));
                    const_cast<libMesh::Elem*>(neighbor)->set_extra_integer(index, 0); //another horible idea but I am gonna blame paraview here
                    ApplyRecursiveClustering(system, parent, const_cast<libMesh::Elem&>(*neighbor), variable_num, index);
                }
            }
        }
    }
}

const unsigned int FindCluster(libMesh::Mesh &mesh, LinearImplicitSystem &system, const std::string &variable_name) {
    const unsigned int variable_num = system.variable_number(variable_name);
    const unsigned int index = mesh.add_elem_integer(variable_name);

    SetElementIDasClusterID(mesh, index);

    for (const auto &elem : mesh.element_ptr_range()) {
        if (elem->get_extra_integer(index) == elem->id()) { // If yes, then they aren't part of any cluster

            for (unsigned int side = 0; side < elem->n_sides(); side++) {
                const libMesh::Elem* neighbor = elem->neighbor_ptr(side);
                if (neighbor) {
                    ApplyRecursiveClustering(system, *elem, const_cast<libMesh::Elem&>(*neighbor), variable_num, index);
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
