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
                              const libMesh::Elem *elem,
                              const unsigned int variable_num) {

    DofMap &dof_map = system.get_dof_map();
    std::vector<dof_id_type> dof_indices;
    dof_map.dof_indices(elem, dof_indices, variable_num);

    std::vector<double> solution_value(1);
    system.solution->get(dof_indices, solution_value);

    return static_cast<double>(solution_value[0]);
}

bool BelongToCluster(libMesh::Elem *elem, libMesh::Elem *neighbor_elem,
                     int variable_num, libMesh::LinearImplicitSystem &system) {

    double element_solution = GetElementDataFromMesh(system, elem, variable_num);
    double neighbor_solution =
            GetElementDataFromMesh(system, neighbor_elem, variable_num);
    if (static_cast<int>(element_solution) ==
        static_cast<int>(neighbor_solution)) {

        return true;
    } else {
        return false;
    }
}

const unsigned int FindCluster(libMesh::Mesh &mesh,
                               libMesh::LinearImplicitSystem &system,
                               const std::string &variable_name) {
    const unsigned int variable_num = system.variable_number(variable_name);
    int not_visited = -1;
    const unsigned int index = mesh.add_elem_integer(variable_name,not_visited);
    std::stack<libMesh::Elem *> neighbor_stack;

    for (const auto &elem : mesh.element_ptr_range()) {
        neighbor_stack.push(elem);
        int cluster_id = elem->id();

        while (!neighbor_stack.empty()) {
            libMesh::Elem *test_elem = neighbor_stack.top();
            neighbor_stack.pop();
            if (test_elem->get_extra_integer(index) == not_visited &&
                BelongToCluster(elem, test_elem, variable_num, system)) {

                test_elem->set_extra_integer(index, cluster_id);
                for (unsigned int s = 0; s < test_elem->n_sides(); s++) {
                    libMesh::Elem *neighbor_elem = test_elem->neighbor_ptr(s);
                    if (neighbor_elem && neighbor_elem->get_extra_integer(index) == not_visited) {
                        neighbor_stack.push(neighbor_elem);
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

int main(int argc, char **argv) {
    LibMeshInit init(argc, argv);

    Mesh mesh(init.comm());
    unsigned int nx = 5;
    unsigned int ny = 5;
    srand(673);

    CreateMesh(mesh, nx, ny);
    EquationSystems equation_systems(mesh);
    LinearImplicitSystem &system =
            equation_systems.add_system<LinearImplicitSystem>(
                    "random_solution_field");
    AddVariablesToSystem(equation_systems, system, "random_field");
    PopulateRandomIntegers(mesh, system);

    const unsigned int index = FindCluster(mesh, system, "random_field");
    CaptureClusterID(mesh, equation_systems, index);
    ExodusII_IO(mesh).write_discontinuous_equation_systems("output_new.e",
                                                           equation_systems);
    return 0;
}
