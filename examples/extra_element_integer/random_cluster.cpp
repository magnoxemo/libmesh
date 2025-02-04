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

#include <cstdlib>
#include <ctime>
#include <random>
#include <vector>

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
  for (const auto &elem : mesh.element_ptr_range()) {

    std::vector<dof_id_type> dof_indices;
    dof_map.dof_indices(elem, dof_indices);
    int random_int = rand() % 5;
    system.solution->set(dof_indices[0], 5 - random_int);
  }
}

void AddVaribalesToSystem(libMesh::EquationSystems &equation_systems,
                          LinearImplicitSystem &system,
                          std::string var_name, bool needs_re_init =false) {
    system.add_variable(var_name, CONSTANT, MONOMIAL);
    if (!needs_re_init){
        equation_systems.init();
    }
    else{
        equation_systems.reinit();
    }
}

void FindCluster(libMesh::Mesh &mesh, libMesh::DofMap &dof_map,
                 libMesh::DofMap &local_dof_map, LinearImplicitSystem &system,
                 LinearImplicitSystem &local_system,
                 const std::string &variable_name,
                 libMesh::EquationSystems &equation_systems,
                 const unsigned int index) {

  const unsigned int variable_num = system.variable_number(variable_name);

  for (const auto &elem : mesh.element_ptr_range()) {

    std::vector<dof_id_type> global_dof_indices;
    dof_map.dof_indices(elem, global_dof_indices, variable_num);

    std::vector<double> solution_value(1);
    system.solution->get(global_dof_indices, solution_value);
    int element_solution = static_cast<int>(solution_value[0]);

    bool belong_to_a_cluster = false;
    for (unsigned int side = 0; side < elem->n_sides(); ++side) {
      const Elem *neighbor = elem->neighbor_ptr(side);
      if (neighbor) {
        std::vector<dof_id_type> neighbor_dof_indices;
        dof_map.dof_indices(neighbor, neighbor_dof_indices, variable_num);

        std::vector<double> neighbor_solution_value(1);
        system.solution->get(neighbor_dof_indices, neighbor_solution_value);
        int neighbor_element_solution =
            static_cast<int>(neighbor_solution_value[0]);

        if (element_solution == neighbor_element_solution) {
          belong_to_a_cluster = true;
          break;
        }
      }
    }

    std::vector<dof_id_type> local_dof_indices;
    local_dof_map.dof_indices(elem, local_dof_indices);
    if (!belong_to_a_cluster) {
      elem->set_extra_integer(index, 0);
      local_system.solution->set(local_dof_indices[0], 0);
    } else {
      /*To do:
       * maybe figure out another way to find out what should be the value of
       * the extra_integer rather than just element solution
       * */
      elem->set_extra_integer(index, static_cast<int>(element_solution));
      local_system.solution->set(local_dof_indices[0], element_solution);
    }
  }

  local_system.solution->close();
}

void CloseSystems(libMesh::Mesh &mesh, LinearImplicitSystem &system,
                  libMesh::EquationSystems &equation_systems,
                  bool print_system_info = true) {
  system.solution->close();
  if (print_system_info) {
    mesh.print_info();
    equation_systems.print_info();
  }
}

int main(int argc, char **argv) {
  LibMeshInit init(argc, argv);

  Mesh mesh(init.comm());
  unsigned int nx = 40;
  unsigned int ny = 40;
  srand(time(0));

  CreateMesh(mesh, nx, ny);

  EquationSystems equation_systems(mesh);

  //creating pesudo solution field
  LinearImplicitSystem &system = equation_systems.add_system<LinearImplicitSystem>("random_solution_field");
  DofMap &dof_map = system.get_dof_map();
  AddVaribalesToSystem(equation_systems, system, "random_field");
  PopulateRandomIntegers(mesh, system);

  //creating cluster_id_field
  LinearImplicitSystem &local_system = equation_systems.add_system<LinearImplicitSystem>("cluster_id");
  DofMap &local_dof_map = local_system.get_dof_map();
  AddVaribalesToSystem(equation_systems, local_system, "cluster_id_field",true);
  const unsigned int index = mesh.add_elem_integer("cluster_id");
  FindCluster(mesh, dof_map, local_dof_map, system, local_system, "random_field", equation_systems, index);

  CloseSystems(mesh, system, equation_systems, false);

  ExodusII_IO(mesh).write_discontinuous_equation_systems("output_rand.e",
                                                         equation_systems);

  return 0;
}