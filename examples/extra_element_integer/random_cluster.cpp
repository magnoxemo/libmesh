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
  system.solution->close();

}

void AddVariablesToSystem(libMesh::EquationSystems &equation_systems,
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

double GetElementDataFromMesh(LinearImplicitSystem &system,const libMesh::Elem& elem,
                                    const unsigned int variable_num){

    DofMap &dof_map = system.get_dof_map();
    std::vector<dof_id_type> dof_indices;
    dof_map.dof_indices(&elem, dof_indices,variable_num);

    std::vector<double> solution_value(1);
    system.solution->get(dof_indices, solution_value);

    return static_cast<double> (solution_value[0]);
}

const unsigned int FindCluster(libMesh::Mesh &mesh, LinearImplicitSystem &system,
                 const std::string & variable_name) {

  const unsigned int variable_num = system.variable_number(variable_name);

  const unsigned int index = mesh.add_elem_integer(variable_name);

  for (const auto &elem : mesh.element_ptr_range()) {

    int element_solution = static_cast<int>(GetElementDataFromMesh(system,*elem,variable_num));

    bool belong_to_a_cluster = false;
    for (unsigned int side = 0; side < elem->n_sides(); ++side) {
      const Elem *neighbor = elem->neighbor_ptr(side);
      if (neighbor) {

        int neighbor_element_solution = static_cast<int>(GetElementDataFromMesh(system,*neighbor,variable_num));

        if (element_solution == neighbor_element_solution) {
          belong_to_a_cluster = true;
          //Danger Zone for future ref
          //cluster_id should never be equal to element_solution
          unsigned int cluster_id =element_solution;
          elem->set_extra_integer(index,cluster_id);
          break;
        }
      }
    }
    if (!belong_to_a_cluster){
        elem->set_extra_integer(index,0);
    }
  }
  return index;

}

void CaptureClusterID(libMesh::Mesh & mesh, libMesh::EquationSystems &equation_systems,const unsigned int index){

    LinearImplicitSystem &local_system = equation_systems.add_system<LinearImplicitSystem>("cluster_id");
    DofMap &local_dof_map = local_system.get_dof_map();
    AddVariablesToSystem(equation_systems, local_system, "cluster_id_field",true);

    for (const auto & elem:mesh.element_ptr_range()){
        std::vector<dof_id_type> local_dof_indices;
        local_dof_map.dof_indices(elem, local_dof_indices);
        int cluster_id = elem->get_extra_integer(index);
        local_system.solution->set(local_dof_indices[0], cluster_id);
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
  AddVariablesToSystem(equation_systems, system, "random_field");
  PopulateRandomIntegers(mesh, system);

  const unsigned int index =  FindCluster(mesh, system, "random_field");
  CaptureClusterID(mesh,equation_systems,index);


  CloseSystems(mesh, system, equation_systems, false);

  ExodusII_IO(mesh).write_discontinuous_equation_systems("output_rand.e",
                                                         equation_systems);

  return 0;
}