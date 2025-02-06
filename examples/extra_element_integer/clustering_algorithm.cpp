#include "clustering_algorithm.h"

using namespace libMesh;

int main(int argc, char **argv) {
    LibMeshInit init(argc, argv);

    Mesh mesh(init.comm());
    unsigned int nx = 4;
    unsigned int ny = 4;
    srand(time(0));

    CreateMesh(mesh, nx, ny);
    EquationSystems equation_systems(mesh);
    LinearImplicitSystem &system =
            equation_systems.add_system<LinearImplicitSystem>(
                    "random_solution_field");
    AddVariablesToSystem(equation_systems, system, "random_field");
    PopulateRandomIntegers(mesh, system);

    const unsigned int index = FindCluster(mesh, system, "random_field");
    CaptureClusterID(mesh, equation_systems, index);
    ExodusII_IO(mesh).write_discontinuous_equation_systems("output_rand.e", equation_systems);

    return 0;
}
