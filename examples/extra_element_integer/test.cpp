/*
 * this example reads a mesh from a path
 * and them declares an extra_interger field called cluster_id
 * iterates through all the elements and prints out the centroid of  the elements
 * lastly prints out the mesh info
 */
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_generation.h"

using namespace libMesh;

int main(int argc, char** argv) {
    LibMeshInit init(argc, argv);

    Mesh mesh(init.comm());

    MeshTools::Generation::build_square (mesh, 50, 50);

    const unsigned int cluster_index = mesh.add_elem_integer("cluster_id");

    for (auto& elem : mesh.element_ptr_range()) {
        elem->set_extra_integer(cluster_index, elem->id() * 2 + 1);
        std::cout<<cluster_index<<" "<<elem->get_extra_integer(cluster_index)<<std::endl;
    }
    mesh.write("mesh_out.e");

    return 0;
}