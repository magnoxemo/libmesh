//
// Created by ebny_walid on 1/26/25.
//

/*
 * this example reads a mesh from a path
 * and them declares an extra_interger field called cluster_id
 * iterates through all the elements and sets extra_element_int to the element
 * again iterates through the elements and gets that extra element integer that we previously assigned
 * lastly prints out the mesh info
 */

#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/elem.h"

using namespace libMesh;

int main(int argc, char** argv) {

    LibMeshInit init(argc, argv);

    Mesh mesh(init.comm());
    libMesh::ExodusII_IO exodus_io(mesh);
    exodus_io.read("mesh_in.e");
    const unsigned int index =0;
    for (const auto& elem : mesh.element_ptr_range()) {
        std::cout << "Element ID: " << elem->id()<<" Elem type: "<<elem->type()<< " cluster_id: " << elem->get_extra_integer(index) << std::endl;
    }
    std::cout<<mesh.n_elem_integers()<<"\n\n";

    mesh.print_info();
    mesh.write("mesh.e");

    return 0;
}
