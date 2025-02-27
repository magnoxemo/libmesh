// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh_logging.h"


// Local Includes
#include "libmesh/libmesh_version.h"
#include "libmesh/equation_systems.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/parallel.h"
#include "libmesh/xdr_cxx.h"
#include "libmesh/mesh_refinement.h"

// C++ Includes
#include <iomanip> // setfill
#include <sstream>
#include <string>

namespace libMesh
{

// Forward Declarations

// Anonymous namespace for implementation details.
namespace {
std::string local_file_name (const unsigned int processor_id,
                             std::string_view basename)
{
  std::ostringstream returnval;

  std::string_view suffix;
  if (basename.size() - basename.rfind(".bz2") == 4)
    {
      basename.remove_suffix(4);
      suffix = ".bz2";
    }
  else if (basename.size() - basename.rfind(".gz") == 3)
    {
      basename.remove_suffix(3);
      suffix = ".gz";
    }

  returnval << basename << '.';
  returnval << std::setfill('0') << std::setw(4);
  returnval << processor_id;
  returnval << suffix;

  return returnval.str();
}
}




// ------------------------------------------------------------
// EquationSystem class implementation
template <typename InValType>
void EquationSystems::read (std::string_view name,
                            const unsigned int read_flags,
                            bool partition_agnostic)
{
  XdrMODE mode = READ;
  if (name.find(".xdr") != std::string::npos)
    mode = DECODE;
  this->read(name, mode, read_flags, partition_agnostic);
}



template <typename InValType>
void EquationSystems::read (std::string_view name,
                            const XdrMODE mode,
                            const unsigned int read_flags,
                            bool partition_agnostic)
{
  // This will unzip a file with .bz2 as the extension, otherwise it
  // simply returns the name if the file need not be unzipped.
  Xdr io ((this->processor_id() == 0) ? std::string(name) : "", mode);

  std::function<std::unique_ptr<Xdr>()> local_io_functor;
  local_io_functor = [this,&name,&mode]() {
    return std::make_unique<Xdr>(local_file_name(this->processor_id(), name), mode); };

  this->read(io, local_io_functor, read_flags, partition_agnostic);
}



template <typename InValType>
void EquationSystems::read (Xdr & io,
                            std::function<std::unique_ptr<Xdr>()> & local_io_functor,
                            const unsigned int read_flags,
                            bool partition_agnostic)
{
  /**
   * This program implements the output of an
   * EquationSystems object.  This warrants some
   * documentation.  The output file essentially
   * consists of 11 sections:
   \verbatim
   1.) A version header (for non-'legacy' formats, libMesh-0.7.0 and greater).
   2.) The number of individual equation systems (unsigned int)

   for each system

   3.)  The name of the system (string)
   4.)  The type of the system (string)

   handled through System::read():

   +-------------------------------------------------------------+
   |  5.) The number of variables in the system (unsigned int)   |
   |                                                             |
   |   for each variable in the system                           |
   |                                                             |
   |    6.) The name of the variable (string)                    |
   |                                                             |
   |    7.) Combined in an FEType:                               |
   |         - The approximation order(s) of the variable (Order |
   |           Enum, cast to int/s)                              |
   |         - The finite element family/ies of the variable     |
   |           (FEFamily Enum, cast to int/s)                    |
   |                                                             |
   |   end variable loop                                         |
   |                                                             |
   | 8.) The number of additional vectors (unsigned int),        |
   |                                                             |
   |    for each additional vector in the equation system object |
   |                                                             |
   |    9.) the name of the additional vector  (string)          |
   +-------------------------------------------------------------+

   end system loop


   for each system, handled through System::read_{serialized,parallel}_data():

   +--------------------------------------------------------------+
   | 10.) The global solution vector, re-ordered to be node-major |
   |     (More on this later.)                                    |
   |                                                              |
   |    for each additional vector in the equation system object  |
   |                                                              |
   |    11.) The global additional vector, re-ordered to be       |
   |         node-major (More on this later.)                     |
   +--------------------------------------------------------------+

   end system loop
   \endverbatim
   *
   * Note that the actual IO is handled through the Xdr class
   * (to be renamed later?) which provides a uniform interface to
   * both the XDR (eXternal Data Representation) interface and standard
   * ASCII output.  Thus this one section of code will read XDR or ASCII
   * files with no changes.
   */

  // Set booleans from the read_flags argument
  const bool read_header          = read_flags & EquationSystems::READ_HEADER;
  const bool read_data            = read_flags & EquationSystems::READ_DATA;
  const bool read_additional_data = read_flags & EquationSystems::READ_ADDITIONAL_DATA;
  const bool read_legacy_format   = read_flags & EquationSystems::READ_LEGACY_FORMAT;
  const bool try_read_ifems       = read_flags & EquationSystems::TRY_READ_IFEMS;
  const bool read_basic_only      = read_flags & EquationSystems::READ_BASIC_ONLY;
  bool read_parallel_files  = false;

  std::vector<std::pair<std::string, System *>> xda_systems;

  libmesh_assert (io.reading());

  {
    // 1.)
    // Read the version header.
    std::string version = "legacy";
    if (!read_legacy_format)
      {
        if (this->processor_id() == 0) io.data(version);
        this->comm().broadcast(version);

        // All processors have the version header, if it does not contain
        // the libMesh_label string then it is a legacy file.
        const std::string libMesh_label = "libMesh-";
        std::string::size_type lm_pos = version.find(libMesh_label);
        if (lm_pos==std::string::npos)
          {
            io.close();

            // Recursively call this read() function but with the
            // EquationSystems::READ_LEGACY_FORMAT bit set.
            this->read (io, local_io_functor, (read_flags | EquationSystems::READ_LEGACY_FORMAT), partition_agnostic);
            return;
          }

        // Figure out the libMesh version that created this file
        std::istringstream iss(version.substr(lm_pos + libMesh_label.size()));
        int ver_major = 0, ver_minor = 0, ver_patch = 0;
        char dot;
        iss >> ver_major >> dot >> ver_minor >> dot >> ver_patch;
        io.set_version(LIBMESH_VERSION_ID(ver_major, ver_minor, ver_patch));


        read_parallel_files = (version.rfind(" parallel") < version.size());

        // If requested that we try to read infinite element information,
        // and the string " with infinite elements" is not in the version,
        // then tack it on.  This is for compatibility reading ifem
        // files written prior to 11/10/2008 - BSK
        if (try_read_ifems)
          if (!(version.rfind(" with infinite elements") < version.size()))
            version += " with infinite elements";

      }
    else
      libmesh_deprecated();

    LOG_SCOPE("read()", "EquationSystems");

    // 2.)
    // Read the number of equation systems
    unsigned int n_sys=0;
    if (this->processor_id() == 0) io.data (n_sys);
    this->comm().broadcast(n_sys);

    for (unsigned int sys=0; sys<n_sys; sys++)
      {
        // 3.)
        // Read the name of the sys-th equation system
        std::string sys_name;
        if (this->processor_id() == 0) io.data (sys_name);
        this->comm().broadcast(sys_name);

        // 4.)
        // Read the type of the sys-th equation system
        std::string sys_type;
        if (this->processor_id() == 0) io.data (sys_type);
        this->comm().broadcast(sys_type);

        if (read_header)
          this->add_system (sys_type, sys_name);

        // 5.) - 9.)
        // Let System::read_header() do the job
        System & new_system = this->get_system(sys_name);
        new_system.read_header (io,
                                version,
                                read_header,
                                read_additional_data,
                                read_legacy_format);

        xda_systems.emplace_back(sys_name, &new_system);

        // If we're only creating "basic" systems, we need to tell
        // each system that before we call init() later.
        if (read_basic_only)
          new_system.set_basic_system_only();
      }
  }



  // Now we are ready to initialize the underlying data
  // structures. This will initialize the vectors for
  // storage, the dof_map, etc...
  if (read_header)
    this->init();

  // 10.) & 11.)
  // Read and set the numeric vector values
  if (read_data)
    {
      std::unique_ptr<Xdr> local_io;

      // the EquationSystems::read() method should look constant from the mesh
      // perspective, but we need to assign a temporary numbering to the nodes
      // and elements in the mesh, which requires that we abuse const_cast
      if (!read_legacy_format && partition_agnostic)
        {
          MeshBase & mesh = const_cast<MeshBase &>(this->get_mesh());
          MeshTools::Private::globally_renumber_nodes_and_elements(mesh);
        }

      for (auto & pr : xda_systems)
        if (read_legacy_format)
          {
            libmesh_deprecated();
#ifdef LIBMESH_ENABLE_DEPRECATED
            pr.second->read_legacy_data (io, read_additional_data);
#endif
          }
        else
          if (read_parallel_files)
            {
              if (!local_io)
              {
                local_io = local_io_functor();
                libmesh_assert(local_io->reading());
              }
              pr.second->read_parallel_data<InValType> (*local_io, read_additional_data);
            }
          else
            pr.second->read_serialized_data<InValType> (io, read_additional_data);


      // Undo the temporary numbering.
      if (!read_legacy_format && partition_agnostic)
        _mesh.fix_broken_node_and_element_numbering();
    }

  // Localize each system's data
  this->update();

  #ifdef LIBMESH_ENABLE_AMR
    MeshRefinement mesh_refine(_mesh);
    mesh_refine.clean_refinement_flags();
  #endif
}



void EquationSystems::write(std::string_view name,
                            const unsigned int write_flags,
                            bool partition_agnostic) const
{
  XdrMODE mode = WRITE;
  if (name.find(".xdr") != std::string::npos)
    mode = ENCODE;
  this->write(name, mode, write_flags, partition_agnostic);
}



void EquationSystems::write(std::string_view name,
                            const XdrMODE mode,
                            const unsigned int write_flags,
                            bool partition_agnostic) const
{
  Xdr io((this->processor_id()==0) ? std::string(name) : "", mode);

  std::unique_ptr<Xdr> local_io;
  // open a parallel buffer if warranted
  if (write_flags & EquationSystems::WRITE_PARALLEL_FILES && write_flags & EquationSystems::WRITE_DATA)
    local_io = std::make_unique<Xdr>(local_file_name(this->processor_id(),name), mode);

  this->write(io, write_flags, partition_agnostic, local_io.get());
}



void EquationSystems::write(Xdr & io,
                            const unsigned int write_flags,
                            bool partition_agnostic,
                            Xdr * const local_io) const
{
  /**
   * This program implements the output of an
   * EquationSystems object.  This warrants some
   * documentation.  The output file essentially
   * consists of 11 sections:
   \verbatim
   1.) The version header.
   2.) The number of individual equation systems (unsigned int)

   for each system

   3.)  The name of the system (string)
   4.)  The type of the system (string)

   handled through System::read():

   +-------------------------------------------------------------+
   |  5.) The number of variables in the system (unsigned int)   |
   |                                                             |
   |   for each variable in the system                           |
   |                                                             |
   |    6.) The name of the variable (string)                    |
   |                                                             |
   |    7.) Combined in an FEType:                               |
   |         - The approximation order(s) of the variable (Order |
   |           Enum, cast to int/s)                              |
   |         - The finite element family/ies of the variable     |
   |           (FEFamily Enum, cast to int/s)                    |
   |                                                             |
   |   end variable loop                                         |
   |                                                             |
   | 8.) The number of additional vectors (unsigned int),        |
   |                                                             |
   |    for each additional vector in the equation system object |
   |                                                             |
   |    9.) the name of the additional vector  (string)          |
   +-------------------------------------------------------------+

   end system loop


   for each system, handled through System::write_{serialized,parallel}_data():

   +--------------------------------------------------------------+
   | 10.) The global solution vector, re-ordered to be node-major |
   |     (More on this later.)                                    |
   |                                                              |
   |    for each additional vector in the equation system object  |
   |                                                              |
   |    11.) The global additional vector, re-ordered to be       |
   |         node-major (More on this later.)                     |
   +--------------------------------------------------------------+

   end system loop
   \endverbatim
   *
   * Note that the actual IO is handled through the Xdr class
   * (to be renamed later?) which provides a uniform interface to
   * both the XDR (eXternal Data Representation) interface and standard
   * ASCII output.  Thus this one section of code will write XDR or ASCII
   * files with no changes.
   */

  // the EquationSystems::write() method should look constant,
  // but we need to assign a temporary numbering to the nodes
  // and elements in the mesh, which requires that we abuse const_cast
  if (partition_agnostic)
    {
      MeshBase & mesh = const_cast<MeshBase &>(this->get_mesh());
      MeshTools::Private::globally_renumber_nodes_and_elements(mesh);
    }

  // set booleans from write_flags argument
  const bool write_data            = write_flags & EquationSystems::WRITE_DATA;
  const bool write_additional_data = write_flags & EquationSystems::WRITE_ADDITIONAL_DATA;

  // always write parallel files if we're instructed to write in
  // parallel
  const bool write_parallel_files  =
    (write_flags & EquationSystems::WRITE_PARALLEL_FILES)
    // Even if we're on a distributed mesh, we may or may not have a
    // consistent way of reconstructing the same mesh partitioning
    // later, but we need the same mesh partitioning if we want to
    // reread the parallel solution safely, so let's write a serial file
    // unless specifically requested not to.
    // ||
    // // but also write parallel files if we haven't been instructed to
    // // write in serial and we're on a distributed mesh
    // (!(write_flags & EquationSystems::WRITE_SERIAL_FILES) &&
    // !this->get_mesh().is_serial())
    ;

  if (write_parallel_files && write_data)
    libmesh_assert(local_io);

  {
    libmesh_assert (io.writing());

    LOG_SCOPE("write()", "EquationSystems");

    const unsigned int proc_id = this->processor_id();

    unsigned int n_sys = 0;
    for (auto & pr : _systems)
      if (!pr.second->hide_output())
        n_sys++;

    // set the version number in the Xdr object
    io.set_version(LIBMESH_VERSION_ID(LIBMESH_MAJOR_VERSION,
                                      LIBMESH_MINOR_VERSION,
                                      LIBMESH_MICRO_VERSION));

    // Only write the header information
    // if we are processor 0.
    if (proc_id == 0)
      {
        std::string comment;

        // 1.)
        // Write the version header
        std::string version("libMesh-" + libMesh::get_io_compatibility_version());
        if (write_parallel_files) version += " parallel";

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
        version += " with infinite elements";
#endif
        io.data (version, "# File Format Identifier");

        // 2.)
        // Write the number of equation systems
        io.data (n_sys, "# No. of Equation Systems");

        for (auto & [sys_name, sys] : _systems)
          {
            // Ignore this system if it has been marked as hidden
            if (sys->hide_output()) continue;

            // 3.)
            // Write the name of the sys_num-th system
            {
              const unsigned int sys_num = sys->number();

              comment =  "# Name, System No. ";
              comment += std::to_string(sys_num);

              // Note: There is no Xdr::data overload taking a "const
              // std::string &" so we need to make a copy.
              std::string copy = sys_name;
              io.data (copy, comment);
            }

            // 4.)
            // Write the type of system handled
            {
              const unsigned int sys_num = sys->number();
              std::string sys_type       = sys->system_type();

              comment =  "# Type, System No. ";
              comment += std::to_string(sys_num);

              io.data (sys_type, comment);
            }

            // 5.) - 9.)
            // Let System::write_header() do the job
            sys->write_header (io, version, write_additional_data);
          }
      }

    // Start from the first system, again,
    // to write vectors to disk, if wanted
    if (write_data)
      {
        for (auto & pr : _systems)
          {
            // Ignore this system if it has been marked as hidden
            if (pr.second->hide_output()) continue;

            // 10.) + 11.)
            if (write_parallel_files)
              pr.second->write_parallel_data (*local_io,write_additional_data);
            else
              pr.second->write_serialized_data (io,write_additional_data);
          }

        if (local_io)
          local_io->close();
      }

    io.close();
  }

  // the EquationSystems::write() method should look constant,
  // but we need to undo the temporary numbering of the nodes
  // and elements in the mesh, which requires that we abuse const_cast
  if (partition_agnostic)
    const_cast<MeshBase &>(_mesh).fix_broken_node_and_element_numbering();
}



// template specialization

template LIBMESH_EXPORT void EquationSystems::read<Number> (Xdr & io, std::function<std::unique_ptr<Xdr>()> & local_io_functor, const unsigned int read_flags, bool partition_agnostic);
template LIBMESH_EXPORT void EquationSystems::read<Number> (std::string_view name, const unsigned int read_flags, bool partition_agnostic);
template LIBMESH_EXPORT void EquationSystems::read<Number> (std::string_view name, const XdrMODE mode, const unsigned int read_flags, bool partition_agnostic);
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
template LIBMESH_EXPORT void EquationSystems::read<Real> (Xdr & io, std::function<std::unique_ptr<Xdr>()> & local_io_functor, const unsigned int read_flags, bool partition_agnostic);
template LIBMESH_EXPORT void EquationSystems::read<Real> (std::string_view name, const unsigned int read_flags, bool partition_agnostic);
template LIBMESH_EXPORT void EquationSystems::read<Real> (std::string_view name, const XdrMODE mode, const unsigned int read_flags, bool partition_agnostic);
#endif

} // namespace libMesh
