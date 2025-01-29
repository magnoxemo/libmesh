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



#ifndef LIBMESH_REMOTE_ELEM_H
#define LIBMESH_REMOTE_ELEM_H

// Local includes
#include "libmesh/elem.h"
#include "libmesh/libmesh_singleton.h"

// C++ includes
#include <limits>

#define remote_elem_error(func_name)                                    \
  do {                                                                  \
    std::stringstream message_stream;                                   \
    message_stream << "RemoteElem::" << func_name << " was called.  A RemoteElem is not a\n"  << \
                      "real Elem, merely a shim accessible on distributed meshes via\n"       << \
                      "Elem::neighbor_ptr(), to indicate when a neighbor exists but is not\n" << \
                      "local or ghosted on this processor.  Calling this method was an\n"     << \
                      "application or library error.  Investigate a stack trace for details.\n"; \
    libMesh::Threads::lock_singleton_spin_mutex();                      \
    libMesh::MacroFunctions::report_error(__FILE__, __LINE__, LIBMESH_DATE, LIBMESH_TIME, message_stream); \
    libMesh::Threads::unlock_singleton_spin_mutex();                    \
    LIBMESH_THROW(libMesh::LogicError(message_stream.str()));           \
  } while (0)


namespace libMesh
{

/**
 * In parallel meshes where a ghost element has neighbors which do
 * not exist on the local processor, the ghost element's neighbors
 * are set to point to the singleton RemoteElement instead.
 * Library code can then distinguish between such elements and
 * boundary elements (with nullptr neighbors).
 *
 * \author Roy H. Stogner
 * \date 2007
 * \brief Used by ParallelMesh to represent an Elem owned by another processor.
 */
class RemoteElem : public Elem,
                   public Singleton
{
public:

  /**
   * A unique \p id to distinguish remote element links
   */
  static const dof_id_type remote_elem_id = static_cast<dof_id_type>(-2);

  /**
   * Constructor. Private to force use of the \p create() member.
   */
private:
  RemoteElem () : Elem(0,
                       0,
                       nullptr,
                       _elemlinks_data,
                       nullptr)
  { this->set_id(remote_elem_id); }

public:

  RemoteElem (RemoteElem &&) = delete;
  RemoteElem (const RemoteElem &) = delete;
  RemoteElem & operator= (const RemoteElem &) = delete;
  RemoteElem & operator= (RemoteElem &&) = delete;

  /**
   * Sets remote_elem to nullptr.
   */
  virtual ~RemoteElem();

  /**
   * Return a reference to the global \p RemoteElem
   * singleton object.
   */
  static const Elem & create ();

  virtual Point master_point (const unsigned int /*i*/) const override
  { remote_elem_error("master_point"); return Point(); }

  virtual Node * & set_node (const unsigned int i) override
  { remote_elem_error("set_node"); return Elem::set_node(i); }

  /**
   * Don't hide Elem::key() defined in the base class.
   */
  using Elem::key;

  virtual dof_id_type key (const unsigned int) const override
  { remote_elem_error("key"); return 0; }

  virtual dof_id_type low_order_key (const unsigned int) const override
  { remote_elem_error("low_order_key"); return 0; }

  virtual unsigned int local_side_node(unsigned int /*side*/,
                                       unsigned int /*side_node*/) const override
  { remote_elem_error("local_side_node"); return 0; }

  virtual unsigned int local_edge_node(unsigned int /*side*/,
                                       unsigned int /*side_node*/) const override
  { remote_elem_error("local_edge_node"); return 0; }

  virtual bool is_remote () const override
  { return true; }

  virtual void connectivity(const unsigned int,
                            const IOPackage,
                            std::vector<dof_id_type> &) const override
  { remote_elem_error("connectivity"); }

  virtual ElemType type () const override
  { return REMOTEELEM; }

  virtual unsigned short dim () const override
  { remote_elem_error("dim"); return 0; }

  virtual unsigned int n_nodes () const override
  { remote_elem_error("n_nodes"); return 0; }

  virtual unsigned int n_sides () const override
  { remote_elem_error("n_sides"); return 0; }

  virtual unsigned int n_vertices () const override
  { remote_elem_error("n_vertices"); return 0; }

  virtual unsigned int n_edges () const override
  { remote_elem_error("n_edges"); return 0; }

  virtual unsigned int n_faces () const override
  { remote_elem_error("n_faces"); return 0; }

  virtual unsigned int n_children () const override
  { remote_elem_error("n_children"); return 0; }

  virtual bool is_vertex(const unsigned int) const override
  { remote_elem_error("is_vertex"); return false; }

  virtual bool is_edge(const unsigned int) const override
  { remote_elem_error("is_edge"); return false; }

  virtual bool is_face(const unsigned int) const override
  { remote_elem_error("is_face"); return false; }

  virtual bool is_node_on_side(const unsigned int,
                               const unsigned int) const override
  { remote_elem_error("is_node_on_side"); return false; }

  virtual std::vector<unsigned int> nodes_on_side(const unsigned int) const override
  { remote_elem_error("nodes_on_side"); return {0}; }

  virtual std::vector<unsigned int> nodes_on_edge(const unsigned int) const override
  { remote_elem_error("nodes_on_edge"); return {0}; }

  virtual std::vector<unsigned int> edges_adjacent_to_node(const unsigned int) const override
  { remote_elem_error("edges_adjacent_to_node"); return {0};}

  virtual std::vector<unsigned int> sides_on_edge(const unsigned int) const override
  { remote_elem_error("sides_on_edge"); return {0}; }

  virtual bool is_child_on_side(const unsigned int,
                                const unsigned int) const override
  { remote_elem_error("is_child_on_side"); return false; }

  virtual bool is_edge_on_side(const unsigned int,
                               const unsigned int) const override
  { remote_elem_error("is_edge_on_side"); return false; }

  virtual bool is_node_on_edge(const unsigned int,
                               const unsigned int) const override
  { remote_elem_error("is_node_on_edge"); return false; }

  virtual unsigned int n_sub_elem () const override
  { remote_elem_error("n_sub_elem"); return 0; }

  virtual std::unique_ptr<Elem> side_ptr (const unsigned int) override
  { remote_elem_error("side_ptr"); return std::unique_ptr<Elem>(); }

  virtual void side_ptr (std::unique_ptr<Elem> &,
                         const unsigned int) override
  { remote_elem_error("side_ptr"); }

  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int,
                                                bool) override
  { remote_elem_error("build_side_ptr"); return std::unique_ptr<Elem>(); }

  virtual void build_side_ptr (std::unique_ptr<Elem> &,
                               const unsigned int) override
  { remote_elem_error("build_side_ptr"); }

  virtual std::unique_ptr<Elem> build_edge_ptr (const unsigned int) override
  { remote_elem_error("build_edge_ptr"); return std::unique_ptr<Elem>(); }

  virtual void build_edge_ptr (std::unique_ptr<Elem> &, const unsigned int) override
  { remote_elem_error("build_edge_ptr"); }

  virtual Order default_order () const override
  { remote_elem_error("default_order"); return static_cast<Order>(1); }

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  virtual bool infinite () const override
  { remote_elem_error("infinite"); return false; }

#endif

#ifdef LIBMESH_ENABLE_AMR

  /**
   * Matrix that transforms the parents nodes into the children's
   * nodes.
   */
  virtual Real embedding_matrix (const unsigned int,
                                 const unsigned int,
                                 const unsigned int) const override
  { remote_elem_error("embedding_matrix"); return 0.; }

  LIBMESH_ENABLE_TOPOLOGY_CACHES;

#endif // LIBMESH_ENABLE_AMR

  virtual unsigned int n_permutations() const override
  { remote_elem_error("n_permutations"); return 0; }

  virtual void permute(unsigned int) override
  { remote_elem_error("permute"); }

  virtual void flip(BoundaryInfo *) override final
  { remote_elem_error("flip"); }

  virtual bool is_flipped() const override final
  { remote_elem_error("is_flipped"); return false; }

  virtual unsigned int center_node_on_side(const unsigned short) const override
  { remote_elem_error("center_node_on_side"); return invalid_uint; }

  virtual ElemType side_type (const unsigned int) const override
  { remote_elem_error("side_type"); return INVALID_ELEM; }

protected:

  /**
   * Data for link to (nullptr!) parent.
   */
  Elem * _elemlinks_data[1];
};

// Singleton RemoteElem
extern const RemoteElem * remote_elem;

} // namespace libMesh

#endif // LIBMESH_REMOTE_ELEM_H
