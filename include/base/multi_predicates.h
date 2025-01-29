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

#ifndef LIBMESH_MULTI_PREDICATES_H
#define LIBMESH_MULTI_PREDICATES_H

// Local includes
#include "libmesh/libmesh.h" // libMesh::invalid_uint
#include "libmesh/single_predicates.h"

// C++ includes
#include <vector>
#include <memory>
#include <iterator> // std::make_move_iterator
#include <initializer_list> // std::initializer_list

namespace libMesh {
class Elem;
}

namespace libMesh
{

// Forward declarations
class BoundaryInfo;

/**
 * This namespace defines several multi_predicates which are used by
 * the element and node iterators.  These classes are not in general
 * used by the user, although they could be.
 *
 * \author John W. Peterson
 * \date 2004
 */
namespace Predicates
{

// Empty place-holder base class for multi_predicates
struct multi_predicate {};

// Alias declaration for the predicate pointer we use
template <typename T>
using pred_ptr = std::unique_ptr<predicate<T>>;

// This class represents a generic combination of more than one predicate.
// It is meant to be derived from to actually be used.
template <typename T>
struct abstract_multi_predicate : multi_predicate
{
  // virtual destructor.
  virtual ~abstract_multi_predicate() = default;

  // operator= (perform deep copy of entries in _predicates vector
  abstract_multi_predicate & operator=(const abstract_multi_predicate & rhs)
  {
    // Copy over the information from the rhs.
    this->deep_copy(rhs);

    return *this;
  }

  // operator() checks all the predicates in the vector.
  virtual bool operator()(const T & it) const
  {
    for (const auto & pred : _predicates)
      {
        libmesh_assert (pred);

        if (!(*pred)(it))
          return false;
      }

    return true;
  }

protected:
  // Do not instantiate the base class.
  abstract_multi_predicate() = default;

  // Construct from a vector of single predicates
  abstract_multi_predicate(std::vector<pred_ptr<T>> && predicates)
    : _predicates(std::move(predicates))
  {}

  // Copy constructor.
  abstract_multi_predicate(const abstract_multi_predicate & rhs)
  {
    this->deep_copy(rhs);
  }

  // The deep_copy function is used by both the op= and
  // copy constructors.  This function uses the default (empty)
  // copy constructor for the predicate class.
  void deep_copy(const abstract_multi_predicate & rhs)
  {
    // First clear out the predicates vector
    _predicates.clear();

    for (const auto & p : rhs._predicates)
      _predicates.push_back(p->clone());
  }

  // Predicates to be evaluated.
  std::vector<pred_ptr<T>> _predicates;
};

/**
 * Helper object for creating a std::vector from a std::initializer_list
 * https://stackoverflow.com/questions/46737054/vectorunique-ptra-using-initialization-list
 */
template<class T>
struct movable_il
{
  /**
   * Construct from rvalue reference of type T
   */
  movable_il(T && in) : t(std::move(in)) {}

  /**
   * Construct from rvalue reference of type U, using forwarding
   */
  template <typename U>
  movable_il(U && in): t(std::forward<U>(in)) {}

  /**
   * Return an rvalue reference to ourself
   */
  operator T() const&& { return std::move(t); }

  mutable T t;
};

/**
 * Helper function that creates a std::vector from an initializer_list of movable_il objects
 */
template<class T>
std::vector<T> make_vec( std::initializer_list< movable_il<T> > il )
{
  std::vector<T> r( std::make_move_iterator(il.begin()), std::make_move_iterator(il.end()) );
  return r;
}


/**
 * Used to iterate over nullptr entries in a container.
 */
template <typename T>
struct IsNull : abstract_multi_predicate<T>
{
  IsNull() : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<is_null<T>>()
  })) {}
};



/**
 * Used to iterate over non-nullptr entries in a container.
 */
template <typename T>
struct NotNull : abstract_multi_predicate<T>
{
  NotNull() : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>()
  })) {}
};



/**
 * Used to iterate over non-nullptr, active entries in a container.
 */
template <typename T>
struct Active : abstract_multi_predicate<T>
{
  Active() : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<active<T>>()
  })) {}
};



/**
 * Used to iterate over non-nullptr, inactive entries in a container.
 */
template <typename T>
struct NotActive : abstract_multi_predicate<T>
{
  NotActive() : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<not_active<T>>()
  })) {}
};




/**
 * Used to iterate over non-nullptr, entries that have children (i.e. are
 * ancestors) in a container.
 */
template <typename T>
struct Ancestor : abstract_multi_predicate<T>
{
  Ancestor() : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<ancestor<T>>()
  })) {}
};




/**
 * Used to iterate over non-nullptr, entries that have no children (i.e. are not
 * ancestors) in a container.
 */
template <typename T>
struct NotAncestor : abstract_multi_predicate<T>
{
  NotAncestor() : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<not_ancestor<T>>()
  })) {}
};




/**
 * Used to iterate over non-nullptr, subactive entries (i.e. has no
 * active children) in a container.
 */
template <typename T>
struct SubActive : abstract_multi_predicate<T>
{
  SubActive() : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<subactive<T>>()
  })) {}
};




/**
 * Used to iterate over non-nullptr, non-subactive entries (i.e. has one
 * or more active children) in a container.
 */
template <typename T>
struct NotSubActive : abstract_multi_predicate<T>
{
  NotSubActive() : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<not_subactive<T>>()
  })) {}
};



/**
 * Used to iterate over non-nullptr, local entries (i.e. owned by the
 * current processor) in a container.
 */
template <typename T>
struct Local : abstract_multi_predicate<T>
{
  Local(processor_id_type my_pid) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<pid<T>>(my_pid)
  })) {}
};


/**
 * Used to iterate over non-nullptr, semi-local entries (i.e. are not
 * subactive and have are owned by an attached processor) in a
 * container.
 *
 * FIXME: This is not currently safe to use on adaptively-refined
 * grids, it should be added back when Elem::is_semilocal() has been
 * patched to not require the Elem to be active.
 */
// template <typename T>
// struct SemiLocal : abstract_multi_predicate<T>
// {
//   SemiLocal(processor_id_type my_pid) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
//   {
//     std::make_unique<not_null<T>>(),
//     std::make_unique<not_subactive<T>>(),
//     std::make_unique<semilocal_pid<T>>(my_pid)
//   })) {}
// };


/**
 * Used to iterate over non-nullptr, active, non sub-active, semi-local
 * elements in a container.
 */
template <typename T>
struct ActiveSemiLocal : abstract_multi_predicate<T>
{
  ActiveSemiLocal(processor_id_type my_pid) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<active<T>>(),
    std::make_unique<not_subactive<T>>(),
    std::make_unique<semilocal_pid<T>>(my_pid),
  })) {}
};


/**
 * Used to iterate over non-nullptr, face-local entries (i.e. are not
 * subactive and are on or have a neighbor on processor my_pid) in a
 * container.
 */
template <typename T>
struct FaceLocal : abstract_multi_predicate<T>
{
  FaceLocal(processor_id_type my_pid) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<not_subactive<T>>(),
    std::make_unique<facelocal_pid<T>>(my_pid)
  })) {}
};



/**
 * Used to iterate over non-nullptr, non-local entries in a
 * container.
 */
template <typename T>
struct NotLocal : abstract_multi_predicate<T>
{
  NotLocal(processor_id_type my_pid) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<not_pid<T>>(my_pid)
  })) {}
};


/**
 * Used to iterate over non-nullptr, active, non-local entries in a
 * container.
 */
template <typename T>
struct ActiveNotLocal : abstract_multi_predicate<T>
{
  ActiveNotLocal(processor_id_type my_pid) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<active<T>>(),
    std::make_unique<not_pid<T>>(my_pid)
  })) {}
};


/**
 * Used to iterate over non-nullptr, elements of a given geometric type.
 */
template <typename T>
struct Type : abstract_multi_predicate<T>
{
  Type(ElemType type) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<elem_type<T>>(type)
  })) {}
};



/**
 * Used to iterate over non-nullptr, active elements of a given geometric type.
 */
template <typename T>
struct ActiveType : abstract_multi_predicate<T>
{
  ActiveType(ElemType type) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<active<T>>(),
    std::make_unique<elem_type<T>>(type)
  })) {}
};



#ifdef LIBMESH_ENABLE_AMR
/**
 * Used to iterate over non-nullptr, elements with a given refinement
 * flag.
 */
template <typename T>
struct Flagged : abstract_multi_predicate<T>
{
  Flagged(unsigned char rflag) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<flagged<T>>(rflag)
  })) {}
};



/**
 * Used to iterate over non-nullptr, elements with a given refinement
 * flag belonging to a given processor.
 */
template <typename T>
struct FlaggedPID : abstract_multi_predicate<T>
{
  FlaggedPID(unsigned char rflag, processor_id_type proc_id) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<flagged<T>>(rflag),
    std::make_unique<pid<T>>(proc_id)
  })) {}
};

#endif // LIBMESH_ENABLE_AMR




/**
 * Used to iterate over non-nullptr, active elements owned by a given
 * processor.
 */
template <typename T>
struct ActivePID : abstract_multi_predicate<T>
{
  ActivePID(processor_id_type proc_id) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<active<T>>(),
    std::make_unique<pid<T>>(proc_id)
  })) {}
};





/**
 * Used to iterate over non-nullptr, active, local elements owned by a
 * given processor.
 */
template <typename T>
struct ActiveLocal : abstract_multi_predicate<T>
{
  ActiveLocal(processor_id_type my_pid) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<active<T>>(),
    std::make_unique<pid<T>>(my_pid)
  })) {}
};





/**
 * Used to iterate over non-nullptr elements owned by a given processor.
 */
template <typename T>
struct PID : abstract_multi_predicate<T>
{
  PID(processor_id_type proc_id) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<pid<T>>(proc_id)
  })) {}
};



/**
 * Used to iterate over non-nullptr elements on the boundary with a given
 * ID.
 */
template <typename T>
struct BID : abstract_multi_predicate<T>
{
  BID(boundary_id_type bndry_id, const BoundaryInfo & bndry_info) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<bid<T>>(bndry_id, bndry_info)
  })) {}
};



/**
 * Used to iterate over non-nullptr elements on the boundary.
 */
template <typename T>
struct BND : abstract_multi_predicate<T>
{
  BND(const BoundaryInfo & bndry_info) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<bnd<T>>(bndry_info)
  })) {}
};



/**
 * Used to iterate over non-nullptr elements *not* owned by a given
 * processor.
 */
template <typename T>
struct NotPID : abstract_multi_predicate<T>
{
  NotPID(processor_id_type proc_id) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<not_pid<T>>(proc_id)
  })) {}
};



/**
 * Used to iterate over non-nullptr elements of a specified (refinement) level.
 */
template <typename T>
struct Level : abstract_multi_predicate<T>
{
  Level(unsigned int l) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<level<T>>(l)
  })) {}
};



/**
 * Used to iterate over non-nullptr elements *not* of a specified
 * (refinement) level.
 */
template <typename T>
struct NotLevel : abstract_multi_predicate<T>
{
  NotLevel(unsigned int l) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<not_level<T>>(l)
  })) {}
};



/**
 * Used to iterate over non-nullptr local elements with a specified
 * (refinement) level.
 */
template <typename T>
struct LocalLevel : abstract_multi_predicate<T>
{
  LocalLevel(processor_id_type my_pid,
             unsigned int l) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<pid<T>>(my_pid),
    std::make_unique<level<T>>(l)
  })) {}
};



/**
 * Used to iterate over non-nullptr local elements *not* of a specified
 * (refinement) level.
 */
template <typename T>
struct LocalNotLevel : abstract_multi_predicate<T>
{
  LocalNotLevel(processor_id_type my_pid,
                unsigned int l) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<pid<T>>(my_pid),
    std::make_unique<not_level<T>>(l)
  })) {}
};



/**
 * Used to iterate over non-nullptr, active elements which are on the
 * boundary.
 */
template <typename T>
struct ActiveOnBoundary : abstract_multi_predicate<T>
{
  ActiveOnBoundary() : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<active<T>>(),
    std::make_unique<null_neighbor<T>>()
  })) {}
};



/**
 * Used to iterate over the sides of an element which are on the
 * boundary of the Mesh.
 */
template <typename T>
struct BoundarySide : abstract_multi_predicate<T>
{
  BoundarySide() : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<boundary_side<T>>()
  })) {}
};



/**
 * Used to iterate over non-nullptr, active elements with a given PID on
 * a given subdomain.
 */
template <typename T>
struct ActiveLocalSubdomain : abstract_multi_predicate<T>
{
  ActiveLocalSubdomain(processor_id_type my_pid,
                       subdomain_id_type subdomain_id) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<active<T>>(),
    std::make_unique<pid<T>>(my_pid),
    std::make_unique<subdomain<T>>(subdomain_id)
  })) {}
};



/**
 * Used to iterate over non-nullptr, active elements on a given
 * subdomain.
 */
template <typename T>
struct ActiveSubdomain : abstract_multi_predicate<T>
{
  ActiveSubdomain(subdomain_id_type subdomain_id) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<active<T>>(),
    std::make_unique<subdomain<T>>(subdomain_id)
  })) {}
};



/**
 * Used to iterate over non-nullptr, active elements whose
 * subdomains are in a user-specified set.
 */
template <typename T>
struct ActiveSubdomainSet : abstract_multi_predicate<T>
{
  ActiveSubdomainSet(std::set<subdomain_id_type> sset) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<active<T>>(),
    std::make_unique<subdomain_set<T>>(sset)
  })) {}
};



/**
 * Used to iterate over non-nullptr, active elements with a given PID
 * whose subdomains are in a user-specified set.
 */
template <typename T>
struct ActiveLocalSubdomainSet : abstract_multi_predicate<T>
{
  ActiveLocalSubdomainSet(processor_id_type my_pid,
                          std::set<subdomain_id_type> sset) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<active<T>>(),
    std::make_unique<pid<T>>(my_pid),
    std::make_unique<subdomain_set<T>>(sset)
  })) {}
};






/**
 * Used to iterate over non-nullptr elements not owned by a given
 * processor but semi-local to that processor, i.e. ghost elements.
 */
template <typename T>
struct Ghost : abstract_multi_predicate<T>
{
  Ghost(processor_id_type my_pid) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<active<T>>(),
    std::make_unique<not_pid<T>>(my_pid),
    std::make_unique<semilocal_pid<T>>(my_pid)
  })) {}
};



/**
 * Used to iterate over elements where solutions indexed by a given
 * DofMap are evaluable for a given variable var_num.
 */
template <typename T>
struct Evaluable: abstract_multi_predicate<T>
{
  Evaluable(const DofMap & dof_map,
            unsigned int var_num = libMesh::invalid_uint) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<active<T>>(),
    std::make_unique<evaluable<T>>(dof_map, var_num)
  })) {}
};



/**
 * Used to iterate over elements where solutions indexed by a given
 * vector of DofMaps are evaluable for all variables
 */
template <typename T>
struct MultiEvaluable: abstract_multi_predicate<T>
{
  MultiEvaluable(const std::vector<const DofMap *> & dof_maps) : abstract_multi_predicate<T>(make_vec<pred_ptr<T>>(
  {
    std::make_unique<not_null<T>>(),
    std::make_unique<active<T>>(),
    std::make_unique<multi_evaluable<T>>(dof_maps)
  })) {}
};

}


} // namespace libMesh

#endif // LIBMESH_MULTI_PREDICATES_H
