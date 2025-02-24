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


// Local includes
#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/number_lookups.h"
#include "libmesh/enum_to_string.h"

namespace
{
using namespace libMesh;

// Compute the static coefficients for an element
void hermite_compute_coefs(const Elem * elem, std::vector<std::vector<Real>> & dxdxi

#ifdef DEBUG
                           , std::vector<Real> & dydxi, std::vector<Real> & dzdeta, std::vector<Real> & dxdzeta,
                           std::vector<Real> & dzdxi, std::vector<Real> & dxdeta, std::vector<Real> & dydzeta
#endif
                           )
{
  const FEFamily mapping_family = FEMap::map_fe_type(*elem);
  const FEType map_fe_type(elem->default_order(), mapping_family);

  // Note: we explicitly don't consider the elem->p_level() when
  // computing the number of mapping shape functions.
  const int n_mapping_shape_functions =
    FEInterface::n_shape_functions (map_fe_type, /*extra_order=*/0, elem);

  static const Point dofpt[2] = {Point(-1,-1,-1), Point(1,1,1)};

  FEInterface::shape_deriv_ptr shape_deriv_ptr =
    FEInterface::shape_deriv_function(map_fe_type, elem);

  for (int p = 0; p != 2; ++p)
    {
      dxdxi[0][p] = 0;
      dxdxi[1][p] = 0;
      dxdxi[2][p] = 0;
#ifdef DEBUG
      dydxi[p] = 0;
      dzdeta[p] = 0;
      dxdzeta[p] = 0;
      dzdxi[p] = 0;
      dxdeta[p] = 0;
      dydzeta[p] = 0;
#endif
      for (int i = 0; i != n_mapping_shape_functions; ++i)
        {
          const Real ddxi = shape_deriv_ptr
            (map_fe_type, elem, i, 0, dofpt[p], /*add_p_level=*/false);
          const Real ddeta = shape_deriv_ptr
            (map_fe_type, elem, i, 1, dofpt[p], /*add_p_level=*/false);
          const Real ddzeta = shape_deriv_ptr
            (map_fe_type, elem, i, 2, dofpt[p], /*add_p_level=*/false);

          // dxdeta, dxdzeta, dydxi, dydzeta, dzdxi, dzdeta should all
          // be 0!
          const Point & point_i = elem->point(i);
          dxdxi[0][p] += point_i(0) * ddxi;
          dxdxi[1][p] += point_i(1) * ddeta;
          dxdxi[2][p] += point_i(2) * ddzeta;
#ifdef DEBUG
          dydxi[p] += point_i(1) * ddxi;
          dzdeta[p] += point_i(2) * ddeta;
          dxdzeta[p] += point_i(0) * ddzeta;
          dzdxi[p] += point_i(2) * ddxi;
          dxdeta[p] += point_i(0) * ddeta;
          dydzeta[p] += point_i(1) * ddzeta;
#endif
        }

      // No singular elements!
      libmesh_assert(dxdxi[0][p]);
      libmesh_assert(dxdxi[1][p]);
      libmesh_assert(dxdxi[2][p]);
      // No non-rectilinear or non-axis-aligned elements!
#ifdef DEBUG
      libmesh_assert_less (std::abs(dydxi[p]), TOLERANCE);
      libmesh_assert_less (std::abs(dzdeta[p]), TOLERANCE);
      libmesh_assert_less (std::abs(dxdzeta[p]), TOLERANCE);
      libmesh_assert_less (std::abs(dzdxi[p]), TOLERANCE);
      libmesh_assert_less (std::abs(dxdeta[p]), TOLERANCE);
      libmesh_assert_less (std::abs(dydzeta[p]), TOLERANCE);
#endif
    }
}



Real hermite_bases_3D (std::vector<unsigned int> & bases1D,
                       const std::vector<std::vector<Real>> & dxdxi,
                       const Order & o,
                       unsigned int i)
{
  bases1D.clear();
  bases1D.resize(3,0);
  Real coef = 1.0;

  unsigned int e = o-2;

  // Nodes
  if (i < 64)
    {
      switch (i / 8)
        {
        case 0:
          break;
        case 1:
          bases1D[0] = 1;
          break;
        case 2:
          bases1D[0] = 1;
          bases1D[1] = 1;
          break;
        case 3:
          bases1D[1] = 1;
          break;
        case 4:
          bases1D[2] = 1;
          break;
        case 5:
          bases1D[0] = 1;
          bases1D[2] = 1;
          break;
        case 6:
          bases1D[0] = 1;
          bases1D[1] = 1;
          bases1D[2] = 1;
          break;
        case 7:
          bases1D[1] = 1;
          bases1D[2] = 1;
          break;
        default:
          libmesh_error_msg("Invalid basis node = " << i/8);
        }

      unsigned int basisnum = i%8;
      switch (basisnum) // DoF type
        {
        case 0: // DoF = value at node
          coef = 1.0;
          break;
        case 1: // DoF = x derivative at node
          coef = dxdxi[0][bases1D[0]];
          bases1D[0] += 2; break;
        case 2: // DoF = y derivative at node
          coef = dxdxi[1][bases1D[1]];
          bases1D[1] += 2; break;
        case 3: // DoF = xy derivative at node
          coef = dxdxi[0][bases1D[0]] * dxdxi[1][bases1D[1]];
          bases1D[0] += 2; bases1D[1] += 2; break;
        case 4: // DoF = z derivative at node
          coef = dxdxi[2][bases1D[2]];
          bases1D[2] += 2; break;
        case 5: // DoF = xz derivative at node
          coef = dxdxi[0][bases1D[0]] * dxdxi[2][bases1D[2]];
          bases1D[0] += 2; bases1D[2] += 2; break;
        case 6: // DoF = yz derivative at node
          coef = dxdxi[1][bases1D[1]] * dxdxi[2][bases1D[2]];
          bases1D[1] += 2; bases1D[2] += 2; break;
        case 7: // DoF = xyz derivative at node
          coef = dxdxi[0][bases1D[0]] * dxdxi[1][bases1D[1]] * dxdxi[2][bases1D[2]];
          bases1D[0] += 2; bases1D[1] += 2; bases1D[2] += 2; break;
        default:
          libmesh_error_msg("Invalid basisnum = " << basisnum);
        }
    }
  // Edges
  else if (i < 64 + 12*4*e)
    {
      unsigned int basisnum = (i - 64) % (4*e);
      switch ((i - 64) / (4*e))
        {
        case 0:
          bases1D[0] = basisnum / 4 + 4;
          bases1D[1] = basisnum % 4 / 2 * 2;
          bases1D[2] = basisnum % 2 * 2;
          if (basisnum % 4 / 2)
            coef *= dxdxi[1][0];
          if (basisnum % 2)
            coef *= dxdxi[2][0];
          break;
        case 1:
          bases1D[0] = basisnum % 4 / 2 * 2 + 1;
          bases1D[1] = basisnum / 4 + 4;
          bases1D[2] = basisnum % 2 * 2;
          if (basisnum % 4 / 2)
            coef *= dxdxi[0][1];
          if (basisnum % 2)
            coef *= dxdxi[2][0];
          break;
        case 2:
          bases1D[0] = basisnum / 4 + 4;
          bases1D[1] = basisnum % 4 / 2 * 2 + 1;
          bases1D[2] = basisnum % 2 * 2;
          if (basisnum % 4 / 2)
            coef *= dxdxi[1][1];
          if (basisnum % 2)
            coef *= dxdxi[2][0];
          break;
        case 3:
          bases1D[0] = basisnum % 4 / 2 * 2;
          bases1D[1] = basisnum / 4 + 4;
          bases1D[2] = basisnum % 2 * 2;
          if (basisnum % 4 / 2)
            coef *= dxdxi[0][0];
          if (basisnum % 2)
            coef *= dxdxi[2][0];
          break;
        case 4:
          bases1D[0] = basisnum % 4 / 2 * 2;
          bases1D[1] = basisnum % 2 * 2;
          bases1D[2] = basisnum / 4 + 4;
          if (basisnum % 4 / 2)
            coef *= dxdxi[0][0];
          if (basisnum % 2)
            coef *= dxdxi[1][0];
          break;
        case 5:
          bases1D[0] = basisnum % 4 / 2 * 2 + 1;
          bases1D[1] = basisnum % 2 * 2;
          bases1D[2] = basisnum / 4 + 4;
          if (basisnum % 4 / 2)
            coef *= dxdxi[0][1];
          if (basisnum % 2)
            coef *= dxdxi[1][0];
          break;
        case 6:
          bases1D[0] = basisnum % 4 / 2 * 2 + 1;
          bases1D[1] = basisnum % 2 * 2 + 1;
          bases1D[2] = basisnum / 4 + 4;
          if (basisnum % 4 / 2)
            coef *= dxdxi[0][1];
          if (basisnum % 2)
            coef *= dxdxi[1][1];
          break;
        case 7:
          bases1D[0] = basisnum % 4 / 2 * 2;
          bases1D[1] = basisnum % 2 * 2 + 1;
          bases1D[2] = basisnum / 4 + 4;
          if (basisnum % 4 / 2)
            coef *= dxdxi[0][0];
          if (basisnum % 2)
            coef *= dxdxi[1][1];
          break;
        case 8:
          bases1D[0] = basisnum / 4 + 4;
          bases1D[1] = basisnum % 4 / 2 * 2;
          bases1D[2] = basisnum % 2 * 2 + 1;
          if (basisnum % 4 / 2)
            coef *= dxdxi[1][0];
          if (basisnum % 2)
            coef *= dxdxi[2][1];
          break;
        case 9:
          bases1D[0] = basisnum % 4 / 2 * 2 + 1;
          bases1D[1] = basisnum / 4 + 4;
          bases1D[2] = basisnum % 2 * 2;
          if (basisnum % 4 / 2)
            coef *= dxdxi[0][1];
          if (basisnum % 2)
            coef *= dxdxi[2][1];
          break;
        case 10:
          bases1D[0] = basisnum / 4 + 4;
          bases1D[1] = basisnum % 4 / 2 * 2 + 1;
          bases1D[2] = basisnum % 2 * 2 + 1;
          if (basisnum % 4 / 2)
            coef *= dxdxi[1][1];
          if (basisnum % 2)
            coef *= dxdxi[2][1];
          break;
        case 11:
          bases1D[0] = basisnum % 4 / 2 * 2;
          bases1D[1] = basisnum / 4 + 4;
          bases1D[2] = basisnum % 2 * 2 + 1;
          if (basisnum % 4 / 2)
            coef *= dxdxi[0][0];
          if (basisnum % 2)
            coef *= dxdxi[2][1];
          break;
        default:
          libmesh_error_msg("Invalid basis node = " << (i - 64) / (4*e));
        }
    }
  // Faces
  else if (i < 64 + 12*4*e + 6*2*e*e)
    {
      unsigned int basisnum = (i - 64 - 12*4*e) % (2*e*e);
      switch ((i - 64 - 12*4*e) / (2*e*e))
        {
        case 0:
          bases1D[0] = square_number_column[basisnum / 2];
          bases1D[1] = square_number_row[basisnum / 2];
          bases1D[2] = basisnum % 2 * 2;
          if (basisnum % 2)
            coef *= dxdxi[2][0];
          break;
        case 1:
          bases1D[0] = square_number_column[basisnum / 2];
          bases1D[1] = basisnum % 2 * 2;
          bases1D[2] = square_number_row[basisnum / 2];
          if (basisnum % 2)
            coef *= dxdxi[1][0];
          break;
        case 2:
          bases1D[0] = basisnum % 2 * 2 + 1;
          bases1D[1] = square_number_column[basisnum / 2];
          bases1D[2] = square_number_row[basisnum / 2];
          if (basisnum % 2)
            coef *= dxdxi[0][1];
          break;
        case 3:
          bases1D[0] = square_number_column[basisnum / 2];
          bases1D[1] = basisnum % 2 * 2 + 1;
          bases1D[2] = square_number_row[basisnum / 2];
          if (basisnum % 2)
            coef *= dxdxi[1][1];
          break;
        case 4:
          bases1D[0] = basisnum % 2 * 2;
          bases1D[1] = square_number_column[basisnum / 2];
          bases1D[2] = square_number_row[basisnum / 2];
          if (basisnum % 2)
            coef *= dxdxi[0][0];
          break;
        case 5:
          bases1D[0] = square_number_column[basisnum / 2];
          bases1D[1] = square_number_row[basisnum / 2];
          bases1D[2] = basisnum % 2 * 2 + 1;
          if (basisnum % 2)
            coef *= dxdxi[2][1];
          break;
        default:
          libmesh_error_msg("Invalid basis node = " << (i - 64 - 12*4*e) / (2*e*e));
        }
    }
  // Interior
  else
    {
      unsigned int basisnum = i - 64 - 12*4*e;
      bases1D[0] = cube_number_column[basisnum] + 4;
      bases1D[1] = cube_number_row[basisnum] + 4;
      bases1D[2] = cube_number_page[basisnum] + 4;
    }

  // No singular elements
  libmesh_assert(coef);
  return coef;
}


} // end anonymous namespace


namespace libMesh
{


LIBMESH_DEFAULT_VECTORIZED_FE(3,HERMITE)


template <>
Real FE<3,HERMITE>::shape(const Elem * elem,
                          const Order order,
                          const unsigned int i,
                          const Point & p,
                          const bool add_p_level)
{
  libmesh_assert(elem);

  std::vector<std::vector<Real>> dxdxi(3, std::vector<Real>(2, 0));

#ifdef DEBUG
  std::vector<Real> dydxi(2), dzdeta(2), dxdzeta(2);
  std::vector<Real> dzdxi(2), dxdeta(2), dydzeta(2);
#endif //DEBUG

  hermite_compute_coefs(elem, dxdxi
#ifdef DEBUG
                        , dydxi, dzdeta, dxdzeta, dzdxi, dxdeta, dydzeta
#endif
                        );

  const ElemType type = elem->type();

  const Order totalorder =
    order + add_p_level*elem->p_level();

  switch (totalorder)
    {
      // 3rd-order tricubic Hermite functions
    case THIRD:
      {
        switch (type)
          {
          case HEX8:
          case HEX20:
          case HEX27:
            {
              libmesh_assert_less (i, 64);

              std::vector<unsigned int> bases1D;

              Real coef = hermite_bases_3D(bases1D, dxdxi, totalorder, i);

              return coef *
                FEHermite<1>::hermite_raw_shape(bases1D[0],p(0)) *
                FEHermite<1>::hermite_raw_shape(bases1D[1],p(1)) *
                FEHermite<1>::hermite_raw_shape(bases1D[2],p(2));
            }
          default:
            libmesh_error_msg("ERROR: Unsupported element type " << Utility::enum_to_string(type));
          }
      }
      // by default throw an error
    default:
      libmesh_error_msg("ERROR: Unsupported polynomial order " << totalorder);
    }
}


template <>
Real FE<3,HERMITE>::shape(const ElemType,
                          const Order,
                          const unsigned int,
                          const Point &)
{
  libmesh_error_msg("Hermite elements require the real element \nto construct gradient-based degrees of freedom.");
  return 0.;
}

template <>
Real FE<3,HERMITE>::shape(const FEType fet,
                          const Elem * elem,
                          const unsigned int i,
                          const Point & p,
                          const bool add_p_level)
{
  return FE<3,HERMITE>::shape(elem, fet.order, i, p, add_p_level);
}



template <>
Real FE<3,HERMITE>::shape_deriv(const Elem * elem,
                                const Order order,
                                const unsigned int i,
                                const unsigned int j,
                                const Point & p,
                                const bool add_p_level)
{
  libmesh_assert(elem);
  libmesh_assert (j == 0 || j == 1 || j == 2);

  std::vector<std::vector<Real>> dxdxi(3, std::vector<Real>(2, 0));

#ifdef DEBUG
  std::vector<Real> dydxi(2), dzdeta(2), dxdzeta(2);
  std::vector<Real> dzdxi(2), dxdeta(2), dydzeta(2);
#endif //DEBUG

  hermite_compute_coefs(elem, dxdxi
#ifdef DEBUG
                        , dydxi, dzdeta, dxdzeta, dzdxi, dxdeta, dydzeta
#endif
                        );

  const ElemType type = elem->type();

  const Order totalorder =
    order + add_p_level*elem->p_level();

  switch (totalorder)
    {
      // 3rd-order tricubic Hermite functions
    case THIRD:
      {
        switch (type)
          {
          case HEX8:
          case HEX20:
          case HEX27:
            {
              libmesh_assert_less (i, 64);

              std::vector<unsigned int> bases1D;

              Real coef = hermite_bases_3D(bases1D, dxdxi, totalorder, i);

              switch (j) // Derivative type
                {
                case 0:
                  return coef *
                    FEHermite<1>::hermite_raw_shape_deriv(bases1D[0],p(0)) *
                    FEHermite<1>::hermite_raw_shape(bases1D[1],p(1)) *
                    FEHermite<1>::hermite_raw_shape(bases1D[2],p(2));
                  break;
                case 1:
                  return coef *
                    FEHermite<1>::hermite_raw_shape(bases1D[0],p(0)) *
                    FEHermite<1>::hermite_raw_shape_deriv(bases1D[1],p(1)) *
                    FEHermite<1>::hermite_raw_shape(bases1D[2],p(2));
                  break;
                case 2:
                  return coef *
                    FEHermite<1>::hermite_raw_shape(bases1D[0],p(0)) *
                    FEHermite<1>::hermite_raw_shape(bases1D[1],p(1)) *
                    FEHermite<1>::hermite_raw_shape_deriv(bases1D[2],p(2));
                  break;
                default:
                  libmesh_error_msg("Invalid shape function derivative j = " << j);
                }

            }
          default:
            libmesh_error_msg("ERROR: Unsupported element type " << Utility::enum_to_string(type));
          }
      }
      // by default throw an error
    default:
      libmesh_error_msg("ERROR: Unsupported polynomial order " << totalorder);
    }
}


template <>
Real FE<3,HERMITE>::shape_deriv(const ElemType,
                                const Order,
                                const unsigned int,
                                const unsigned int,
                                const Point &)
{
  libmesh_error_msg("Hermite elements require the real element \nto construct gradient-based degrees of freedom.");
  return 0.;
}


template <>
Real FE<3,HERMITE>::shape_deriv(const FEType fet,
                                const Elem * elem,
                                const unsigned int i,
                                const unsigned int j,
                                const Point & p,
                                const bool add_p_level)
{
  return FE<3,HERMITE>::shape_deriv(elem, fet.order, i, j, p, add_p_level);
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES


template <>
Real FE<3,HERMITE>::shape_second_deriv(const Elem * elem,
                                       const Order order,
                                       const unsigned int i,
                                       const unsigned int j,
                                       const Point & p,
                                       const bool add_p_level)
{
  libmesh_assert(elem);

  std::vector<std::vector<Real>> dxdxi(3, std::vector<Real>(2, 0));

#ifdef DEBUG
  std::vector<Real> dydxi(2), dzdeta(2), dxdzeta(2);
  std::vector<Real> dzdxi(2), dxdeta(2), dydzeta(2);
#endif //DEBUG

  hermite_compute_coefs(elem, dxdxi
#ifdef DEBUG
                        , dydxi, dzdeta, dxdzeta, dzdxi, dxdeta, dydzeta
#endif
                        );

  const ElemType type = elem->type();

  const Order totalorder =
    order + add_p_level*elem->p_level();

  switch (totalorder)
    {
      // 3rd-order tricubic Hermite functions
    case THIRD:
      {
        switch (type)
          {
          case HEX8:
          case HEX20:
          case HEX27:
            {
              libmesh_assert_less (i, 64);

              std::vector<unsigned int> bases1D;

              Real coef = hermite_bases_3D(bases1D, dxdxi, totalorder, i);

              switch (j) // Derivative type
                {
                case 0:
                  return coef *
                    FEHermite<1>::hermite_raw_shape_second_deriv(bases1D[0],p(0)) *
                    FEHermite<1>::hermite_raw_shape(bases1D[1],p(1)) *
                    FEHermite<1>::hermite_raw_shape(bases1D[2],p(2));
                  break;
                case 1:
                  return coef *
                    FEHermite<1>::hermite_raw_shape_deriv(bases1D[0],p(0)) *
                    FEHermite<1>::hermite_raw_shape_deriv(bases1D[1],p(1)) *
                    FEHermite<1>::hermite_raw_shape(bases1D[2],p(2));
                  break;
                case 2:
                  return coef *
                    FEHermite<1>::hermite_raw_shape(bases1D[0],p(0)) *
                    FEHermite<1>::hermite_raw_shape_second_deriv(bases1D[1],p(1)) *
                    FEHermite<1>::hermite_raw_shape(bases1D[2],p(2));
                  break;
                case 3:
                  return coef *
                    FEHermite<1>::hermite_raw_shape_deriv(bases1D[0],p(0)) *
                    FEHermite<1>::hermite_raw_shape(bases1D[1],p(1)) *
                    FEHermite<1>::hermite_raw_shape_deriv(bases1D[2],p(2));
                  break;
                case 4:
                  return coef *
                    FEHermite<1>::hermite_raw_shape(bases1D[0],p(0)) *
                    FEHermite<1>::hermite_raw_shape_deriv(bases1D[1],p(1)) *
                    FEHermite<1>::hermite_raw_shape_deriv(bases1D[2],p(2));
                  break;
                case 5:
                  return coef *
                    FEHermite<1>::hermite_raw_shape(bases1D[0],p(0)) *
                    FEHermite<1>::hermite_raw_shape(bases1D[1],p(1)) *
                    FEHermite<1>::hermite_raw_shape_second_deriv(bases1D[2],p(2));
                  break;
                default:
                  libmesh_error_msg("Invalid shape function derivative j = " << j);
                }

            }
          default:
            libmesh_error_msg("ERROR: Unsupported element type " << Utility::enum_to_string(type));
          }
      }
      // by default throw an error
    default:
      libmesh_error_msg("ERROR: Unsupported polynomial order " << totalorder);
    }
}


template <>
Real FE<3,HERMITE>::shape_second_deriv(const ElemType,
                                       const Order,
                                       const unsigned int,
                                       const unsigned int,
                                       const Point &)
{
  libmesh_error_msg("Hermite elements require the real element \nto construct gradient-based degrees of freedom.");
  return 0.;
}

template <>
Real FE<3,HERMITE>::shape_second_deriv(const FEType fet,
                                       const Elem * elem,
                                       const unsigned int i,
                                       const unsigned int j,
                                       const Point & p,
                                       const bool add_p_level)
{
  return FE<3,HERMITE>::shape_second_deriv(elem, fet.order, i, j, p, add_p_level);
}

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

} // namespace libMesh
