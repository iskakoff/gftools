//
// Created by iskakoff on 06/01/17.
//

#ifndef GFTOOLS_TAIL_HPP
#define GFTOOLS_TAIL_HPP

#include <cassert>
#include <array>

#include "defaults.hpp"
#include "grid_object.hpp"
#include "index_grid.hpp"
#include "matsubara_grid.hpp"

/**
 * @brief tail class
 *
 * @author iskakoff
 */
namespace gftools {

  template <size_t D, class ArrayType>
  struct tail_container_base: container_base<real_type, D+1, ArrayType> {
    using typename container_base<real_type , D+1, ArrayType>::boost_t;
    tail_container_base(const tail_container_base<D,ArrayType> &in) : tail_container_base(in, in.min_tail_order(), in.tail_orders()) {}
    tail_container_base(const tail_container_base<D,ArrayType> &in, int min_tail_order, int tail_orders) : container_base<real_type, D+1, ArrayType>(in),
                                                                                                           _min_tail_order(min_tail_order), _tail_orders(tail_orders) {}
    tail_container_base(const boost_t& in) : tail_container_base(in, -1, 0){}
    tail_container_base(const boost_t& in, int min_tail_order, int tail_orders) : container_base<real_type, D+1, ArrayType>(in),
                                                                                  _min_tail_order(min_tail_order), _tail_orders(tail_orders) {}
    int min_tail_order() const {
      return _min_tail_order;
    }
    int tail_orders() const {
      return _tail_orders;
    }
  private:
    /// lowest non-zero tail order
    int _min_tail_order;
    /// largest non-zero tail order
    int _tail_orders;

  };

  template<size_t D>
  using tail_container_ref = tail_container_base<D, boost::multi_array_ref<real_type, D+1> >;

  /**
   * @tparam D tail object grid dimension. Since we need additional grid for tail coefficient base class grid has actually D+1 dimension
   */
  template <size_t D>
  struct tail_container: tail_container_base<D, boost::multi_array<real_type, D+1> > {
    using tail_container_base<D, boost::multi_array<real_type, D+1> >::boost_container_;
    template <typename CT>
    tail_container(const tail_container_base<D, CT>& in) : tail_container_base<D,typename boost::multi_array<real_type, D+1>>(in.boost_container_(), in.min_tail_order(), in.tail_orders()) {};
  };

  /**
   *
   * @tparam GridTypes - definition of grid for high frequency tail object
   */
  template<typename value_type, typename LeadingType, typename ... GridTypes>
  class tail : public grid_object_base<tail_container<sizeof...(GridTypes)>, index_grid, GridTypes...> {
  public:
    /// Since we have additional index grid for tail coefficients we need to redefine grid tuple.
    typedef std::tuple<GridTypes...> grid_tuple;
    typedef std::tuple<LeadingType, GridTypes...> call_grid_tuple;
    typedef std::tuple<typename LeadingType::point, typename GridTypes::point...> function_tuple;
//    typedef tools::grid_tuple_traits<call_grid_tuple> trs;
    typedef tools::grid_tuple_traits<call_grid_tuple> trs2;
    /// A typedef for a tuple of grid points.
    typedef typename trs2::point_tuple call_point_tuple;
    /// A typedef for a tuple of grid point values.
    typedef typename trs2::arg_tuple call_arg_tuple;
    /// A typedef for an array of grid indices.
    typedef typename trs2::indices call_indices_t;

    ///
    using base_class = grid_object_base<tail_container<sizeof...(GridTypes)>, index_grid, GridTypes...>;
    /// import type declarations from base class
    using typename base_class::container_type;
    using base_value_type = typename base_class::value_type;
    using typename base_class::point_tuple;
    /// import field and method definition from base class
    using base_class::N;
    using base_class::get_indices;
    using base_class::grids;
    /// constructors
    tail(GridTypes&...grids) : tail(std::forward_as_tuple(index_grid(0), grids...)) {};
    tail(const std::tuple<GridTypes...>& grids) : grid_object_base<container_type, index_grid, GridTypes...>(std::tuple_cat(std::make_tuple(index_grid(0)), grids)){};

    tail(const std::tuple<GridTypes...>& grids, container_type& c) : grid_object_base<container_type, index_grid, GridTypes...>(std::tuple_cat(std::make_tuple(index_grid(c.tail_orders())), grids), c){};

    template <typename ...ArgTypes>
    typename std::enable_if<std::is_convertible<std::tuple<ArgTypes...>, call_point_tuple>::value && sizeof...(ArgTypes) != 1 && sizeof...(ArgTypes)==N
      ,value_type>::type
    operator()(ArgTypes... in) const {
      const gftools::index_grid &c_grid = std::get<0>(grids());
      std::vector<base_value_type> v;
      /// get tail coefficients
      for(auto & c_index: c_grid.points()) {
        base_value_type x = operator()(c_index, in...);
        v.push_back(x);
      }
      const std::tuple <ArgTypes...> &tuple = std::make_tuple(in...);
      return f(std::get<0>(tuple), v);
    }

  private:
    template<typename ...RestArgs>
    base_value_type operator()(const index_grid::point& c, const typename LeadingType::point& m, RestArgs...rest) const {
//      std::cout << typeid(rest).name() << '\n';
      base_value_type res = base_class::operator()(c, rest...);
      return res;
    }
    /// Matsubara frequency tail
    value_type f(const point_base<complex_type>& iwn, const std::vector<base_value_type>& c) const {
      value_type iwnsq=iwn.value()*iwn.value();
      return (c.size()>0 ? c[0]/iwn.value() + (c.size()>1 ? c[1]/(iwnsq) + (c.size()>2 ? c[2]/(iwn.value()*iwnsq) : 0.0) : 0.0 ) : 0.0);
    };
    /// TODO: implement
    /// Imaginary time tail
    real_type f_tau(real_type beta, real_type tau, real_type c1, real_type c2, real_type c3) const {
      return -0.5*c1 + (c2*0.25)*(-beta+2.*tau) + (c3*0.25)*(beta*tau-tau*tau);
    };
  };
}

#endif //GFTOOLS_TAIL_HPP
