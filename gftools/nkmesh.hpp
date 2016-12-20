//
// Created by iskakoff on 16/12/16.
//

#ifndef GFTOOLS_NKMESH_H
#define GFTOOLS_NKMESH_H

#include <vector>
#include <array>


#include "almost_equal.hpp"
#include "grid_base.hpp"
#include "exceptions.hpp"

#include <boost/math/special_functions.hpp>

namespace gftools {

  /** A typedef for a point in the Brillouin zone. */
  template<int D>
  class BZPoint {
  public:
    BZPoint() : _domain_len(2*M_PI) {};
    BZPoint(const BZPoint& point): _domain_len(point._domain_len), _points(point._points){};
    BZPoint(std::array<real_type,D> points, double domain_len = 2*M_PI) : _points(points), _domain_len(domain_len){}
    BZPoint<D>& operator+(const BZPoint<D>& point) {
      for (int j = 0; j < D; ++j) {
        real_type out;
        out = _points[j] + real_type(point.points()[j]);
        out-= _domain_len*(almost_equal(out, _domain_len, num_io<double>::tolerance())?1.0:std::floor(out/_domain_len));
        _points[j] = out;
      }
      return *this;
    }

    const std::array < real_type, D > &points() const {
      return _points;
    }

    std::array < real_type, D > &points() {
      return _points;
    }
    double domain_len() const {return _domain_len;}
  private:
    std::array<real_type,D> _points;
    double _domain_len;
  };

  template <int D>
  BZPoint<D>& operator+(const BZPoint<D>& in, const BZPoint<D>& shift) {
    assert(in.domain_len() == shift.domain_len());
    BZPoint<D> p = in;
    for (int j = 0; j < D; ++j) {
      real_type out;
      out = p.points()[j] + real_type(shift.points()[j]);
      out-= in.domain_len()*(almost_equal(out, in.domain_len(), num_io<double>::tolerance())?1.0:std::floor(out/in.domain_len()));
      p.points()[j] = out;
    }
    return std::move(p);
  }

  /**
   * @brief N-dimensional k-mesh
   *
   * @author iskakoff
   */
  template<int N>
  class nkmesh : public grid_base < BZPoint<N> , nkmesh<N> > {
  public:
    ///typedef for the underlying grid
    typedef grid_base<BZPoint<N>, nkmesh<N> > base;
    ///typedef of the point type
    typedef BZPoint<N> point_type;
    typedef typename grid_base<BZPoint<N>, nkmesh<N> >::point point;
    ///this corresponds to this->vals_
    using grid_base<BZPoint<N>, nkmesh<N> >::vals_;
    ///constructor: expects number of discretization points and extent len of interval
    nkmesh(size_t n_points, real_type len = 2.0*M_PI);
    ///default constructor makes empty grid with zero points
    nkmesh():npoints_(0){};
    ///constructor from a vector of regularly spaced ints
    nkmesh(const std::vector<point_type> & in);


//    template <class Obj> auto integrate(const Obj &in) const ->decltype(in(vals_[0]));
//    template <class Obj> auto eval(Obj &in, real_type x) const ->decltype(in[0]);
//    template <class Obj> auto eval(Obj &in, point x) const ->decltype(in[0]) { return base::eval(in,x); }
//
    ///shift implements a periodic operator+ in three variants
    real_type shift(real_type in,real_type shift_arg) const {throw gftools::unimplemented_function_exception("Real value shift is unimplemented for multi-dimensional k-mesh;");};
    point shift(point in, real_type shift_arg) const {throw gftools::unimplemented_function_exception("Real value shift is unimplemented for multi-dimensional k-mesh;");}
    point shift(point in, point shift_arg) const { return in + shift_arg;}
  private:
    ///number of equidistantly spaced points in k-space
    size_t npoints_;
    ///extent of domain in k-space, usually 2*PI
    real_type domain_len_ = 2.0*M_PI;

    void get_k_point(std::array < real_type , N >& array, size_t i, int k) {
      array[k] = (i % npoints_)*domain_len_/npoints_;
      if (i>=npoints_) get_k_point(array, i/npoints_, k+1);
    }
  };

  template<int N>
  nkmesh<N>::nkmesh(size_t n_points, real_type len) : npoints_(n_points), domain_len_(len) {
    size_t total_points = boost::math::pow<N>(n_points);
    for (size_t i = 0; i < total_points; ++i) {
      std::array<real_type, N> points;
      get_k_point(points, i, 0);
      BZPoint<N> p(points);
      vals_.push_back(point(p, i));
    }
  }
  template<int N>
  nkmesh<N>::nkmesh(const std::vector < nkmesh::point_type > &in) : base(in), npoints_(in.size()) {}

}

#endif //GFTOOLS_NKMESH_H
