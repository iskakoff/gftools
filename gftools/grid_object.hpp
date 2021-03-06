#pragma once

#include "defaults.hpp"
#include "tools.hpp"
#include "grid_base.hpp"
#include "grid_tools.hpp"
#include "container.hpp"

namespace gftools {

/** A grid_object_base is a wrapper over container class, that stores data,
 * defined on multiple different grids.
 */
template <typename ContainerType, typename ...GridTypes> 
class grid_object_base;

template <typename ValueType, typename ... GridTypes>
using grid_object = grid_object_base<container<ValueType, sizeof...(GridTypes)>, GridTypes...>;

template <typename ValueType, typename ... GridTypes>
using grid_object_ref = grid_object_base<container_ref<ValueType, sizeof...(GridTypes)>, GridTypes...>;

 
template <typename ContainerType, typename ...GridTypes> 
class grid_object_base 
{
public:
    /// A typedef for the container. 
    typedef ContainerType container_type;
    /// A typedef for the values stored in the container. 
    typedef typename ContainerType::value_type value_type;
    /// A typedef for a tuple of grids. 
    typedef std::tuple<GridTypes...> grid_tuple;
    /// Tuple of grids traits. 
    typedef tools::grid_tuple_traits<grid_tuple> trs;
    /// Total number of grids. 
    static constexpr size_t N = sizeof...(GridTypes);

    /// A typedef for a function that gives the analytical value of the object, when it's not stored. 
    typedef typename std::remove_cv<typename tools::GridArgTypeExtractor<value_type, std::tuple<GridTypes...> >::f_type>::type function_type;
    /// A typedef for a function that gives the analytical value of the object, when it's not stored. 
    typedef typename tools::GridPointExtractor<value_type, std::tuple<GridTypes...> >::f_type point_function_type;
    /// A typedef for a tuple of grid points. 
    typedef typename trs::point_tuple point_tuple;
    /// A typedef for a tuple of grid point values. 
    typedef typename trs::arg_tuple arg_tuple;
    /// A typedef for an array of grid indices. 
    typedef typename trs::indices indices_t;

    class exPointMismatch : public std::exception { virtual const char* what() const throw() { return "Index mismatch."; }; };
    class exIOProblem : public std::exception { virtual const char* what() const throw(){return "IO problem.";} }; 
    class ex_wrong_index : public std::exception { virtual const char* what() const throw(){return "Index out of bounds";}}; 
protected:
    /// Grids on which the data is defined. 
    const std::tuple<GridTypes...> grids_;
    /// Cache data dimensions. 
    const indices_t dims_;
    /// Actual data - can be a container (data allocated) or a view (proxy to some other data). 
    container_type data_;
    /// This function returns the value of the object when the point is not in container. 
    function_type tail_;

public:
    // Constructors
    /// Constructs a grid object out of a tuple containing various grids. 
    grid_object_base( const std::tuple<GridTypes...> &grids);
    //grid_object_base( const std::tuple<GridTypes...> &grids, ContainerType&& in);
    grid_object_base( GridTypes... grids):grid_object_base(std::forward_as_tuple(grids...)){};
    /// Constructor of grids and data. 

    template <typename CType>
        grid_object_base( const std::tuple<GridTypes...> &grids, CType& data);
    grid_object_base( const std::tuple<GridTypes...> &grids, ContainerType&& data);

    /// Copy constructor. 
    grid_object_base( const grid_object_base<ContainerType, GridTypes...>& rhs);
    template <typename CType>
        grid_object_base( const grid_object_base<CType, GridTypes...>& rhs);
    /// Move constructor. 
    grid_object_base( grid_object_base<ContainerType, GridTypes...>&& rhs);
// Assignments
    grid_object_base& operator= (const grid_object_base & rhs);
    template <typename CType>
        grid_object_base& operator= (const grid_object_base<CType,GridTypes...> & rhs);
    grid_object_base& operator= (grid_object_base && rhs);
    grid_object_base& operator= (const value_type & rhs);

// Properties of gridobjects - dimensions, grids, etc
    std::tuple<GridTypes...> const& grids() const { return grids_; }
    /// Returns an Mth grid in grids_. 
    template<size_t M = 0> 
        auto grid() const -> const typename std::tuple_element<M, std::tuple<GridTypes...>>::type&
           { return std::get<M>(grids_); };
    size_t size() const { return data_.size(); }
    point_tuple points(indices_t in) const { return trs::points(in,grids_); }
    arg_tuple get_args(indices_t in) const { return trs::get_args(in, grids_); }
    arg_tuple get_args(point_tuple in) const { return trs::get_args(in, grids_); }
    indices_t get_indices(point_tuple in) const { return trs::get_indices(in, grids_); }
    /// Returns the data_ container. 
    container_type& data(){return data_;}
    container_type const& data() const {return data_;}
    function_type const& tail(){return tail_;}

    void set_tail(function_type&& f){tail_ = f; }
    void set_tail(function_type const& f){tail_ = f; }

// Global operations - reductions, shifts 
    /// Returns the complex conjugate of this object, if it's complex valued. 
    grid_object_base<ContainerType, GridTypes...> conj() const { return grid_object_base(grids_, data_.conj()); }
    /// Returns a norm of difference between two objects. 
    template <typename CT>
        real_type diff (const grid_object_base<CT, GridTypes...> &rhs, bool norm = true) const { 
            if (!trs::is_equal(this->grids(), rhs.grids())) throw (gftools::ex_generic("Objects are defined on different grids"));
            return data_.diff(rhs.data(), norm); } 
    /// Returns the sum of all elements in the container. 
    value_type sum() const { return data_.sum(); };
    /// Returns an object with arguments, shifted by the given values.
    template <typename ...ArgTypes> 
        typename std::enable_if<
            (std::is_convertible<std::tuple<typename std::remove_reference<ArgTypes>::type...>, arg_tuple>::value || std::is_convertible<std::tuple<ArgTypes...>, point_tuple>::value), 
            grid_object<value_type, GridTypes...>>::type shift(ArgTypes... args) const { return this->shift(std::forward_as_tuple(args...)); }
    template <typename ...ArgTypes> 
        typename std::enable_if<
            (std::is_convertible<std::tuple<ArgTypes...>, arg_tuple>::value || std::is_same<std::tuple<ArgTypes...>, point_tuple>::value), 
            grid_object<value_type, GridTypes...>>::type shift(const std::tuple<ArgTypes...>& arg_tuple) const ; 

// IO
    /// Save the data to the txt file. 
    void savetxt(const std::string& fname) const;
    /// Loads the data to the txt file. 
    void loadtxt(const std::string& fname, real_type tol = 1e-8);
    /// Dumps the object to the stream. 
    template <typename CT, class ...GridTypes2> friend std::ostream& operator<<(std::ostream& lhs, const grid_object_base<CT,GridTypes2...> &in);

// Access operators
    /// Returns element number i, which corresponds to (*_grid)[i]. 
    auto operator[](size_t i)->decltype(data_[i]) { return data_[i]; };
    auto operator[](size_t i) const -> const typename container_type::under_ref_type { return data_[i]; };
    auto operator[](typename std::tuple_element<0,grid_tuple>::type::point in)->decltype(data_[0]) { assert(in.index() < grid<0>().size()); return data_[in.index()]; };
    auto operator[](typename std::tuple_element<0,grid_tuple>::type::point in) const -> const typename container_type::under_ref_type { assert(in.index() < grid<0>().size()); return data_[in.index()]; };
    //template <size_t M> value_type& operator[](const std::array<size_t,M>& in);
    /// Returns tail_(in). 
    value_type tail_eval(const arg_tuple& in) const { return tuple_tools::unfold_tuple(tail_, in); }
    template <typename ... ArgTypes>
        typename std::enable_if<std::is_convertible<std::tuple<ArgTypes...>, arg_tuple>::value, value_type>::type 
            tail_eval(ArgTypes...in) { return tail_(in...); };
    /// Return the value by grid values. 
    value_type& get(const point_tuple& in) { return data_(get_indices(in)); }

    value_type& operator()(const point_tuple& in) { return data_(get_indices(in)); }
    value_type operator()(const point_tuple& in) const { 
            try { return data_(get_indices(in)); } 
            catch (gftools::ex_generic) { return this->tail_eval(in); };
        }

    template <int M=N>
    typename std::enable_if<(M>1 ), value_type>::type  operator()(arg_tuple in) const
    	{  
            
            point_tuple x = trs::find_nearest(in, grids_); 
            if (gftools::tools::is_float_equal<arg_tuple>(x,in)) {  return (*this)(x); } // FIXME with expression templates 
            else { return this->tail_eval(in); };
        // FIXME # warning grid_object::operator() doesn't provide interpolation for D>=2 
        }
    template <int M=N>
    typename std::enable_if<(M==1 ), value_type>::type  operator()(arg_tuple in) const
    	{ try { return std::get<0>(grids_).eval(data_, std::get<0>(in)); } catch (gftools::gftools_exception) { return this->tail_eval(in); } }

    template <int M=N>
    typename std::enable_if<(M==1 ), value_type>::type
    operator()(typename std::tuple_element<0,grid_tuple>::type::point in) const { 
        if (in.index() < grid().size() && tools::is_float_equal(in.value(), grid().points()[in.index()])) { return data_[in.index()]; } 
        return this->operator()(in.value()); 
        };

    template <int M=N>
    typename std::enable_if<(M==1 ), value_type&>::type
    operator()(typename std::tuple_element<0,grid_tuple>::type::point in)
        { return data_[in.index()]; };


    template <typename ...ArgTypes>
    	typename std::enable_if<std::is_convertible<std::tuple<ArgTypes...>, point_tuple>::value && sizeof...(ArgTypes) != 1 && sizeof...(ArgTypes)==N
    							, value_type&>::type
        operator()(ArgTypes... in) { return (*this)(point_tuple(std::forward_as_tuple(in...))); }

    template <typename ...ArgTypes>
    	typename std::enable_if<std::is_convertible<std::tuple<ArgTypes...>, point_tuple>::value && sizeof...(ArgTypes) != 1 && sizeof...(ArgTypes)==N
    							,value_type>::type
        operator()(ArgTypes... in) const { return (*this)(point_tuple(std::forward_as_tuple(in...))); }



    template <typename ...ArgTypes>
        	typename std::enable_if<!std::is_convertible<std::tuple<ArgTypes...>, point_tuple>::value &&
        							sizeof...(ArgTypes)!=1 && (sizeof...(ArgTypes)==N)
        							, value_type>::type
            operator()(ArgTypes... in) const { return (*this)(arg_tuple(std::forward_as_tuple(in...))); }

    template <typename ArgType>
    	typename std::enable_if<std::is_same<std::tuple<ArgType>, arg_tuple>::value, value_type>::type
    	 operator()(ArgType in) const { try { return std::get<0>(grids_).eval(data_, in); } 
            catch (gftools::gftools_exception) { return tail_(in); }; }

    /*template <typename ArgType>
        	typename std::enable_if<std::is_same<std::tuple<ArgType>, arg_tuple>::value, value_type>::type
        	 operator()(ArgType in) { return std::get<0>(grids_).eval(data_, in); }
*/
    /// Return value of grid_object. eval methods of grids are used, so interpolation is done if provided with grids


// Fill values
    /// Fills the container with a provided function. 
    void fill_function(const function_type &in);
    void fill_point_function(const point_function_type &in);
    template <typename F> 
        typename std::enable_if<std::is_convertible<typename std::remove_cv<F>::type, point_function_type>::value, void>::type fill(const F& f) {
            this->fill_point_function(f); }
    template <typename F> 
        typename std::enable_if<std::is_convertible<typename std::remove_cv<F>::type, function_type>::value && 
                                !std::is_convertible<typename std::remove_cv<F>::type, point_function_type>::value, void>::type fill(const F& f) { 
            this->fill_function(f); }
    void fill(const std::function<value_type(arg_tuple)>& in) { this->fill_function(tools::extract_tuple_f(in)); }
    void fill(const std::function<value_type(point_tuple)>& in) { this->fill_point_function(tools::extract_tuple_f(in)); }
    /// A shortcut for fill method. 
    //template <typename ...ArgTypes> grid_object_base& operator= (const std::function<value_type(ArgTypes...)> &);
    /// Same as operator=, but allows for non-equal grids. Slow. Uses analytic function to provide missing values. 
    template <typename CType2>
    grid_object_base& copy_interpolate(const grid_object_base<CType2, GridTypes...> &rhs);

    typedef grid_object<value_type, GridTypes...> gobj_t;

// Math (should be removed to an external algebra class). 
    grid_object_base& operator*= (const grid_object_base & rhs);
    grid_object_base& operator*= (const value_type& rhs);
    grid_object<value_type, GridTypes...> operator* (const grid_object_base & rhs) const { gobj_t out(this->grids()); out.data() = this->data(); out*=rhs; return out; }  
    grid_object<value_type, GridTypes...> operator* (const value_type & rhs) const { gobj_t out(*this); out*=rhs; return out; }
    grid_object_base& operator+= (const grid_object_base & rhs);
    grid_object_base& operator+= (const value_type& rhs);
    grid_object<value_type, GridTypes...> operator+ (const grid_object_base & rhs) const { gobj_t out(*this); out+=rhs; return out; }
    grid_object<value_type, GridTypes...> operator+ (const value_type & rhs) const { gobj_t out(*this); out+=rhs; return out; }
    grid_object_base& operator-= (const grid_object_base & rhs);
    grid_object_base& operator-= (const value_type& rhs);
    grid_object<value_type, GridTypes...> operator- (const grid_object_base & rhs) const { gobj_t out(*this); out-=rhs; return out; }
    grid_object<value_type, GridTypes...> operator- (const value_type & rhs) const { gobj_t out(*this); out-=rhs; return out; }
    grid_object_base& operator/= (const grid_object_base & rhs);
    grid_object_base& operator/= (const value_type& rhs);
    grid_object<value_type, GridTypes...> operator/ (const grid_object_base & rhs) const { gobj_t out(*this); out/=rhs; return out; }
    grid_object<value_type, GridTypes...> operator/ (const value_type & rhs) const { gobj_t out(*this); out/=rhs; return out; }
    friend grid_object<value_type, GridTypes...> operator* (const value_type & lhs, const grid_object_base & rhs) {return rhs*lhs;};
    friend grid_object<value_type, GridTypes...> operator+ (const value_type & lhs, const grid_object_base & rhs) {return rhs+lhs;};
    friend grid_object<value_type, GridTypes...> operator- (const value_type & lhs, const grid_object_base & rhs) {return rhs*(-1.0)+lhs;};
    friend grid_object<value_type, GridTypes...> operator/ (const value_type & lhs, const grid_object_base & rhs) {grid_object_base out(rhs); out=lhs; return out/rhs;};
};

template <typename GridObjectType>
GridObjectType loadtxt(std::string const& fname, double tol = 1e-8);


/*
/// A helper recursive template utility to extract and set data from the container.
    template <size_t Nc, typename CT, typename ArgType1, typename ...ArgTypes> struct containerExtractor { 
        typedef typename CT::value_type value_type;
        /// Gets the data by values.
        static value_type get(CT &data, const std::tuple<GridTypes...> &grids, const ArgType1& arg1, const ArgTypes&... args);
        static value_type& get_ref(CT &data, const std::tuple<GridTypes...> &grids, const ArgType1& arg1, const ArgTypes&... args);
             };
    /// Specialization of containerExtractor for 1-dim container.
    template <typename CT, typename ArgType1> struct containerExtractor<1,CT,ArgType1> {
        typedef typename CT::value_type value_type;
        static value_type get(CT &data, const std::tuple<GridTypes...> &grids, const ArgType1& arg1); 
        static value_type get(CT &data, const std::tuple<GridTypes...> &grids, const std::tuple<ArgType1>& arg1); 

        static value_type& get_ref(CT &data, const std::tuple<GridTypes...> &grids, const ArgType1& arg1); 
        static value_type& get_ref(CT &data, const std::tuple<GridTypes...> &grids, const std::tuple<ArgType1>& arg1); 
        };
*/

} // end of namespace gftools

#include "grid_object.hxx"

