#ifndef ZZ_CAPG_MATHS_CUBIC_SPLINE_2008032719__H___
#define ZZ_CAPG_MATHS_CUBIC_SPLINE_2008032719__H___


#include "./SMath.h"
#include "./SVector3.h"
#include <vector>

namespace CAPG
{
;

/**
 * @brief Cubic Spline interpolation
 *
 * 	It refers to William H. Press's Numerical Recipes in C++ 2nd Edition
 *	Usage:  1 - input data by calling push_back#() or set#()
 *		    2 - do prepare()
 *		    3 - call Y = line_spline(X) to get the interpolated value
 */

class CubicSpline3
{
public:
    typedef Real                        ValueType;
    typedef std::vector<ValueType>      ValueList;
public:
    void clear();
    void reserve(size_t size);
    
    // input points
    void push_back(const Vector3 & p);
    template <typename ItType>
    void push_back(ItType begin, ItType end);

    void setupValue(const std::vector<Vector3> & pList);
    template <typename ItType>
    void setupValue(ItType begin, ItType end);

    void setValue(size_t offset, const Vector3 & p);

    // input time
    void push_backT(const ValueType t);
    template <typename ItType>
    void push_backT(ItType begin, ItType end);

    void setupT(const std::vector<ValueType> & tList);
    template <typename ItType>
    void setupT(ItType begin, ItType end);

    void setT(size_t offset, ValueType t);
public:
    void prepare();
    Vector3 cubicSplinePoint(ValueType t) const;
    Vector3 derivative(ValueType t) const;
    Vector3 derivative2nd(ValueType t) const;

    Vector3 velocity(ValueType t) const {return derivative(t);}
    Vector3 acceleration(ValueType t) const {return derivative2nd(t);}
private:
    void compute_d2(const ValueList & tList, const ValueList & vList,
                    ValueList & d2List, ValueType d1, ValueType dn);
    ValueType do_spline(const ValueList & tList, const ValueList & vList,
                        const ValueList & d2, ValueType t) const;
    ValueType do_derivative(const ValueList & tList, const ValueList & vList,
                            const ValueList & d2, ValueType t) const;
    ValueType do_2nd_derivative(const ValueList & tList, const ValueList & vList,
                                const ValueList & d2, ValueType t) const;
    void find_a_b_h_low_high(const ValueList & tList, const ValueList & vList,
                             const ValueList & d2, ValueType t,
                             ValueType & a, ValueType & b, ValueType & h,
                             size_t & low, size_t & high) const;
private:
    ValueList           m_xList, m_yList, m_zList;
    ValueList           m_xd2, m_yd2, m_zd2;
    ValueList           m_tList;
}; // class CubicSpline3_T

template <typename ItType>
void CubicSpline3::push_back(ItType begin, ItType end)
{
    for (ItType iv = begin; iv != end; ++iv) {
        push_back(*iv);
    }
} // push_back(begin, end)

template <typename ItType>
void CubicSpline3::setupValue(ItType begin, ItType end)
{
    clear();
    push_back(begin, end);
} // setupValue(begin, end)

template <typename ItType>
void CubicSpline3::push_backT(ItType begin, ItType end)
{
    for (ItType it = begin; it != end; ++it) {
        push_backT(*it);
    }
} // push_backT(begin, end)

template <typename ItType>
void CubicSpline3::setupT(ItType begin, ItType end)
{
    m_tList.clear();
    push_backT(begin, end);
} // setupT(begin, end)

} // namespace CAPG

#endif

