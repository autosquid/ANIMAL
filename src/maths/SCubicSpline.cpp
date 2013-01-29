// ===========================================================================
//    Description:  Cubic Spline
// ===========================================================================

#include "./SCubicSpline.h"
#include <iostream>

namespace CAPG
{
;


void CubicSpline3::
clear()
{
    m_xList.clear();
    m_yList.clear();
    m_zList.clear();

    m_xd2.clear();
    m_yd2.clear();
    m_zd2.clear();

    m_tList.clear();
} // clear()

void CubicSpline3::
reserve(size_t size)
{
    m_xList.reserve(size);
    m_yList.reserve(size);
    m_zList.reserve(size);
    m_tList.reserve(size);
} // reserve(size)

void CubicSpline3::
push_back(const Vector3 & p)
{
    m_xList.push_back(p.x());
    m_yList.push_back(p.y());
    m_zList.push_back(p.z());
} // push_back(p)

void CubicSpline3::
setupValue(const std::vector<Vector3> & pList)
{
    clear();
    //reserve(pList.size());
    push_back(pList.begin(), pList.end());
} // setupValue(pList)

void CubicSpline3::
setValue(size_t offset, const Vector3 & p)
{
    if (offset >= m_xList.size())  {
        return ;
    }

    m_xList[offset] = p.x();
    m_yList[offset] = p.y();
    m_zList[offset] = p.z();
} // setValue(offset, p)

void CubicSpline3::
push_backT(const CubicSpline3::ValueType t)
{
    m_tList.push_back(t);
} // push_backT(t)

void CubicSpline3::
setupT(const std::vector<CubicSpline3::ValueType> & tList)
{
    m_tList = tList;
} // setupT(tList)

void CubicSpline3::
setT(size_t offset, CubicSpline3::ValueType t)
{
    if (offset >= m_tList.size()) {
        return ;
    }

    m_tList[offset] = t;
} // setT(offset, t)

void CubicSpline3::
prepare()
{
    size_t s = m_xList.size();
    if (s < 2){
        return ;
    }

    m_xd2.resize(s);
    m_yd2.resize(s);
    m_zd2.resize(s);
    // generate t list, if needed.
    if (m_tList.size() != s) {
        m_tList.resize(s);
        ValueType delta = 1.0f/(s-1);
        for (size_t i = 0; i < s; ++i) {
            m_tList[i] = i*delta;
        }
        m_tList[0]   = 0;
        m_tList[s-1] = 1;
    }
    // compute d2
    compute_d2(m_tList, m_xList, m_xd2, 0, 0);
    compute_d2(m_tList, m_yList, m_yd2, 0, 0);
    compute_d2(m_tList, m_zList, m_zd2, 0, 0);

    return ;
} // prepare()

Vector3 CubicSpline3::
cubicSplinePoint(CubicSpline3::ValueType t) const
{
    if (t < 0  || t > 1) {
        return Vector3::ZERO;
    }

    if (m_xList.size() < 2 || m_xList.size() != m_tList.size()) {
        return Vector3::ZERO;
    }

    ValueType x = do_spline(m_tList, m_xList, m_xd2, t);
    ValueType y = do_spline(m_tList, m_yList, m_yd2, t);
    ValueType z = do_spline(m_tList, m_zList, m_zd2, t);

    return Vector3(x, y, z);
} // cubicSplinePoint(t)

Vector3 CubicSpline3::
derivative(CubicSpline3::ValueType t) const
{
    if (t < 0  || t > 1) {
        return Vector3::ZERO;
    }

    if (m_xList.size() < 2 || m_xList.size() != m_tList.size()) {
        return Vector3::ZERO;
    }

    ValueType x = do_derivative(m_tList, m_xList, m_xd2, t);
    ValueType y = do_derivative(m_tList, m_yList, m_yd2, t);
    ValueType z = do_derivative(m_tList, m_zList, m_zd2, t);

    return Vector3(x, y, z);
} // derivative(t)

Vector3 CubicSpline3::
derivative2nd(CubicSpline3::ValueType t) const
{
    if (t < 0  || t > 1) {
        return Vector3::ZERO;
    }

    if (m_xList.size() < 2 || m_xList.size() != m_tList.size()) {
        return Vector3::ZERO;
    }

    ValueType x = do_2nd_derivative(m_tList, m_xList, m_xd2, t);
    ValueType y = do_2nd_derivative(m_tList, m_yList, m_yd2, t);
    ValueType z = do_2nd_derivative(m_tList, m_zList, m_zd2, t);

    return Vector3(x, y, z);
} // derivative2nd(t)

void CubicSpline3::
compute_d2(const CubicSpline3::ValueList & tList,
           const CubicSpline3::ValueList & vList,
           CubicSpline3::ValueList & d2,
           CubicSpline3::ValueType d1,
           CubicSpline3::ValueType dn)
{
    size_t size = tList.size();

    ValueList u(size-1);
    d2.resize(size);

    if(d1 > 0.99e99){
        d2[0] = u[0] = 0;
    }else{
        d2[0] = -0.5;
        u[0] = static_cast<ValueType>( 
                ( 3.0/(tList[1]-tList[0]) ) * ( (vList[1]-vList[0])/(tList[1]-tList[0]) - d1 ) );
    }

    for(size_t i = 1; i < size-1; ++i){
        ValueType sig = (tList[i] - tList[i-1]) / (tList[i+1] - tList[i]);
        ValueType p = sig * d2[i-1] + 2;
        d2[i] = static_cast<ValueType>( (sig - 1.0) / p );
        u[i] = (vList[i+1]-vList[i])/(tList[i+1]-tList[i]) - (vList[i]-vList[i-1])/(tList[i]-tList[i-1]);
        u[i] = static_cast<ValueType>( (6.0*u[i]/(tList[i+1]-tList[i-1]) - sig*u[i-1]) / p );
    }

    ValueType qn(0.5), un(0);
    if(dn > 0.99e99){
        qn = un = 0;
    }else{
        qn = 0.5;
        un = static_cast<ValueType>( ( 3.0/(tList[size-1]-tList[size-2]) ) 
                * ( dn - (vList[size-1]-vList[size-2])/(tList[size-1]-tList[size-2]) ) );
    }

    d2[size-1] = static_cast<ValueType>( (un - qn*u[size-2]) / (qn*d2[size-2] + 1.0) );
    for(int i = size-2; i >= 0; --i){
        d2[i] = d2[i]*d2[i+1] + u[i];
    }

    return ;
} // compute_d2(t, v, d2, d1, dn)

CubicSpline3::ValueType CubicSpline3::
do_spline(const CubicSpline3::ValueList & tList,
        const CubicSpline3::ValueList & vList,
        const CubicSpline3::ValueList & d2,
        CubicSpline3::ValueType t) const
{
    ValueType a, b, h;
    size_t low, high;

    find_a_b_h_low_high(tList, vList, d2, t, a, b, h, low, high);

    return 	static_cast<ValueType>(
            a * vList[low]
            + b * vList[high]
            + ( (a*a*a - a)*d2[low] + (b*b*b - b)*d2[high] ) * (h*h)/6.0 );
} // do_spline(tList, vList, d2, t)

CubicSpline3::ValueType CubicSpline3::
do_derivative(const CubicSpline3::ValueList & tList,
              const CubicSpline3::ValueList & vList,
              const CubicSpline3::ValueList & d2,
              CubicSpline3::ValueType t) const
{
    ValueType a, b, h;
    size_t low, high;

    find_a_b_h_low_high(tList, vList, d2, t, a, b, h, low, high);

    ValueType tdy = vList[high] - vList[low];

    return static_cast<ValueType>(
            tdy/h + (h/2.0) * ((b*b-1/3.0)*d2[high]-(a*a-1/3.0)*d2[low]) );
} // do_derivative(tList, vList, d2, t)

CubicSpline3::ValueType CubicSpline3::
do_2nd_derivative(const CubicSpline3::ValueList & tList,
                  const CubicSpline3::ValueList & vList,
                  const CubicSpline3::ValueList & d2,
                  CubicSpline3::ValueType t) const
{
    ValueType a, b, h;
    size_t low, high;

    find_a_b_h_low_high(tList, vList, d2, t, a, b, h, low, high);

    return static_cast<ValueType> (a*d2[low] + b*d2[high]);
} // do_2nd_derivative(tList, vList, d2, t)

void CubicSpline3::
find_a_b_h_low_high(const CubicSpline3::ValueList & tList,
                    const CubicSpline3::ValueList & vList,
                    const CubicSpline3::ValueList & d2,
                    CubicSpline3::ValueType t,
                    CubicSpline3::ValueType & a,
                    CubicSpline3::ValueType & b,
                    CubicSpline3::ValueType & h,
                    size_t & low, size_t & high) const
{
    size_t size = tList.size();

    low = 0;
    high = size - 1;
    size_t mid(0);
    while( high - low > 1 ){
        mid = (high + low) >> 1;
        if(tList[mid] > t)
            high = mid;
        else
            low = mid;
    }
    h = tList[high] - tList[low];
    if(0.0 == h) throw ( "bad input" );
    a = (tList[high] - t) / h;
    b = (t - tList[low]) / h;

    return ;
} // do_derivative(tList, vList, d2, t)



} // namespace CAPG


