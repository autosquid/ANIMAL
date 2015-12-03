//////////////////////////////////////////////////////
// file: SEuler.cpp @ 2008-8-28 by Zhang Xiang
// defines of the class SEuler
// SEuler is a class ...
//////////////////////////////////////////////////////

// INCLUDES //////////////////////////////////////////
#include "SEuler.h"
#include "SMatrix3.h"

// DECLARES //////////////////////////////////////////

// DEFINES ///////////////////////////////////////////
namespace CAPG{

	const Euler Euler::ZERO(0.0, 0.0, 0.0);

	void Euler::fromRotateMatrix(const Matrix3& rotMat, RotateOrder _rotateOrder){
		this->rotateOrder = _rotateOrder;
		switch(rotateOrder)
		{
		case ZXY:
			m_fx = Math::ASin(rotMat.arr()[7]);
			m_fy = Math::ATan2(-rotMat.arr()[6], rotMat.arr()[8]);
			m_fz = Math::ATan2(-rotMat.arr()[1], rotMat.arr()[4]);
			break;
		case XYZ:
            assert(0);
			break;
        case XZY:
            assert(0);
            break;
        case YXZ:
            assert(0);
            break;
        case YZX:
            assert(0);
            break;
        case ZYX:
            assert(0);
            break;
		default:
			break;
		}
	}

	void Euler::toRotateMatrix(Matrix3& rmRot) const
	{
		switch(rotateOrder)
		{
		case ZXY:
			rmRot.fromEulerAnglesZXY(m_fz, m_fx, m_fy);
			break;
		case XYZ:
            rmRot.fromEulerAnglesXYZ(m_fx, m_fy, m_fz);
			break;
        case XZY:
            rmRot.fromEulerAnglesXZY(m_fx, m_fz, m_fy);
            break;
        case YXZ:
            rmRot.fromEulerAnglesYXZ(m_fy, m_fx, m_fz);
            break;
        case YZX:
            rmRot.fromEulerAnglesYZX(m_fy, m_fz, m_fx);
            break;
        case ZYX:
            rmRot.fromEulerAnglesZYX(m_fz, m_fy, m_fx);
            break;
		default:
			break;
		}
		
	}

} // namespace CAPG
