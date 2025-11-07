#include "Matrix3x3.hpp"
#include <stdexcept>

#define TOL 1e-6
#define PI 3.14159265358979323846

// ------------------ Vec3 -------------------------

double Vec3::Dot(const Vec3& a, const Vec3& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 Vec3::Cross(const Vec3& a, const Vec3& b)
{
    return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
}

double Vec3::Norm() const
{
    return std::sqrt(Dot(*this, *this));
}

Vec3 Vec3::Normalize() const
{
    double n = Norm();
    if (n == 0) throw std::invalid_argument("normalize: zero vector");
    return { x / n, y / n, z / n };
}

// ------------------ Matrix3x3 ---------------------

Matrix3x3 Matrix3x3::Identity()
{
    Matrix3x3 I;
    I.At(0, 0) = 1; I.At(1, 1) = 1; I.At(2, 2) = 1;
    return I;
}

Vec3 Matrix3x3::Multiply(const Vec3& x, OpsCounter* op) const
{
    // y = A * x
    Vec3 y;
    y.x = At(0, 0) * x.x + At(0, 1) * x.y + At(0, 2) * x.z;
    y.y = At(1, 0) * x.x + At(1, 1) * x.y + At(1, 2) * x.z;
    y.z = At(2, 0) * x.x + At(2, 1) * x.y + At(2, 2) * x.z;
    if (op) { op->IncMul(9); op->IncAdd(6); }
    return y;
}

Matrix3x3 Matrix3x3::Multiply(const Matrix3x3& B, OpsCounter* op) const
{
    Matrix3x3 C{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double s = 0.0;
            for (int k = 0; k < 3; ++k) {
                s += At(i, k) * B.At(k, j);
            }
            C.At(i, j) = s;
        }
    }
    if (op) { op->IncMul(27); op->IncAdd(18); }
    return C;
}

double Matrix3x3::Det() const
{
    double pos= m[0] * m[4] * m[8] + m[3] * m[7] * m[2] + m[1] * m[5] * m[6];
   
    double ret= m[2] * m[4] * m[6] + m[1] * m[3] * m[8] + m[5] * m[7] * m[0];
 
    return pos-ret;
}

Matrix3x3 Matrix3x3::Transposed() const
{
    Matrix3x3 Transposed;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            Transposed.At(j, i) = this->At(i, j);

        }
   }
    return Transposed;
}

double Matrix3x3::Trace() const
{
    double trace =At(0,0)+At(1,1)+At(2,2);
    return trace;
}

bool Matrix3x3::IsRotation() const
{
    
    if (std::abs(this->Det()- 1.0)>TOL) {
        return false;
    }  
    return true;
}

Matrix3x3 Matrix3x3::RotationAxisAngle(const Vec3& u_in, double phi)
{
    //TODO
    return {};
}

Vec3 Matrix3x3::Rotate(const Vec3& v, OpsCounter* op) const
{
    //TODO
    return {};
}

void Matrix3x3::ToAxisAngle(Vec3& axis, double& angle) const
{
    //TODO
}

Matrix3x3 Matrix3x3::FromEulerZYX(double yaw, double pitch, double roll)
{
    //TODO
    return {};
}

void Matrix3x3::ToEulerZYX(double& yaw, double& pitch, double& roll) const
{
    //TODO
}