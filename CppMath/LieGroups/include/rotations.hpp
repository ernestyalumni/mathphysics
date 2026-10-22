#pragma once
// Goldstein §§4.2–4.7. Active right-handed rotations, column vectors.
// A frame coordinate conversion body -> inertial uses the same R; Goldstein A=R^T.
#include "lie_groups.hpp"
#include <algorithm>
#include <numbers>

namespace mathphysics {
template<class T> Matrix<T,3,3> axis_rotation(std::size_t axis,const T& c,const T& s) {
  if(axis>=3)throw std::out_of_range("rotation axis");
  auto m=Matrix<T,3,3>::identity();const auto j=(axis+1)%3,k=(axis+2)%3;
  m(j,j)=c;m(k,k)=c;m(j,k)=-s;m(k,j)=s;return m;
}
template<std::floating_point T> So3<T> axis_rotation(std::size_t axis,T angle) {
  if(!std::isfinite(angle))throw std::invalid_argument("nonfinite angle");
  return So3<T>::from_matrix(axis_rotation(axis,std::cos(angle),std::sin(angle)));
}
// Chronological intrinsic rotations (about the successively moving axes).
template<std::floating_point T> So3<T> from_euler_zxz(T phi,T theta,T psi) {
  return axis_rotation(2,phi).compose(axis_rotation(0,theta)).compose(axis_rotation(2,psi));
}
template<std::floating_point T> So3<T> from_euler_zyz(T phi,T theta,T psi) {
  return axis_rotation(2,phi).compose(axis_rotation(1,theta)).compose(axis_rotation(2,psi));
}
template<std::floating_point T> So3<T> from_yaw_pitch_roll(T psi,T theta,T phi) {
  return axis_rotation(2,psi).compose(axis_rotation(1,theta)).compose(axis_rotation(0,phi));
}
// Half-angle lifts retain the quaternion sign (unlike a principal inverse).
template<std::floating_point T> Su2<T> quaternion_from_euler_zxz(T phi,T theta,T psi) {
  if(!std::isfinite(phi)||!std::isfinite(theta)||!std::isfinite(psi))throw std::invalid_argument("nonfinite Euler angle");
  const T c=std::cos(theta/T{2}),s=std::sin(theta/T{2}),ca=std::cos(phi/T{2}),sa=std::sin(phi/T{2}),cb=std::cos(psi/T{2}),sb=std::sin(psi/T{2});
  return Su2<T>::from_unit_quaternion({c*(ca*cb-sa*sb),s*(ca*cb+sa*sb),s*(sa*cb-ca*sb),c*(sa*cb+ca*sb)});
}
template<std::floating_point T> Su2<T> quaternion_from_euler_zyz(T phi,T theta,T psi) {
  if(!std::isfinite(phi)||!std::isfinite(theta)||!std::isfinite(psi))throw std::invalid_argument("nonfinite Euler angle");
  const T c=std::cos(theta/T{2}),s=std::sin(theta/T{2}),ca=std::cos(phi/T{2}),sa=std::sin(phi/T{2}),cb=std::cos(psi/T{2}),sb=std::sin(psi/T{2});
  return Su2<T>::from_unit_quaternion({c*(ca*cb-sa*sb),-s*(sa*cb-ca*sb),s*(ca*cb+sa*sb),c*(sa*cb+ca*sb)});
}
template<std::floating_point T> Su2<T> quaternion_from_yaw_pitch_roll(T psi,T theta,T phi) {
  if(!std::isfinite(phi)||!std::isfinite(theta)||!std::isfinite(psi))throw std::invalid_argument("nonfinite Euler angle");
  const T ca=std::cos(psi/T{2}),sa=std::sin(psi/T{2}),cb=std::cos(theta/T{2}),sb=std::sin(theta/T{2}),cc=std::cos(phi/T{2}),sc=std::sin(phi/T{2});
  return Su2<T>::from_unit_quaternion({ca*cb*cc+sa*sb*sc,ca*cb*sc-sa*sb*cc,ca*sb*cc+sa*cb*sc,sa*cb*cc-ca*sb*sc});
}

// Canonical representative: outer angles atan2 in [-pi,pi], third=0 at a
// singularity. At singularity individual outer angles cannot be recovered.
template<std::floating_point T> struct EulerAngles {
  T first{},middle{},last{};bool singular{};
};
template<std::floating_point T> EulerAngles<T> to_euler_zxz(const So3<T>& r,T tol=default_tolerance<T>()) {
  validate_tolerance(tol);const auto& m=r.matrix();const T s=std::hypot(m(0,2),m(1,2));
  if(s<=tol)return {std::atan2(m(1,0),m(0,0)),m(2,2)>=T{}?T{}:std::numbers::pi_v<T>,T{},true};
  return {std::atan2(m(0,2),-m(1,2)),std::atan2(s,m(2,2)),std::atan2(m(2,0),m(2,1)),false};
}
template<std::floating_point T> EulerAngles<T> to_euler_zyz(const So3<T>& r,T tol=default_tolerance<T>()) {
  validate_tolerance(tol);const auto& m=r.matrix();const T s=std::hypot(m(0,2),m(1,2));
  if(s<=tol){
    const T sign=m(2,2)>=T{}?T{1}:T{-1};
    return {std::atan2(sign*m(1,0),sign*m(0,0)),sign>T{}?T{}:std::numbers::pi_v<T>,T{},true};
  }
  return {std::atan2(m(1,2),m(0,2)),std::atan2(s,m(2,2)),std::atan2(m(2,1),-m(2,0)),false};
}
template<std::floating_point T> EulerAngles<T> to_yaw_pitch_roll(const So3<T>& r,T tol=default_tolerance<T>()) {
  validate_tolerance(tol);const auto& m=r.matrix();const T c=std::hypot(m(0,0),m(1,0));
  if(c<=tol)return {std::atan2(-m(0,1),m(1,1)),std::copysign(std::numbers::pi_v<T>/T{2},-m(2,0)),T{},true};
  return {std::atan2(m(1,0),m(0,0)),std::atan2(-m(2,0),c),std::atan2(m(2,1),m(2,2)),false};
}

// Largest quaternion-component pivot avoids division by sin(phi) near pi,
// and atan2 below avoids acos(trace) loss of small angles near identity.
// Normalize the extracted lift for roundoff; its sign is chosen w>=0.
// This does not change an independently imported Su2 lift or identify q with -q.
template<std::floating_point T> Su2<T> quaternion_from_rotation(const So3<T>& r) {
  const auto& m=r.matrix();const T tr=m.trace();
  const std::array<T,4> t{T{1}+tr,T{1}+T{2}*m(0,0)-tr,T{1}+T{2}*m(1,1)-tr,T{1}+T{2}*m(2,2)-tr};
  const auto p=std::size_t(std::max_element(t.begin(),t.end())-t.begin());
  std::array<T,4> q{};q[p]=std::sqrt(std::max(T{},t[p]))/T{2};const T d=T{4}*q[p];
  if(p==0){q[1]=(m(2,1)-m(1,2))/d;q[2]=(m(0,2)-m(2,0))/d;q[3]=(m(1,0)-m(0,1))/d;}
  else {
    const auto i=p-1,j=(i+1)%3,k=(i+2)%3;
    q[0]=(m(k,j)-m(j,k))/d;q[j+1]=(m(j,i)+m(i,j))/d;q[k+1]=(m(k,i)+m(i,k))/d;
  }
  const T length=std::hypot(std::hypot(q[0],q[1]),std::hypot(q[2],q[3]));
  for(auto& v:q)v/=length;
  if(q[0]<T{})for(auto& v:q)v=-v;
  return Su2<T>::from_unit_quaternion(q);
}
template<std::floating_point T> struct AxisAngle {
  Vector3<T> axis;T angle;bool axis_line_unique; // principal [0,pi]; at pi the oriented axis sign is ambiguous.
};
template<std::floating_point T> AxisAngle<T> to_axis_angle(const So3<T>& r) {
  const auto q=quaternion_from_rotation(r).quaternion();Vector3<T> v{q[1],q[2],q[3]};const T n=norm(v);
  if(n==T{})return {unit_axis<T>(0),T{},false};
  for(auto& x:v)x/=n;
  return {v,T{2}*std::atan2(n,q[0]),true};
}
template<std::floating_point T> So3Algebra<T> rotation_log(const So3<T>& r) {
  auto a=to_axis_angle(r);for(auto& x:a.axis)x*=a.angle;
  return So3Algebra<T>(a.axis);
}

// Goldstein's CK matrix K=U(q)^dagger; its SO(3) image is A_G=R(q)^T.
// alpha=w+iz, beta=y+ix, gamma=-y+ix=-beta*, delta=w-iz=alpha*.
template<class T> Matrix<Complex<T>,2,2> cayley_klein_matrix(const std::array<T,4>& q) {
  return adjoint(quaternion_matrix(q));
}
template<std::floating_point T> class CayleyKlein {
public:
  static CayleyKlein from_quaternion(Su2<T> q){return CayleyKlein(std::move(q));}
  static CayleyKlein from_parameters(Complex<T> alpha,Complex<T> beta,Complex<T> gamma,Complex<T> delta,T tol=default_tolerance<T>()) {
    const auto k=Matrix<Complex<T>,2,2>::from_rows({{{alpha,beta},{gamma,delta}}});
    return CayleyKlein(Su2<T>::from_matrix(adjoint(k),tol));
  }
  Matrix<Complex<T>,2,2> matrix()const{return cayley_klein_matrix(q_.quaternion());}
  const Su2<T>& quaternion()const{return q_;}
  So3<T> goldstein_rotation()const{return q_.to_so3().inverse();}
private:
  explicit CayleyKlein(Su2<T> q):q_(std::move(q)){}
  Su2<T> q_;
};

// Passive coordinate relabeling f'=Bf, g'=Bg => A'=BAB^T.
// A is arbitrary; only B is constrained to be an orthonormal coordinate change.
template<std::floating_point T> Matrix<T,3,3> coordinate_change(const So3<T>& b,const Matrix<T,3,3>& a) {
  for(const auto& row:a.rows())for(auto v:row)if(!std::isfinite(v))throw std::invalid_argument("nonfinite operator");
  const auto result=b.matrix()*a*b.matrix().transpose();
  for(const auto& row:result.rows())for(auto v:row)if(!std::isfinite(v))throw std::overflow_error("coordinate change overflow");
  return result;
}
template<std::floating_point T> So3Algebra<T> adjoint_action(const So3<T>& b,const So3Algebra<T>& x) {
  return So3Algebra<T>(b.apply(x.coordinates()));
}
// General 3x3 Cayley-Hamilton polynomial; no SO(3) assumption and no division.
template<class T> Matrix<T,3,3> cayley_hamilton_residual(const Matrix<T,3,3>& m) {
  const T c2=m(0,0)*m(1,1)-m(0,1)*m(1,0)+m(0,0)*m(2,2)-m(0,2)*m(2,0)+m(1,1)*m(2,2)-m(1,2)*m(2,1);
  const auto m2=m*m;
  return m2*m-m2.scaled(m.trace())+m.scaled(c2)-Matrix<T,3,3>::identity().scaled(m.determinant());
}
} // namespace mathphysics
