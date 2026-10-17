#pragma once
// Independently derived Frankel conventions. Zero-based, row-major, column vectors.
// General matrices and validated numerical groups are intentionally separate types.
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>

namespace mathphysics {
template<std::size_t N>
constexpr int levi_civita(const std::array<std::size_t,N>& indices) {
  for (auto i:indices) if(i>=N) throw std::out_of_range("epsilon index");
  int sign=1;
  for(std::size_t i=0;i<N;++i) for(std::size_t j=i+1;j<N;++j) {
    if(indices[i]==indices[j]) return 0;
    if(indices[i]>indices[j]) sign=-sign;
  }
  return sign;
}
constexpr int epsilon3(int i,int j,int k) {
  if(i<0 || j<0 || k<0) throw std::out_of_range("negative epsilon index");
  return levi_civita(std::array<std::size_t,3>{std::size_t(i),std::size_t(j),std::size_t(k)});
}

template<class T> struct Complex {
  T re{},im{};
  Complex()=default;
  explicit Complex(T r):re(std::move(r)){}
  Complex(T r,T i):re(std::move(r)),im(std::move(i)){}
  Complex conjugate() const {return {re,-im};}
  T norm_squared() const {return re*re+im*im;}
  bool operator==(const Complex&) const=default;
};
template<class T> Complex<T> operator+(const Complex<T>& a,const Complex<T>& b){return {a.re+b.re,a.im+b.im};}
template<class T> Complex<T> operator-(const Complex<T>& a){return {-a.re,-a.im};}
template<class T> Complex<T> operator-(const Complex<T>& a,const Complex<T>& b){return a+(-b);}
template<class T> Complex<T> operator*(const Complex<T>& a,const Complex<T>& b){return {a.re*b.re-a.im*b.im,a.re*b.im+a.im*b.re};}

template<class T,std::size_t R,std::size_t C> class Matrix {
public:
  using Rows=std::array<std::array<T,C>,R>;
  Matrix()=default;
  static Matrix from_rows(Rows entries){Matrix result;result.entries_=std::move(entries);return result;}
  static Matrix identity() requires(R==C) {
    Matrix result;for(std::size_t i=0;i<R;++i) result(i,i)=T{1};return result;
  }
  T& operator()(std::size_t i,std::size_t j){return entries_.at(i).at(j);}
  const T& operator()(std::size_t i,std::size_t j)const{return entries_.at(i).at(j);}
  const Rows& rows()const{return entries_;}
  Matrix<T,C,R> transpose()const {
    Matrix<T,C,R> r;for(std::size_t i=0;i<R;++i) for(std::size_t j=0;j<C;++j) r(j,i)=(*this)(i,j);return r;
  }
  Matrix scaled(const T& a)const {
    Matrix r;for(std::size_t i=0;i<R;++i) for(std::size_t j=0;j<C;++j) r(i,j)=(*this)(i,j)*a;return r;
  }
  std::array<T,R> apply(const std::array<T,C>& v)const {
    std::array<T,R> r{};for(std::size_t i=0;i<R;++i) for(std::size_t j=0;j<C;++j) r[i]=r[i]+(*this)(i,j)*v[j];return r;
  }
  T trace()const requires(R==C) {T s{};for(std::size_t i=0;i<R;++i)s=s+(*this)(i,i);return s;}
  T determinant()const requires(R==2 && C==2){return (*this)(0,0)*(*this)(1,1)-(*this)(0,1)*(*this)(1,0);}
  T determinant()const requires(R==3 && C==3){
    return (*this)(0,0)*((*this)(1,1)*(*this)(2,2)-(*this)(1,2)*(*this)(2,1))
      -(*this)(0,1)*((*this)(1,0)*(*this)(2,2)-(*this)(1,2)*(*this)(2,0))
      +(*this)(0,2)*((*this)(1,0)*(*this)(2,1)-(*this)(1,1)*(*this)(2,0));
  }
  bool operator==(const Matrix&)const=default;
private:
  Rows entries_{};
};
template<class T,std::size_t R,std::size_t C>
Matrix<T,R,C> operator+(const Matrix<T,R,C>& a,const Matrix<T,R,C>& b) {
  Matrix<T,R,C> r;for(std::size_t i=0;i<R;++i)for(std::size_t j=0;j<C;++j)r(i,j)=a(i,j)+b(i,j);return r;
}
template<class T,std::size_t R,std::size_t C>
Matrix<T,R,C> operator-(const Matrix<T,R,C>& a,const Matrix<T,R,C>& b) {
  Matrix<T,R,C> r;for(std::size_t i=0;i<R;++i)for(std::size_t j=0;j<C;++j)r(i,j)=a(i,j)-b(i,j);return r;
}
template<class T,std::size_t R,std::size_t C,std::size_t K>
Matrix<T,R,K> operator*(const Matrix<T,R,C>& a,const Matrix<T,C,K>& b) {
  Matrix<T,R,K> r;
  for(std::size_t i=0;i<R;++i)for(std::size_t j=0;j<K;++j)for(std::size_t k=0;k<C;++k)r(i,j)=r(i,j)+a(i,k)*b(k,j);
  return r;
}
template<class T,std::size_t N>
Matrix<T,N,N> commutator(const Matrix<T,N,N>& a,const Matrix<T,N,N>& b){return a*b-b*a;}
template<class T,std::size_t R,std::size_t C>
Matrix<Complex<T>,C,R> adjoint(const Matrix<Complex<T>,R,C>& a) {
  Matrix<Complex<T>,C,R> r;for(std::size_t i=0;i<R;++i)for(std::size_t j=0;j<C;++j)r(j,i)=a(i,j).conjugate();return r;
}
template<class T> using Vector3=std::array<T,3>;
template<class T> T dot(const Vector3<T>& a,const Vector3<T>& b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
template<class T> Vector3<T> cross(const Vector3<T>& a,const Vector3<T>& b) {
  return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};
}
template<class T> Matrix<T,3,3> hat(const Vector3<T>& v) {
  return Matrix<T,3,3>::from_rows({{{T{},-v[2],v[1]},{v[2],T{},-v[0]},{-v[1],v[0],T{}}}});
}
template<class T> Vector3<T> unit_axis(std::size_t index) {
  if(index>=3)throw std::out_of_range("basis index");
  Vector3<T> v{};v[index]=T{1};return v;
}
template<class T> Matrix<Complex<T>,2,2> pauli(std::size_t index) {
  using Z=Complex<T>;const Z z{},o{T{1}},i{T{},T{1}};
  switch(index){
    case 0:return Matrix<Z,2,2>::from_rows({{{z,o},{o,z}}});
    case 1:return Matrix<Z,2,2>::from_rows({{{z,-i},{i,z}}});
    case 2:return Matrix<Z,2,2>::from_rows({{{o,z},{z,-o}}});
    default:throw std::out_of_range("Pauli index");
  }
}
enum class Su2Basis{HalfPauliOverI,PauliOverI};
template<class T> Matrix<Complex<T>,2,2> generator(std::size_t index,Su2Basis basis) {
  const T coefficient=basis==Su2Basis::HalfPauliOverI? T{-1}/T{2}:T{-1};
  return pauli<T>(index).scaled(Complex<T>{T{},coefficient});
}
template<class T> std::array<T,4> hamilton(const std::array<T,4>& a,const std::array<T,4>& b) {
  const Vector3<T> av{a[1],a[2],a[3]},bv{b[1],b[2],b[3]};const auto c=cross(av,bv);
  return {a[0]*b[0]-dot(av,bv),a[0]*b[1]+b[0]*a[1]+c[0],a[0]*b[2]+b[0]*a[2]+c[1],a[0]*b[3]+b[0]*a[3]+c[2]};
}
template<class T> Matrix<Complex<T>,2,2> quaternion_matrix(const std::array<T,4>& q) {
  using Z=Complex<T>;
  return Matrix<Z,2,2>::from_rows({{{Z{q[0],-q[3]},Z{-q[2],-q[1]}},{Z{q[2],-q[1]},Z{q[0],q[3]}}}});
}
template<class T> Matrix<T,3,3> quaternion_rotation(const std::array<T,4>& q) {
  const Vector3<T> v{q[1],q[2],q[3]};
  auto r=Matrix<T,3,3>::identity().scaled(T{q[0]*q[0]-dot(v,v)})+hat(v).scaled(T{2}*q[0]);
  for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j)r(i,j)=r(i,j)+T{2}*v[i]*v[j];
  return r;
}
template<std::floating_point T> T default_tolerance(){return T{32}*std::numeric_limits<T>::epsilon();}
template<std::floating_point T> void validate_tolerance(T t) {
  if(!std::isfinite(t)||t<T{}||t>=T{0.125})throw std::invalid_argument("tolerance outside [0,1/8)");
}
template<std::floating_point T> T norm(const Vector3<T>& v) {
  for(auto x:v)if(!std::isfinite(x))throw std::invalid_argument("nonfinite vector");
  const auto r=std::hypot(v[0],v[1],v[2]);
  if(!std::isfinite(r))throw std::invalid_argument("nonfinite norm");
  return r;
}
template<std::floating_point T> class So3;
template<std::floating_point T> class Su2;

template<class T> class So3Algebra {
public:
  explicit So3Algebra(Vector3<T> v):v_(std::move(v)){}
  static So3Algebra basis(std::size_t i){return So3Algebra(unit_axis<T>(i));}
  const Vector3<T>& coordinates()const{return v_;}
  Matrix<T,3,3> matrix()const{return hat(v_);}
  So3Algebra bracket(const So3Algebra& b)const{return So3Algebra(cross(v_,b.v_));}
  static So3Algebra from_matrix(const Matrix<T,3,3>& m) {
    if(m.transpose()!=m.scaled(T{-1}))throw std::invalid_argument("matrix not skew");
    for(std::size_t i=0;i<3;++i)if(m(i,i)!=T{})throw std::invalid_argument("nonzero skew diagonal");
    return So3Algebra({m(2,1),m(0,2),m(1,0)});
  }
  auto exp()const requires std::floating_point<T>;
private:
  Vector3<T> v_;
};
template<std::floating_point T> class So3 {
public:
  static So3 identity(){return So3(Matrix<T,3,3>::identity(),Valid{});}
  const Matrix<T,3,3>& matrix()const{return matrix_;}
  static So3 from_matrix(Matrix<T,3,3> m,T tol=default_tolerance<T>()){
    validate_tolerance(tol);
    for(const auto& row:m.rows())for(auto x:row)if(!std::isfinite(x))throw std::invalid_argument("nonfinite rotation");
    const auto residual=m.transpose()*m-Matrix<T,3,3>::identity();const T det=m.determinant();
    for(const auto& row:residual.rows())for(auto x:row)if(!std::isfinite(x)||std::abs(x)>tol)throw std::invalid_argument("not orthogonal");
    if(!std::isfinite(det)||std::abs(det-T{1})>tol)throw std::invalid_argument("determinant not +1");
    return So3(std::move(m),Valid{});
  }
  static So3 from_axis_angle(Vector3<T> axis,T angle){
    if(!std::isfinite(angle))throw std::invalid_argument("nonfinite angle");
    const T n=norm(axis);if(n==T{})throw std::invalid_argument("zero axis");
    for(auto& x:axis)x=(x/n)*angle;
    return So3Algebra<T>(axis).exp();
  }
  So3 compose(const So3& b)const{return from_matrix(matrix_*b.matrix_);}
  So3 inverse()const{return from_matrix(matrix_.transpose());}
  Vector3<T> apply(const Vector3<T>& v)const{
    for(auto x:v)if(!std::isfinite(x))throw std::invalid_argument("nonfinite vector");
    const auto r=matrix_.apply(v);
    for(auto x:r)if(!std::isfinite(x))throw std::invalid_argument("nonfinite rotated vector");
    return r;
  }
private:
  struct Valid{};
  So3(Matrix<T,3,3> m,Valid):matrix_(std::move(m)){}
  Matrix<T,3,3> matrix_;
};
template<class T> auto So3Algebra<T>::exp()const requires std::floating_point<T>{
  const T r=norm(v_);Matrix<T,3,3> m;
  if(std::abs(r)<T{0.125}){
    const T x=r*r;
    const T a=T{1}+x*(T{-1}/T{6}+x*(T{1}/T{120}+x*(T{-1}/T{5040}+x/T{362880})));
    const T b=T{1}/T{2}+x*(T{-1}/T{24}+x*(T{1}/T{720}+x*(T{-1}/T{40320}+x/T{3628800})));
    const auto k=matrix();m=Matrix<T,3,3>::identity()+k.scaled(a)+(k*k).scaled(b);
  }else{
    auto n=v_;for(auto& x:n)x=x/r;const auto k=hat(n);
    m=Matrix<T,3,3>::identity()+k.scaled(std::sin(r))+(k*k).scaled(T{1}-std::cos(r));
  }
  return So3<T>::from_matrix(m);
}
template<class T> class Su2Algebra {
public:
  explicit Su2Algebra(Vector3<T> v):v_(std::move(v)){}
  static Su2Algebra basis(std::size_t i){return Su2Algebra(unit_axis<T>(i));}
  static Su2Algebra from_matrix(const Matrix<Complex<T>,2,2>& m){
    if(adjoint(m)!=m.scaled(Complex<T>{T{-1}})||m.trace()!=Complex<T>{})throw std::invalid_argument("not traceless anti-Hermitian");
    const Su2Algebra result({-(m(0,1).im+m(1,0).im),m(1,0).re-m(0,1).re,m(1,1).im-m(0,0).im});
    if(result.matrix()!=m)throw std::invalid_argument("not real su2 coefficients");
    return result;
  }
  static Su2Algebra from_full_basis_coordinates(Vector3<T> v){for(auto& x:v)x=T{2}*x;return Su2Algebra(v);}
  const Vector3<T>& coordinates()const{return v_;}
  Matrix<Complex<T>,2,2> matrix()const{
    Matrix<Complex<T>,2,2> r;
    for(std::size_t i=0;i<3;++i)r=r+generator<T>(i,Su2Basis::HalfPauliOverI).scaled(Complex<T>{v_[i]});
    return r;
  }
  Su2Algebra bracket(const Su2Algebra& b)const{return Su2Algebra(cross(v_,b.v_));}
  So3Algebra<T> to_so3()const{return So3Algebra<T>(v_);}
  auto exp()const requires std::floating_point<T>;
private:
  Vector3<T> v_;
};
template<std::floating_point T> class Su2 {
public:
  static Su2 identity(){return Su2({T{1},T{},T{},T{}},Valid{});}
  const std::array<T,4>& quaternion()const{return q_;}
  static Su2 from_unit_quaternion(std::array<T,4> q,T tol=default_tolerance<T>()){
    validate_tolerance(tol);T sum{};
    for(auto x:q){if(!std::isfinite(x))throw std::invalid_argument("nonfinite quaternion");sum+=x*x;}
    if(!std::isfinite(sum)||std::abs(sum-T{1})>tol)throw std::invalid_argument("quaternion not unit");
    return Su2(q,Valid{});
  }
  static Su2 from_matrix(const Matrix<Complex<T>,2,2>& m,T tol=default_tolerance<T>()){
    validate_tolerance(tol);
    for(const auto& row:m.rows())for(auto x:row)if(!std::isfinite(x.re)||!std::isfinite(x.im))throw std::invalid_argument("nonfinite matrix");
    const auto residual=adjoint(m)*m-Matrix<Complex<T>,2,2>::identity();
    const auto det=m.determinant()-Complex<T>{T{1}};
    const auto good=[tol](const Complex<T>& z){return std::isfinite(z.re)&&std::isfinite(z.im)&&std::abs(z.re)<=tol&&std::abs(z.im)<=tol;};
    for(const auto& row:residual.rows())for(auto z:row)if(!good(z))throw std::invalid_argument("not unitary");
    if(!good(det))throw std::invalid_argument("not special unitary");
    const T h=T{1}/T{2};
    auto result=from_unit_quaternion({(m(0,0).re+m(1,1).re)*h,-(m(0,1).im+m(1,0).im)*h,(m(1,0).re-m(0,1).re)*h,(m(1,1).im-m(0,0).im)*h},tol);
    const auto difference=m-result.matrix();
    for(const auto& row:difference.rows())for(auto z:row)if(!good(z))throw std::invalid_argument("matrix outside quaternion form");
    return result;
  }
  Matrix<Complex<T>,2,2> matrix()const{return quaternion_matrix(q_);}
  So3<T> to_so3()const{return So3<T>::from_matrix(quaternion_rotation(q_));}
  Su2 inverse()const{return Su2({q_[0],-q_[1],-q_[2],-q_[3]},Valid{});}
  Su2 antipode()const{return Su2({-q_[0],-q_[1],-q_[2],-q_[3]},Valid{});}
  Su2 compose(const Su2& b)const{return from_unit_quaternion(hamilton(q_,b.q_));}
  static Su2 exp_pauli_over_i(const Vector3<T>& v){
    const T r=norm(v);std::array<T,4> q{std::cos(r),T{},T{},T{}};
    if(r!=T{})for(std::size_t i=0;i<3;++i)q[i+1]=(v[i]/r)*std::sin(r);
    return from_unit_quaternion(q);
  }
private:
  struct Valid{};
  Su2(std::array<T,4> q,Valid):q_(q){}
  std::array<T,4> q_;
};
template<class T> auto Su2Algebra<T>::exp()const requires std::floating_point<T>{
  const T r=norm(v_),h=r/T{2};std::array<T,4> q{std::cos(h),T{},T{},T{}};
  if(std::abs(r)<T{0.125}){
    const T x=h*h;
    const T sinc=T{1}+x*(T{-1}/T{6}+x*(T{1}/T{120}-x/T{5040}));
    for(std::size_t i=0;i<3;++i)q[i+1]=(v_[i]/T{2})*sinc;
  }else for(std::size_t i=0;i<3;++i)q[i+1]=(v_[i]/r)*std::sin(h);
  return Su2<T>::from_unit_quaternion(q);
}
} // namespace mathphysics
