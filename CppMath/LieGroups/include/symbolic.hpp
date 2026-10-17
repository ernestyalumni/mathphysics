#pragma once
#include "lie_groups.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <map>
#include <ostream>
#include <string>
#include <vector>

namespace mathphysics {
using Exact=boost::multiprecision::cpp_rational;
// Independent sparse multivariate polynomial over exact rationals.
// Coefficients are arbitrary precision; degree/term work bounds remain explicit.
class Polynomial {
public:
  Polynomial()=default;
  Polynomial(int n):Polynomial(Exact{n}){}
  explicit Polynomial(Exact n){insert({},std::move(n));}
  static Polynomial var(std::string name){Polynomial p;p.insert({std::move(name)},Exact{1});return p;}
  bool is_zero()const{return terms_.empty();}
  Exact constant_value()const{
    if(terms_.empty())return Exact{};
    if(terms_.size()!=1||!terms_.begin()->first.empty())throw std::invalid_argument("not constant");
    return terms_.begin()->second;
  }
  Exact evaluate(const std::map<std::string,Exact>& values)const{
    Exact result{};
    for(const auto& [monomial,c]:terms_){Exact term=c;for(const auto& v:monomial)term*=values.at(v);result+=term;}
    return result;
  }
  bool operator==(const Polynomial&)const=default;
  friend Polynomial operator+(Polynomial a,const Polynomial& b){for(const auto& [m,c]:b.terms_)a.insert(m,c);return a;}
  friend Polynomial operator-(Polynomial a){for(auto& [m,c]:a.terms_){(void)m;c=-c;}return a;}
  friend Polynomial operator-(const Polynomial& a,const Polynomial& b){return a+(-b);}
  friend Polynomial operator*(const Polynomial& a,const Polynomial& b){
    Polynomial r;
    for(const auto& [m,c]:a.terms_)for(const auto& [n,d]:b.terms_){
      auto key=m;key.insert(key.end(),n.begin(),n.end());std::sort(key.begin(),key.end());r.insert(std::move(key),Exact{c*d});
    }return r;
  }
  // Only constant divisors; NOT a general polynomial field operation.
  friend Polynomial operator/(Polynomial a,const Polynomial& b){
    const Exact d=b.constant_value();if(d==0)throw std::domain_error("zero polynomial divisor");
    for(auto& [m,c]:a.terms_){(void)m;c/=d;}return a;
  }
  friend std::ostream& operator<<(std::ostream& os,const Polynomial& p){
    if(p.is_zero())return os<<'0';
    bool first=true;for(const auto& [m,c]:p.terms_){if(!first)os<<" + ";first=false;os<<c;for(const auto& v:m)os<<'*'<<v;}return os;
  }
private:
  using Monomial=std::vector<std::string>;
  std::map<Monomial,Exact> terms_;
  void insert(Monomial m,Exact c){
    if(m.size()>64)throw std::length_error("polynomial degree budget");
    auto found=terms_.find(m);
    if(found!=terms_.end()){c+=found->second;terms_.erase(found);}
    if(c!=0)terms_.emplace(std::move(m),std::move(c));
    if(terms_.size()>8192)throw std::length_error("polynomial term budget");
  }
};
// Exact Q[s,c]/(s^2+c^2-1); one common formal angle only.
class Trig {
public:
  Trig()=default;
  Trig(int n):a_(n){}
  explicit Trig(Exact n):a_(std::move(n)){}
  static Trig sine(){return Trig(Polynomial::var("s"),Polynomial{});}
  static Trig cosine(){return Trig(Polynomial{},Polynomial{1});}
  Exact at(const Exact& s,const Exact& c)const{
    if(s*s+c*c!=1)throw std::invalid_argument("not a unit circle");
    const std::map<std::string,Exact> v{{"s",s}};return a_.evaluate(v)+c*b_.evaluate(v);
  }
  bool operator==(const Trig&)const=default;
  friend Trig operator+(const Trig& a,const Trig& b){return Trig(a.a_+b.a_,a.b_+b.b_);}
  friend Trig operator-(const Trig& a){return Trig(-a.a_,-a.b_);}
  friend Trig operator-(const Trig& a,const Trig& b){return a+(-b);}
  friend Trig operator*(const Trig& a,const Trig& b){
    const auto s=Polynomial::var("s");
    return Trig(a.a_*b.a_+(Polynomial{1}-s*s)*a.b_*b.b_,a.a_*b.b_+a.b_*b.a_);
  }
  friend Trig operator/(const Trig& a,const Trig& b){
    if(!b.b_.is_zero())throw std::invalid_argument("nonconstant trig divisor");
    return Trig(a.a_/b.a_,a.b_/b.a_);
  }
  friend std::ostream& operator<<(std::ostream& os,const Trig& v){return os<<'('<<v.a_<<") + c*("<<v.b_<<')';}
private:
  Trig(Polynomial a,Polynomial b):a_(std::move(a)),b_(std::move(b)){}
  Polynomial a_,b_;
};
inline void validate_unit_axis(const Vector3<Exact>& axis){if(dot(axis,axis)!=1)throw std::invalid_argument("axis not exactly unit");}
inline Matrix<Trig,3,3> so3_exponential(const Vector3<Exact>& axis){
  validate_unit_axis(axis);Vector3<Trig> n{Trig{axis[0]},Trig{axis[1]},Trig{axis[2]}};const auto k=hat(n);
  return Matrix<Trig,3,3>::identity()+k.scaled(Trig::sine())+(k*k).scaled(Trig{1}-Trig::cosine());
}
struct Exponential2 {Matrix<Complex<Trig>,2,2> matrix;Exact argument_scale;};
inline Exponential2 su2_exponential(const Vector3<Exact>& axis,Su2Basis basis){
  validate_unit_axis(axis);const auto s=Trig::sine();
  return {quaternion_matrix(std::array<Trig,4>{Trig::cosine(),Trig{axis[0]}*s,Trig{axis[1]}*s,Trig{axis[2]}*s}),
    basis==Su2Basis::HalfPauliOverI?Exact{1}/2:Exact{1}};
}
// Read-only denominator invariant: never manufacture a false 0/0 equality.
template<class T,std::size_t N> class RationalMatrix {
public:
  RationalMatrix(Matrix<T,N,N> numerator,Polynomial denominator)
    :numerator_(std::move(numerator)),denominator_(std::move(denominator)){
    if(denominator_.is_zero())throw std::domain_error("identically zero denominator");
  }
  RationalMatrix(const RationalMatrix&)=default;
  // Moved-from values remain defined zero matrices with denominator one.
  RationalMatrix(RationalMatrix&& other):numerator_{},denominator_{1}{
    using std::swap;swap(numerator_,other.numerator_);swap(denominator_,other.denominator_);
  }
  RationalMatrix& operator=(RationalMatrix other){
    using std::swap;swap(numerator_,other.numerator_);swap(denominator_,other.denominator_);return *this;
  }
  const Matrix<T,N,N>& numerator()const{return numerator_;}
  const Polynomial& denominator()const{return denominator_;}
  bool equals(const RationalMatrix& b)const{return numerator_.scaled(T{b.denominator_})==b.numerator_.scaled(T{denominator_});}
private:
  Matrix<T,N,N> numerator_;Polynomial denominator_;
};
enum class Chart{Plus,Minus};
class SymbolicSu2 {
public:
  SymbolicSu2(const SymbolicSu2&)=default;
  // Leave a valid identity behind instead of a 0/0 quaternion after moves.
  SymbolicSu2(SymbolicSu2&& other):q_{Polynomial{1},Polynomial{},Polynomial{},Polynomial{}},d_{1}{
    using std::swap;q_.swap(other.q_);swap(d_,other.d_);
  }
  SymbolicSu2& operator=(SymbolicSu2 other){
    using std::swap;q_.swap(other.q_);swap(d_,other.d_);return *this;
  }
  static SymbolicSu2 stereographic(const Vector3<Polynomial>& t,Chart chart){
    const auto square=dot(t,t);const Polynomial sign{chart==Chart::Plus?1:-1};
    return SymbolicSu2({(Polynomial{1}-square)*sign,Polynomial{2}*t[0]*sign,Polynomial{2}*t[1]*sign,Polynomial{2}*t[2]*sign},Polynomial{1}+square);
  }
  static SymbolicSu2 identity(){return SymbolicSu2({Polynomial{1},Polynomial{},Polynomial{},Polynomial{}},Polynomial{1});}
  Polynomial norm_residual()const{Polynomial r;for(const auto& v:q_)r=r+v*v;return r-d_*d_;}
  RationalMatrix<Complex<Polynomial>,2> matrix()const{return {quaternion_matrix(q_),d_};}
  RationalMatrix<Polynomial,3> so3_matrix()const{return {quaternion_rotation(q_),d_*d_};}
  SymbolicSu2 inverse()const{return SymbolicSu2({q_[0],-q_[1],-q_[2],-q_[3]},d_);}
  SymbolicSu2 antipode()const{return SymbolicSu2({-q_[0],-q_[1],-q_[2],-q_[3]},d_);}
  SymbolicSu2 compose(const SymbolicSu2& b)const{return SymbolicSu2(hamilton(q_,b.q_),d_*b.d_);}
  bool equals(const SymbolicSu2& b)const{
    for(std::size_t i=0;i<4;++i)if(q_[i]*b.d_!=b.q_[i]*d_)return false;
    return true;
  }
private:
  SymbolicSu2(std::array<Polynomial,4> q,Polynomial d):q_(std::move(q)),d_(std::move(d)){}
  std::array<Polynomial,4> q_;Polynomial d_;
};
class SymbolicSo3 {
public:
  static SymbolicSo3 stereographic(const Vector3<Polynomial>& t){return SymbolicSo3(SymbolicSu2::stereographic(t,Chart::Plus));}
  explicit SymbolicSo3(SymbolicSu2 lift):lift_(std::move(lift)){}
  RationalMatrix<Polynomial,3> matrix()const{return lift_.so3_matrix();}
  SymbolicSo3 inverse()const{return SymbolicSo3(lift_.inverse());}
  SymbolicSo3 compose(const SymbolicSo3& b)const{return SymbolicSo3(lift_.compose(b.lift_));}
  bool equals(const SymbolicSo3& b)const{return matrix().equals(b.matrix());}
private:
  SymbolicSu2 lift_;
};
} // namespace mathphysics
