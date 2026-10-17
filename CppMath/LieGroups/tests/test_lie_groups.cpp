#include "lie_groups.hpp"
#include "symbolic.hpp"
#include <functional>
#include <iostream>
#include <numbers>
#include <string>
#include <type_traits>
using namespace mathphysics;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
template<class F>void rejects(F f){bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}require(rejected,"invalid input accepted");}
template<class T>void close(T a,long double b,long double tol){require(std::abs(static_cast<long double>(a)-b)<=tol,"numeric mismatch");}
void epsilon_contract(){
  static_assert(epsilon3(0,1,2)==1);
  for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int l=0;l<3;++l)for(int m=0;m<3;++m){
    int s=0;for(int k=0;k<3;++k)s+=epsilon3(i,j,k)*epsilon3(k,l,m);
    require(s==int(i==l)*int(j==m)-int(i==m)*int(j==l),"epsilon contraction sign");
  }
  rejects([]{epsilon3(0,1,3);});rejects([]{epsilon3(-1,0,0);});
}
void matrices(){
  const auto a=Matrix<int,3,3>::from_rows({{{12,8,4},{3,17,14},{9,8,10}}});
  const auto b=Matrix<int,3,3>::from_rows({{{5,19,3},{6,15,9},{7,8,16}}});
  require(a*b==Matrix<int,3,3>::from_rows({{{136,380,172},{215,424,386},{163,371,259}}}),"legacy product");
  const auto r=Matrix<int,2,3>::from_rows({{{1,2,3},{4,5,6}}});
  require(r*r.transpose()==Matrix<int,2,2>::from_rows({{{14,32},{32,77}}}),"rectangular product");
  static_assert(!std::is_constructible_v<So3<double>,Matrix<double,3,3>>);
  rejects([&]{(void)a(3,0);});
}
void exact_bases(){
  using T=Exact;
  for(std::size_t i=0;i<3;++i){
    const auto a=So3Algebra<T>::basis(i).matrix();
    require(a.transpose()==a.scaled(T{-1}),"skew basis");
    const auto s=pauli<T>(i),t=generator<T>(i,Su2Basis::HalfPauliOverI),f=generator<T>(i,Su2Basis::PauliOverI);
    require(adjoint(s)==s,"Pauli Hermitian");require(s*s==Matrix<Complex<T>,2,2>::identity(),"Pauli square");
    require(s.determinant()==Complex<T>{T{-1}},"Pauli not SU2");
    require(adjoint(t)==t.scaled(Complex<T>{T{-1}}),"T anti-Hermitian");
    require(Su2Algebra<T>::from_matrix(t).matrix()==t,"su2 algebra conversion");
    rejects([&]{Su2Algebra<T>::from_matrix(s);});
    require(f==t.scaled(Complex<T>{T{2}}),"basis scaling");
    for(std::size_t j=0;j<3;++j){
      const auto b=So3Algebra<T>::basis(j).matrix();
      Matrix<T,3,3> target;
      Matrix<Complex<T>,2,2> pauli_target,t_target;
      for(std::size_t k=0;k<3;++k){
        const auto e=levi_civita(std::array<std::size_t,3>{i,j,k});
        target=target+So3Algebra<T>::basis(k).matrix().scaled(T{e});
        pauli_target=pauli_target+pauli<T>(k).scaled(Complex<T>{T{},T{2*e}});
        t_target=t_target+generator<T>(k,Su2Basis::HalfPauliOverI).scaled(Complex<T>{T{e}});
      }
      require(commutator(a,b)==target,"so3 structure constants");
      require((a*b).trace()==(i==j?T{-2}:T{}),"trace normalization");
      require(commutator(s,pauli<T>(j))==pauli_target,"Pauli commutator");
      require(commutator(t,generator<T>(j,Su2Basis::HalfPauliOverI))==t_target,"T commutator");
      require(commutator(f,generator<T>(j,Su2Basis::PauliOverI))==t_target.scaled(Complex<T>{T{4}}),"F factor two");
    }
  }
  const auto e=So3Algebra<T>::basis(2);
  require(e.matrix().apply(unit_axis<T>(0))==unit_axis<T>(1),"hat orientation");
  rejects([]{So3Algebra<T>::basis(3);});
  rejects([]{So3Algebra<T>::from_matrix(Matrix<T,3,3>::identity());});
}
template<class T>void numeric_groups(){
  const long double tol=static_cast<long double>(default_tolerance<T>())*4;
  const T angle=std::numbers::pi_v<T>/T{2};
  const auto z=So3<T>::from_axis_angle({T{},T{},T{1}},angle);
  const auto image=z.apply({T{1},T{},T{}});
  close(image[0],0,tol);close(image[1],1,tol);close(image[2],0,tol);
  for(int sample=0;sample<32;++sample){
    const Vector3<T> v{T(0.11*(sample-10)),T(0.07*(sample+1)),T(-0.03*(sample+3))};
    const auto r=So3Algebra<T>(v).exp();const auto u=Su2Algebra<T>(v).exp();
    const auto q=u.to_so3(),inverse=r.inverse().compose(r);
    for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j){
      close(q.matrix()(i,j),r.matrix()(i,j),tol);
      close(inverse.matrix()(i,j),i==j?1:0,tol);
    }
    require(q.matrix()==u.antipode().to_so3().matrix(),"antipodal covering");
    const auto round=Su2<T>::from_matrix(u.matrix());
    for(std::size_t i=0;i<4;++i)close(round.quaternion()[i],u.quaternion()[i],tol);
  }
  const auto negative=Su2Algebra<T>({T{},T{},T{2}*std::numbers::pi_v<T>}).exp();
  close(negative.quaternion()[0],-1,tol);
  require(negative.quaternion()[0]<T{},"spinor 2pi is -I");
  const auto full=Su2<T>::exp_pauli_over_i({T{},T{},angle});
  close(full.quaternion()[3],1,tol);
  const auto a=z.compose(So3<T>::from_axis_angle({T{},T{1},T{}},angle)).apply({T{},T{},T{1}});
  close(a[1],1,tol);
  const auto small=So3Algebra<T>({T{},T{},T(1e-6)}).exp();
  close(small.matrix()(1,0),1e-6,tol);
}
void invalid_membership(){
  rejects([]{So3<double>::from_matrix(Matrix<double,3,3>{});});
  auto reflection=Matrix<double,3,3>::identity();reflection(0,0)=-1;
  rejects([&]{So3<double>::from_matrix(reflection);});
  auto bad=Matrix<double,3,3>::identity();bad(0,0)=std::numeric_limits<double>::quiet_NaN();
  rejects([&]{So3<double>::from_matrix(bad);});
  rejects([]{So3<double>::from_matrix(Matrix<double,3,3>::identity(),1.0);});
  rejects([]{So3<double>::from_axis_angle({0.,0.,0.},0.);});
  rejects([]{Su2<double>::from_unit_quaternion({0.,0.,0.,0.});});
  rejects([]{Su2<double>::from_matrix(pauli<double>(2));});
  rejects([]{So3Algebra<double>({std::numeric_limits<double>::infinity(),0.,0.}).exp();});
}
void abstract_polynomials(){
  const Vector3<Polynomial> a{Polynomial::var("a0"),Polynomial::var("a1"),Polynomial::var("a2")};
  const Vector3<Polynomial> b{Polynomial::var("b0"),Polynomial::var("b1"),Polynomial::var("b2")};
  require(hat(a).apply(b)==cross(a,b),"symbolic hat cross");
  require(commutator(hat(a),hat(b))==hat(cross(a,b)),"symbolic bracket");
  const Su2Algebra<Polynomial> x(a),y(b);
  require(commutator(x.matrix(),y.matrix())==x.bracket(y).matrix(),"symbolic su2 bracket");
}
void symbolic_exponentials(){
  const Vector3<Exact> n{Exact{3}/5,Exact{4}/5,Exact{}};
  const auto r=so3_exponential(n);
  require(r.transpose()*r==Matrix<Trig,3,3>::identity(),"symbolic Rodrigues orthogonal");
  require(r.determinant()==Trig{1},"symbolic determinant");
  const auto full=su2_exponential(n,Su2Basis::PauliOverI),half=su2_exponential(n,Su2Basis::HalfPauliOverI);
  require(full.argument_scale==1&&half.argument_scale==Exact{1}/2,"phase metadata");
  require(adjoint(full.matrix)*full.matrix==Matrix<Complex<Trig>,2,2>::identity(),"symbolic SU2 unitary");
  require(full.matrix.determinant()==Complex<Trig>{Trig{1}},"symbolic SU2 determinant");
  require(Trig::sine()*Trig::sine()+Trig::cosine()*Trig::cosine()==Trig{1},"circle quotient");
  require(Trig::sine().at(Exact{3}/5,Exact{4}/5)==Exact{3}/5,"circle evaluation");
  rejects([]{Trig::sine().at(Exact{1},Exact{1});});
  rejects([]{so3_exponential(Vector3<Exact>{Exact{1},Exact{1},Exact{1}});});
}
void symbolic_charts(){
  const Vector3<Polynomial> t{Polynomial::var("t0"),Polynomial::var("t1"),Polynomial::var("t2")};
  const auto u=SymbolicSu2::stereographic(t,Chart::Plus);
  require(u.norm_residual().is_zero(),"chart norm");
  rejects([]{RationalMatrix<Polynomial,2> invalid(Matrix<Polynomial,2,2>{},Polynomial{});});
  const auto m=u.matrix();const auto r=u.so3_matrix();
  require(adjoint(m.numerator())*m.numerator()==Matrix<Complex<Polynomial>,2,2>::identity().scaled(Complex<Polynomial>{m.denominator()*m.denominator()}),"chart unitarity");
  require(m.numerator().determinant()==Complex<Polynomial>{m.denominator()*m.denominator()},"chart det SU2");
  require(r.numerator().transpose()*r.numerator()==Matrix<Polynomial,3,3>::identity().scaled(r.denominator()*r.denominator()),"chart orthogonality");
  require(r.numerator().determinant()==r.denominator()*r.denominator()*r.denominator(),"chart det SO3");
  require(u.compose(u.inverse()).equals(SymbolicSu2::identity()),"chart inverse");
  require(!u.equals(u.antipode()),"opposite SU2 values");
  require(u.so3_matrix().equals(u.antipode().so3_matrix()),"same SO3 values");
  require(SymbolicSu2::stereographic({Polynomial{},Polynomial{},Polynomial{}},Chart::Minus).equals(SymbolicSu2::identity().antipode()),"second chart");
}
void polynomial_homomorphism(){
  const std::array<Polynomial,4> a{Polynomial::var("p0"),Polynomial::var("p1"),Polynomial::var("p2"),Polynomial::var("p3")};
  const std::array<Polynomial,4> b{Polynomial::var("q0"),Polynomial::var("q1"),Polynomial::var("q2"),Polynomial::var("q3")};
  require(quaternion_matrix(hamilton(a,b))==quaternion_matrix(a)*quaternion_matrix(b),"matrix Hamilton homomorphism");
  require(quaternion_rotation(hamilton(a,b))==quaternion_rotation(a)*quaternion_rotation(b),"cover homomorphism");
}
void checked_inverse_and_move_invariants(){
  const auto m=Matrix<double,3,3>::from_rows({{{0.625,0.625,0.625},{0.6875,-0.6875,0.0},{0.375,0.375,-0.75}}});
  const auto admitted=So3<double>::from_matrix(m,0.12);
  rejects([&]{admitted.inverse();});
  const auto identity=RationalMatrix<Polynomial,2>(Matrix<Polynomial,2,2>::identity(),Polynomial{1});
  auto source=RationalMatrix<Polynomial,2>(Matrix<Polynomial,2,2>{},Polynomial{1});
  auto destination=std::move(source);
  require(!source.equals(identity),"moved-from zero must not become false 0/0 identity");
  require(!destination.equals(identity),"move value preserved");
  auto assigned=identity;assigned=std::move(destination);
  require(!assigned.equals(identity),"move assignment value preserved");
  require(!destination.denominator().is_zero(),"moved-from denominator invariant");
  auto u=SymbolicSu2::stereographic({Polynomial::var("t0"),Polynomial{},Polynomial{}},Chart::Plus);
  auto moved=std::move(u);
  require(!u.matrix().denominator().is_zero()&&u.norm_residual().is_zero(),"moved-from SU2 valid");
  u=std::move(moved);
  require(!moved.matrix().denominator().is_zero()&&moved.norm_residual().is_zero(),"move-assigned SU2 source valid");
}
void independent_taylor_exponential(){
  using M=Matrix<double,3,3>;
  const auto algebra=So3Algebra<double>({0.3,-0.7,0.2});
  auto term=M::identity(),sum=term;
  for(int k=1;k<=50;++k){term=(term*algebra.matrix()).scaled(1.0/k);sum=sum+term;}
  const auto closed=algebra.exp().matrix();
  for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j)close(sum(i,j),closed(i,j),2e-14);
}
int main(){
  const std::pair<const char*,std::function<void()>> tests[]={
    {"epsilon",epsilon_contract},{"matrices",matrices},{"exact_bases",exact_bases},
    {"numeric_float",numeric_groups<float>},{"numeric_double",numeric_groups<double>},
    {"numeric_long_double",numeric_groups<long double>},{"invalid_membership",invalid_membership},
    {"abstract_polynomials",abstract_polynomials},{"symbolic_exponentials",symbolic_exponentials},
    {"symbolic_charts",symbolic_charts},{"polynomial_homomorphism",polynomial_homomorphism},
    {"independent_taylor",independent_taylor_exponential},{"inverse_and_moves",checked_inverse_and_move_invariants}};
  int failed=0;for(const auto& [name,test]:tests){try{test();std::cout<<"PASS "<<name<<'\n';}catch(const std::exception& e){++failed;std::cerr<<"FAIL "<<name<<": "<<e.what()<<'\n';}}
  std::cout<<"test_groups="<<std::size(tests)<<" failed="<<failed<<'\n';return failed?1:0;
}
