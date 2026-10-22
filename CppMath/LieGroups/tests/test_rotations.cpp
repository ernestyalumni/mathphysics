#include "rotations.hpp"
#include "symbolic.hpp"
#include <functional>
#include <iostream>
#include <string>
#include <vector>

using namespace mathphysics;
namespace {
void require(bool condition,const std::string& message){if(!condition)throw std::runtime_error(message);}
template<class F> void rejects(F&& f){bool caught=false;try{f();}catch(const std::exception&){caught=true;}require(caught,"expected rejection");}
template<std::floating_point T,std::size_t N,std::size_t M>
T error(const Matrix<T,N,M>& a,const Matrix<T,N,M>& b){T e{};for(std::size_t i=0;i<N;++i)for(std::size_t j=0;j<M;++j)e=std::max(e,std::abs(a(i,j)-b(i,j)));return e;}
template<std::floating_point T,std::size_t N,std::size_t M>
void close(const Matrix<T,N,M>& a,const Matrix<T,N,M>& b,T tol){require(error(a,b)<=tol,"matrix mismatch");}
template<std::floating_point T>
void close_vec(const Vector3<T>& a,const Vector3<T>& b,T tol){for(std::size_t i=0;i<3;++i)require(std::abs(a[i]-b[i])<=tol,"vector mismatch");}
template<std::floating_point T> Matrix<T,3,3> exp_reference(const Matrix<T,3,3>& a){
  auto sum=Matrix<T,3,3>::identity(),term=sum;
  for(int k=1;k<=100;++k){term=(term*a).scaled(T{1}/static_cast<T>(k));sum=sum+term;}
  return sum;
}
template<std::floating_point T> Vector3<T> geometric_rotation(Vector3<T> n,T phi,const Vector3<T>& v){
  const T length=norm(n);for(auto& x:n)x/=length;
  const T projection=dot(n,v),c=std::cos(phi),s=std::sin(phi);const auto tangent=cross(n,v);Vector3<T> out{};
  for(std::size_t j=0;j<3;++j)out[j]=n[j]*projection+c*(v[j]-n[j]*projection)+s*tangent[j];
  return out;
}
template<std::floating_point T> void scalar_contracts(){
  const T pi=std::numbers::pi_v<T>,eps=std::numeric_limits<T>::epsilon(),tol=T{1024}*eps;
  for(const auto& angles:std::vector<std::array<T,3>>{{T{0.4},T{0.7},T{-0.3}},{T{-2.1},T{-0.5},T{2.8}},{T{},T{},T{}},{T{0.4},pi,T{-0.9}}}){
    const auto [a,b,c]=angles;
    const auto zxz=from_euler_zxz(a,b,c),zyz=from_euler_zyz(a,b,c),zyx=from_yaw_pitch_roll(a,b,c);
    close(zxz.matrix(),quaternion_from_euler_zxz(a,b,c).to_so3().matrix(),tol);
    close(zyz.matrix(),quaternion_from_euler_zyz(a,b,c).to_so3().matrix(),tol);
    close(zyx.matrix(),quaternion_from_yaw_pitch_roll(a,b,c).to_so3().matrix(),tol);
    close(zxz.matrix(),exp_reference(hat(Vector3<T>{T{},T{},a}))*exp_reference(hat(Vector3<T>{b,T{},T{}}))*exp_reference(hat(Vector3<T>{T{},T{},c})),tol);
    close(zyz.matrix(),exp_reference(hat(Vector3<T>{T{},T{},a}))*exp_reference(hat(Vector3<T>{T{},b,T{}}))*exp_reference(hat(Vector3<T>{T{},T{},c})),tol);
    close(zyx.matrix(),exp_reference(hat(Vector3<T>{T{},T{},a}))*exp_reference(hat(Vector3<T>{T{},b,T{}}))*exp_reference(hat(Vector3<T>{c,T{},T{}})),tol);
    for(const auto& r:{zxz,zyz,zyx}){
      const auto aa=to_axis_angle(r);close(So3<T>::from_axis_angle(aa.axis,aa.angle).matrix(),r.matrix(),tol);
      close(quaternion_from_rotation(r).to_so3().matrix(),r.matrix(),tol);
      close(rotation_log(r).exp().matrix(),r.matrix(),tol);
    }
  }
  const Vector3<T> n{T{1},T{2},T{-3}},v{T{-2},T{4},T{5}};
  close_vec(So3<T>::from_axis_angle(n,T{0.8}).apply(v),geometric_rotation(n,T{0.8},v),tol*T{8});
}
void inverse_angles(){
  const double pi=std::numbers::pi;
  for(int k=0;k<97;++k){
    // All signs and noncanonical middle angles; compare represented matrices,
    // never angle triples that differ by valid convention branches.
    const double a=-3.7+0.071*k,b=-4.1+0.097*k,c=2.8-0.053*k;
    const auto x=from_euler_zxz(a,b,c),y=from_euler_zyz(a,b,c),z=from_yaw_pitch_roll(a,b,c);
    const auto u=to_euler_zxz(x),v=to_euler_zyz(y),w=to_yaw_pitch_roll(z);
    require(u.middle>=0&&u.middle<=pi&&v.middle>=0&&v.middle<=pi,"proper Euler branch");
    require(w.middle>=-pi/2&&w.middle<=pi/2,"Tait-Bryan branch");
    close(x.matrix(),from_euler_zxz(u.first,u.middle,u.last).matrix(),2e-14);
    close(y.matrix(),from_euler_zyz(v.first,v.middle,v.last).matrix(),2e-14);
    close(z.matrix(),from_yaw_pitch_roll(w.first,w.middle,w.last).matrix(),2e-14);
  }
}
void singular_angles(){
  const double pi=std::numbers::pi;
  for(double first:{-2.7,-0.4,0.8,2.9})for(double last:{-1.1,0.0,2.1}){
    for(double b:{0.0,pi}){
      const auto x=from_euler_zxz(first,b,last),y=from_euler_zyz(first,b,last);
      const auto u=to_euler_zxz(x),v=to_euler_zyz(y);
      require(u.singular&&v.singular&&u.last==0&&v.last==0,"proper Euler degeneracy");
      close(x.matrix(),from_euler_zxz(u.first,u.middle,u.last).matrix(),2e-14);
      close(y.matrix(),from_euler_zyz(v.first,v.middle,v.last).matrix(),2e-14);
    }
    for(double b:{-pi/2,pi/2}){
      const auto r=from_yaw_pitch_roll(first,b,last);const auto a=to_yaw_pitch_roll(r);
      require(a.singular&&a.last==0,"gimbal lock metadata");
      close(r.matrix(),from_yaw_pitch_roll(a.first,a.middle,a.last).matrix(),2e-14);
    }
    for(double b:{1e-10,pi-1e-10}){
      const auto r=from_euler_zxz(first,b,last);const auto a=to_euler_zxz(r);
      require(!a.singular,"resolvable near-singularity was collapsed");
      close(r.matrix(),from_euler_zxz(a.first,a.middle,a.last).matrix(),2e-14);
    }
  }
}
void axis_angle_endpoints(){
  const double pi=std::numbers::pi;
  const auto zero=to_axis_angle(So3<double>::identity());require(!zero.axis_line_unique&&zero.angle==0,"identity has no selected axis");
  for(Vector3<double> n:std::vector<Vector3<double>>{{1,0,0},{0,-1,0},{0,0,-1},{1,2,3},{-1,1,0}}){
    for(double p:{1e-15,1e-10,0.5,pi-1e-10,pi,pi+1e-10,5.8,-pi,-0.3}){
      const auto r=So3<double>::from_axis_angle(n,p);const auto a=to_axis_angle(r);
      require(a.angle>=0&&a.angle<=pi,"principal axis angle");
      require(a.axis_line_unique,"nonidentity has unique axis line");
      close(So3<double>::from_axis_angle(a.axis,a.angle).matrix(),r.matrix(),2e-14);
      close_vec(r.apply(a.axis),a.axis,2e-14);
      require(std::abs(r.matrix().trace()-(1+2*std::cos(a.angle)))<2e-14,"trace formula");
      require(std::abs(a.angle-std::abs(p))<1e-14||std::abs(p)>pi,"small-angle precision");
    }
  }
  const auto tiny=to_axis_angle(So3<double>::from_axis_angle({1,2,3},1e-15));
  require(std::abs(tiny.angle/1e-15-1)<1e-14,"trace-only acos would lose tiny angle");
  const auto half_turn=So3<double>::from_axis_angle({1,2,3},pi);
  const auto h=to_axis_angle(half_turn);auto negative_axis=h.axis;for(auto& x:negative_axis)x=-x;
  close(half_turn.matrix(),So3<double>::from_axis_angle(negative_axis,h.angle).matrix(),2e-14);
  require(h.axis_line_unique,"pi rotation has one axis line, two oriented axes");
  auto approximate=Matrix<double,3,3>::identity();approximate(0,0)+=1e-4;
  const auto loose=So3<double>::from_matrix(approximate,1e-3);
  require(quaternion_from_rotation(loose).quaternion()==std::array<double,4>{1,0,0,0},"extracted lift normalizes imported roundoff");
  const auto q=quaternion_from_euler_zxz(0.4,0.7,-0.3).quaternion();
  const auto antipode=quaternion_from_euler_zxz(0.4+2*pi,0.7,-0.3).quaternion();
  for(std::size_t j=0;j<4;++j)require(std::abs(q[j]+antipode[j])<1e-14,"2pi lift changes sign");
  for(double huge:{1e20,std::numeric_limits<double>::max()}){
    close(from_euler_zxz(huge,0.7,-0.3).matrix(),quaternion_from_euler_zxz(huge,0.7,-0.3).to_so3().matrix(),2e-14);
    close(from_euler_zyz(huge,-huge,0.8).matrix(),quaternion_from_euler_zyz(huge,-huge,0.8).to_so3().matrix(),2e-14);
  }
}
void ck_dictionary(){
  for(int k=0;k<25;++k){
    const auto q=quaternion_from_rotation(from_euler_zxz(0.11*k,0.07*k,-0.09*k));const auto ck=CayleyKlein<double>::from_quaternion(q);const auto m=ck.matrix();
    require(ck.quaternion().quaternion()==q.quaternion(),"CK preserves specified SU2 lift");
    require(CayleyKlein<double>::from_quaternion(q.antipode()).quaternion().quaternion()!=ck.quaternion().quaternion(),"CK does not identify opposite lifts");
    require(m(1,0)==-m(0,1).conjugate()&&m(1,1)==m(0,0).conjugate(),"CK conjugacy");
    const auto p=q.quaternion();require(m(0,0)==Complex<double>{p[0],p[3]}&&m(0,1)==Complex<double>{p[2],p[1]},"CK wxyz dictionary");
    const auto imported=CayleyKlein<double>::from_parameters(m(0,0),m(0,1),m(1,0),m(1,1));
    close(imported.goldstein_rotation().matrix(),q.to_so3().matrix().transpose(),2e-14);
    close(CayleyKlein<double>::from_quaternion(q.antipode()).goldstein_rotation().matrix(),ck.goldstein_rotation().matrix(),2e-14);
  }
}
void conjugation_and_ch(){
  const auto b=from_yaw_pitch_roll(0.4,-0.8,0.3);const auto a=Matrix<double,3,3>::from_rows({{{1,2,3},{-4,5,1},{-1,3,2}}});
  const Vector3<double> f{2,-1,3};const auto changed=coordinate_change(b,a);
  close_vec(changed.apply(b.apply(f)),b.apply(a.apply(f)),2e-13);
  require(std::abs(changed.trace()-a.trace())<2e-14,"similarity trace");
  const So3Algebra<double> x({0.7,-0.5,1.2});
  close(adjoint_action(b,x).matrix(),coordinate_change(b,x.matrix()),2e-14);
  close(adjoint_action(b,x).exp().matrix(),b.matrix()*x.exp().matrix()*b.matrix().transpose(),2e-14);
  close(cayley_hamilton_residual(a),Matrix<double,3,3>{},1e-12);
  close(cayley_hamilton_residual(b.matrix()),Matrix<double,3,3>{},2e-14);
  const auto e=hat(Vector3<double>{1,2,-3});close(e*e*e,e.scaled(-14.0),1e-14);
  require(error(a*a*a-a*a+ a,Matrix<double,3,3>{})>1,"wrong characteristic polynomial witness");
  // Same displayed matrix, two meanings: passive same-vector conversion and
  // active vector rotation. They agree numerically, their domains differ.
  const auto body=Vector3<double>{0.3,0.7,-0.2};const auto inertial=b.apply(body);
  close_vec(b.inverse().apply(inertial),body,2e-14);
  require(norm(Vector3<double>{inertial[0]-body[0],inertial[1]-body[1],inertial[2]-body[2]})>0.01,"nontrivial active witness");
}
void exact_polynomial_and_trig(){
  Matrix<Polynomial,3,3> m;
  for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j)m(i,j)=Polynomial::var("m"+std::to_string(i)+std::to_string(j));
  require(cayley_hamilton_residual(m)==Matrix<Polynomial,3,3>{},"exact general characteristic identity");
  const Vector3<Polynomial> x{Polynomial::var("x"),Polynomial::var("y"),Polynomial::var("z")};const auto k=hat(x);
  require(k*k*k==k.scaled(-dot(x,x)),"exact hat cubic");
  const auto r=so3_exponential({Exact{1}/3,Exact{2}/3,Exact{2}/3});
  require(cayley_hamilton_residual(r)==Matrix<Trig,3,3>{},"exact all-angle rotation characteristic identity");
  const std::array<Polynomial,4> q{Polynomial::var("w"),x[0],x[1],x[2]};
  Polynomial norm2;for(const auto& v:q)norm2=norm2+v*v;const auto ck=cayley_klein_matrix(q);
  require(adjoint(ck)*ck==Matrix<Complex<Polynomial>,2,2>::identity().scaled(Complex<Polynomial>{norm2}),"exact CK norm");
  require(quaternion_rotation(q).determinant()==norm2*norm2*norm2,"exact quaternion determinant");
}
void invalid_inputs(){
  const double inf=std::numeric_limits<double>::infinity(),nan=std::numeric_limits<double>::quiet_NaN();
  rejects([&]{(void)from_euler_zxz(nan,0.0,0.0);});rejects([&]{(void)from_euler_zyz(0.0,inf,0.0);});rejects([&]{(void)from_yaw_pitch_roll(0.0,0.0,nan);});
  rejects([&]{(void)quaternion_from_euler_zxz(nan,0.0,0.0);});rejects([&]{(void)quaternion_from_euler_zyz(0.0,inf,0.0);});rejects([&]{(void)quaternion_from_yaw_pitch_roll(0.0,0.0,nan);});
  rejects([]{(void)axis_rotation<double>(3,0.3);});
  const auto r=So3<double>::identity();
  for(double tol:{-1.0,0.125,inf,nan}){rejects([&]{(void)to_euler_zxz(r,tol);});rejects([&]{(void)to_euler_zyz(r,tol);});rejects([&]{(void)to_yaw_pitch_roll(r,tol);});}
  rejects([]{(void)CayleyKlein<double>::from_parameters({1,0},{0,0},{0,1},{1,0});});
  rejects([]{(void)CayleyKlein<double>::from_parameters({2,0},{0,0},{0,0},{2,0});});
  rejects([&]{(void)CayleyKlein<double>::from_parameters({nan,0},{0,0},{0,0},{1,0});});
  auto m=Matrix<double,3,3>::identity();m(0,1)=nan;rejects([&]{(void)coordinate_change(r,m);});
  auto reflection=Matrix<double,3,3>::identity();reflection(0,0)=-1;rejects([&]{(void)So3<double>::from_matrix(reflection);});
  // Negative controls distinguish Goldstein A from R and noncommuting orders.
  const auto z=axis_rotation<double>(2,0.4),x=axis_rotation<double>(0,0.7);
  require(error(z.compose(x).matrix(),x.compose(z).matrix())>0.1,"order negative control");
  require(error(z.matrix(),z.matrix().transpose())>0.5,"passive sign negative control");
}
}
int main(){
  const std::vector<std::pair<std::string,std::function<void()>>> groups{
    {"float forward/axis-angle",scalar_contracts<float>},{"double forward/axis-angle",scalar_contracts<double>},{"long double forward/axis-angle",scalar_contracts<long double>},
    {"inverse Euler branches",inverse_angles},{"Euler singular endpoints",singular_angles},{"axis-angle zero/pi robustness",axis_angle_endpoints},
    {"Goldstein CK quaternion dictionary",ck_dictionary},{"operator conjugation, Ad and Cayley-Hamilton",conjugation_and_ch},{"exact polynomial/trig certificates",exact_polynomial_and_trig},{"invalid inputs and sign/order controls",invalid_inputs}};
  for(const auto& [name,f]:groups){try{f();std::cout<<"PASS "<<name<<'\n';}catch(const std::exception& e){std::cerr<<"FAIL "<<name<<": "<<e.what()<<'\n';return 1;}}
  std::cout<<groups.size()<<" Goldstein test groups passed\n";
}
