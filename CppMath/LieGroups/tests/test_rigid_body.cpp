#include "rigid_body.hpp"
#include "rotations.hpp"
#include <functional>
#include <iostream>
#include <vector>
using namespace mathphysics;
using namespace mathphysics::rigid_body;
namespace {
void check(bool x,const char* msg){if(!x)throw std::runtime_error(msg);}
void close(const Vec& a,const Vec& b,double tol){for(std::size_t i=0;i<3;++i)check(std::abs(a[i]-b[i])<=tol,"vector mismatch");}
double error(const Mat& a,const Mat& b){double x=0;for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j)x=std::max(x,std::abs(a(i,j)-b(i,j)));return x;}
template<class F> void rejects(F f){bool yes=false;try{f();}catch(const std::exception&){yes=true;}check(yes,"expected explicit failure");}
Mat diagonal(double a,double b,double c){return Mat::from_rows({{{a,0,0},{0,b,0},{0,0,c}}});}
Mat offdiag(){return Mat::from_rows({{{2,.2,.1},{.2,3,.3},{.1,.3,4}}});}
State evolve(const Model& m,State s,double dt,std::size_t n,const Wrench& w={}){for(std::size_t i=0;i<n;++i)s=m.step(s,dt,w).state;return s;}
Mat attitude(const State& s){return Su2<double>::from_unit_quaternion(s.quaternion).to_so3().matrix();}
void translation(){
  const Model m(2,diagonal(1,2,2.5));State s;s.position={.1,-.2,.3};s.velocity={.3,.1,-.2};const auto initial=s;
  s=evolve(m,s,.01,100);for(std::size_t i=0;i<3;++i)check(std::abs(s.position[i]-initial.position[i]-initial.velocity[i])<1e-14,"exact free translation");
  close(s.velocity,initial.velocity,0);check(s.quaternion==initial.quaternion,"unforced rest attitude");
}
void sphere_and_principal(){
  for(bool sphere:{false,true}){
    const Model m(2,sphere?diagonal(2,2,2):diagonal(1,2,2.5));State s;s.omega_body=sphere?Vec{.2,-.4,.7}:Vec{0,0,.7};
    const auto w=s.omega_body;const double n=norm(w),dt=.01;const std::size_t count=100;
    s=evolve(m,s,dt,count);close(s.omega_body,w,2e-14);
    const auto expected=So3<double>::from_axis_angle(w,static_cast<double>(count)*2*std::atan(dt*n/2));
    check(error(attitude(s),expected.matrix())<5e-14,"constant-omega Cayley reference");
  }
}
void full_inertia_conservation(){
  const Model m(2,offdiag());State s;s.omega_body={.3,-.4,.8};s.velocity={.2,-.1,.3};s.quaternion=quaternion_from_euler_zxz(.4,.7,-.3).quaternion();
  const auto h=m.momentum_inertial(s);const double energy=m.kinetic_energy(s),hnorm=norm(m.momentum_body(s));
  for(int k=0;k<2000;++k){s=m.step(s,.005,{},StepOptions{1e-14,20}).state;close(m.momentum_inertial(s),h,2e-10);check(std::abs(m.kinetic_energy(s)-energy)<2e-10,"free energy invariant");check(std::abs(norm(m.momentum_body(s))-hnorm)<2e-10,"body momentum norm invariant");}
  close(m.solve_inertia(m.inertia().apply({.3,-.7,1.1})),Vec{.3,-.7,1.1},2e-15);
}
void intermediate_axis_flip(){
  const Model m(1,diagonal(1,2,2.5));State s;s.omega_body={.01,1,.01};const auto h=m.momentum_inertial(s);const double energy=m.kinetic_energy(s);
  s=evolve(m,s,.005,4800);check(s.omega_body[1]<-.99,"unforced intermediate-axis flip");
  close(m.momentum_inertial(s),h,2e-9);check(std::abs(m.kinetic_energy(s)-energy)<2e-9,"flip conserved energy");
}
Mat symmetric_top_reference(const State& s,double time,double inertia,double axial){
  const Model m(1,diagonal(inertia,inertia,axial));const auto h=m.momentum_inertial(s);const double a=(axial-inertia)*s.omega_body[2]/inertia;
  return So3<double>::from_axis_angle(h,norm(h)*time/inertia).matrix()*attitude(s)*axis_rotation<double>(2,-a*time).matrix();
}
void symmetric_top_convergence(){
  const Model m(1,diagonal(2,2,3));State s;s.omega_body={.3,-.5,.8};s.quaternion=quaternion_from_yaw_pitch_roll(.4,-.2,.3).quaternion();
  const auto exact=symmetric_top_reference(s,1,2,3);double previous=0;
  for(std::size_t n:{25u,50u,100u}){
    const auto out=evolve(m,s,1.0/static_cast<double>(n),n);const double e=error(attitude(out),exact);
    const double a=.4;const Vec omega{.3*std::cos(a)+.5*std::sin(a),.3*std::sin(a)-.5*std::cos(a),.8};
    close(out.omega_body,omega,2e-5);
    if(previous>0)check(previous/e>3.8&&previous/e<4.2,"second-order attitude convergence");
    previous=e;
  }
}
void forcing_and_rhs(){
  const Model m(2,diagonal(2,2,2));State s;Wrench w;w.force_body={2,-1,.5};s=evolve(m,s,.01,100,w);
  close(s.velocity,{1,-.5,.25},2e-14);close(s.position,{.5,-.25,.125},2e-14);
  s={};w={};w.torque_body={.1,-.2,.3};s=evolve(m,s,.01,100,w);close(s.omega_body,{.05,-.1,.15},2e-14);
  State tilted;tilted.quaternion=quaternion_from_yaw_pitch_roll(std::numbers::pi/2,0.0,0.0).quaternion();w={};w.force_body={2,0,0};
  close(m.rhs(tilted,w).velocity,{0,1,0},1e-14);
  tilted.omega_body={.2,-.4,.7};w={};w.force_inertial={2,-1,.5};s=evolve(m,tilted,.01,100,w);close(s.velocity,{1,-.5,.25},2e-14);close(s.position,{.5,-.25,.125},2e-14);
  const Model asymmetric(1,diagonal(1,2,2.5));s={};s.omega_body={1,2,3};close(asymmetric.rhs(s).omega_body,{-3,2.25,-.8},2e-14);
  const auto dq=asymmetric.rhs(s).quaternion;check(dq==std::array<double,4>{0,.5,1,1.5},"Hamilton quaternion derivative");
}
void lift_and_discrete_work(){
  const Model m(2,offdiag());State a;a.quaternion=quaternion_from_euler_zxz(.2,.5,-.4).quaternion();a.omega_body={.3,-.4,.8};a.velocity={.2,.1,-.3};State b=a;for(auto& x:b.quaternion)x=-x;
  Wrench w;w.force_body={.2,-.1,.3};w.torque_body={.03,-.02,.01};
  const double initial=m.kinetic_energy(a),dt=.02;const auto out=m.step(a,dt,w).state;
  Vec vm{},wm{};for(std::size_t i=0;i<3;++i){vm[i]=(a.velocity[i]+out.velocity[i])/2;wm[i]=(a.omega_body[i]+out.omega_body[i])/2;}
  auto qm=a.quaternion;for(std::size_t i=0;i<4;++i)qm[i]+=out.quaternion[i];qm=mathphysics::rigid_body::detail::normalized(qm);
  const auto force=Su2<double>::from_unit_quaternion(qm).to_so3().apply(w.force_body);
  const double work=dt*(dot(vm,force)+dot(wm,w.torque_body));check(std::abs(m.kinetic_energy(out)-initial-work)<2e-12,"midpoint wrench work identity");
  a=evolve(m,a,.002,100,w);b=evolve(m,b,.002,100,w);close(a.position,b.position,0);close(a.velocity,b.velocity,0);close(a.omega_body,b.omega_body,0);
  for(std::size_t i=0;i<4;++i)check(std::abs(a.quaternion[i]+b.quaternion[i])<1e-15,"q/-q equivariance");
}
void invalid(){
  const auto j=diagonal(1,2,2.5);const double nan=std::numeric_limits<double>::quiet_NaN(),inf=std::numeric_limits<double>::infinity();
  for(double mass:{0.,-1.,nan,inf})rejects([&]{(void)Model(mass,j);});
  rejects([]{(void)Model(1,diagonal(1,1,3));});rejects([]{(void)Model(1,diagonal(0,1,1));});rejects([]{(void)Model(1,diagonal(-1,2,2));});
  rejects([]{(void)Model(1,Mat::from_rows({{{1,2,0},{2,1,0},{0,0,1}}}));});
  const auto rotation=from_yaw_pitch_roll(.4,.7,-.3).matrix();
  rejects([&]{(void)Model(1,rotation*diagonal(1,1,3)*rotation.transpose());});
  (void)Model(1,rotation*diagonal(1,1,2)*rotation.transpose());
  (void)Model(1,diagonal(1,1,2));auto nonsymmetric=j;nonsymmetric(0,1)=.1;rejects([&]{(void)Model(1,nonsymmetric);});
  auto nonfinite=j;nonfinite(0,0)=nan;rejects([&]{(void)Model(1,nonfinite);});const Model m(1,j);State s;
  for(double dt:{0.,-1.,nan,inf})rejects([&]{(void)m.step(s,dt);});
  s.quaternion={2,0,0,0};rejects([&]{(void)m.rhs(s);});rejects([&]{(void)m.step(s,.01);});
  s={};s.velocity[0]=nan;rejects([&]{(void)m.rhs(s);});s={};Wrench w;w.torque_body[0]=inf;rejects([&]{(void)m.step(s,.01,w);});
  rejects([&]{(void)m.step(s,.01,{},StepOptions{0,20});});rejects([&]{(void)m.step(s,.01,{},StepOptions{1e-12,0});});
  rejects([&]{(void)m.step(s,.01,{},StepOptions{1e-3,20});});rejects([&]{(void)m.step(s,.01,{},StepOptions{1e-12,101});});
  s.omega_body={3,2,4};rejects([&]{(void)m.step(s,1.2,{},StepOptions{1e-14,1});});
  s={};s.quaternion[0]+=1e-12;check(m.step(s,.01).state.quaternion==std::array<double,4>{1,0,0,0},"accepted quaternion roundoff is normalized");
}
}
int main(){const std::vector<std::pair<const char*,std::function<void()>>> tests{{"free translation",translation},{"spherical and principal rotations",sphere_and_principal},{"full SPD inertia conservation",full_inertia_conservation},{"intermediate-axis flip",intermediate_axis_flip},{"symmetric-top reference and second-order convergence",symmetric_top_convergence},{"body/inertial forcing and RHS signs",forcing_and_rhs},{"lift equivariance and discrete work",lift_and_discrete_work},{"invalid inputs and solver failure",invalid}};
for(const auto& [name,f]:tests){try{f();std::cout<<"PASS "<<name<<'\n';}catch(const std::exception& e){std::cerr<<"FAIL "<<name<<": "<<e.what()<<'\n';return 1;}}std::cout<<tests.size()<<" rigid-body groups passed\n";}
