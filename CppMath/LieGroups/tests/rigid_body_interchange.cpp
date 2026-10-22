#include "rigid_body.hpp"
#include <iomanip>
#include <iostream>
#include <string>
using namespace mathphysics;
using namespace mathphysics::rigid_body;
namespace {
template<std::size_t N> void array(const std::array<double,N>& v){std::cout<<'[';for(std::size_t i=0;i<N;++i){if(i)std::cout<<',';std::cout<<v[i];}std::cout<<']';}
void state(const State& s){std::cout<<"{\"r\":";array(s.position);std::cout<<",\"v\":";array(s.velocity);std::cout<<",\"q\":";array(s.quaternion);std::cout<<",\"omega\":";array(s.omega_body);std::cout<<'}';}
}
int main(){std::cout<<std::setprecision(17);
for(int k=0;k<8;++k){
  const std::array<std::string,8> names{"free_sphere","free_asymmetric","free_offdiag","principal","forced_translation","forced_torque","simultaneous","nonidentity"};
  Mat j=Mat::from_rows({{{2,0,0},{0,2,0},{0,0,2}}});
  if(k==1||k==3)j=Mat::from_rows({{{1,0,0},{0,2,0},{0,0,2.5}}});
  if(k==2||k==6||k==7)j=Mat::from_rows({{{2,.2,.1},{.2,3,.3},{.1,.3,4}}});
  State s;s.position={.1,-.2,.3};s.velocity={.3,.1,-.2};s.omega_body={.2,-.4,.7};
  if(k==1)s.omega_body={.2,1,.3};
  if(k==2||k==6||k==7)s.omega_body={.3,-.4,.8};
  if(k==3)s.omega_body={0,0,.7};
  if(k==4||k==5)s.omega_body={0,0,0};
  if(k==7)s.quaternion={std::cos(.3),std::sin(.3),0,0};
  Wrench w;
  if(k==4)w.force_body={2,-1,.5};
  if(k==5)w.torque_body={.1,-.2,.3};
  if(k==6){w.force_body={.2,-.1,.3};w.torque_body={.03,-.02,.01};}
  const double dt=k==1?.005:k==6?.002:.01;const int count=k==1?400:100;const auto initial=s;const Model m(2,j);
  for(int i=0;i<count;++i)s=m.step(s,dt,w).state;
  std::cout<<"{\"name\":\""<<names[std::size_t(k)]<<"\",\"mass\":2,\"inertia\":[";
  for(std::size_t i=0;i<3;++i){if(i)std::cout<<',';array(j.rows()[i]);}
  std::cout<<"],\"dt\":"<<dt<<",\"steps\":"<<count<<",\"initial\":";state(initial);
  std::cout<<",\"wrench\":{\"force\":";array(w.force_body);std::cout<<",\"torque\":";array(w.torque_body);std::cout<<",\"force_inertial\":";array(w.force_inertial);std::cout<<"},\"final\":";state(s);
  std::cout<<",\"energy\":"<<m.kinetic_energy(s)<<",\"H_I\":";array(m.momentum_inertial(s));std::cout<<"}\n";
}}
