#include "rotations.hpp"
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
using namespace mathphysics;
int main(){
  std::cout<<std::setprecision(17);const double pi=std::numbers::pi;
  for(const std::string convention:{"zxz","zyz","zyx"}){
    const bool tb=convention=="zyx";
    const std::vector<std::array<double,3>> fixtures{{0.4,0.7,-0.3},{-2.1,-0.5,2.8},{0.4,0.0,-0.3},{0.4,tb?pi/2:pi,-0.3},{0.4,tb?-pi/2:1e-10,-0.3},{0.4,tb?pi/2-1e-10:pi-1e-10,-0.3}};
    for(std::size_t i=0;i<fixtures.size();++i){
      const auto [a,b,c]=fixtures[i];
      const auto r=convention=="zxz"?from_euler_zxz(a,b,c):convention=="zyz"?from_euler_zyz(a,b,c):from_yaw_pitch_roll(a,b,c);
      const auto q=quaternion_from_rotation(r).quaternion();const auto aa=to_axis_angle(r);
      const auto e=convention=="zxz"?to_euler_zxz(r):convention=="zyz"?to_euler_zyz(r):to_yaw_pitch_roll(r);
      const auto ck=cayley_klein_matrix(q);
      std::cout<<"{\"name\":\""<<convention<<'_'<<i<<"\",\"convention\":\""<<convention<<"\",\"angles\":["<<a<<','<<b<<','<<c<<"],\"matrix\":[";
      for(std::size_t j=0;j<3;++j){if(j)std::cout<<',';std::cout<<'[';for(std::size_t k=0;k<3;++k){if(k)std::cout<<',';std::cout<<r.matrix()(j,k);}std::cout<<']';}
      std::cout<<"],\"quaternion\":[";for(std::size_t j=0;j<4;++j){if(j)std::cout<<',';std::cout<<q[j];}
      std::cout<<"],\"axis\":["<<aa.axis[0]<<','<<aa.axis[1]<<','<<aa.axis[2]<<"],\"angle\":"<<aa.angle<<",\"extracted\":["<<e.first<<','<<e.middle<<','<<e.last<<"],\"singular\":"<<(e.singular?"true":"false")<<",\"ck\":[";
      for(std::size_t j=0;j<2;++j)for(std::size_t k=0;k<2;++k){if(j||k)std::cout<<',';std::cout<<'['<<ck(j,k).re<<','<<ck(j,k).im<<']';}
      std::cout<<"]}\n";
    }
  }
}
