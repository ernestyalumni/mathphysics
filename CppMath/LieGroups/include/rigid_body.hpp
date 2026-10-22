#pragma once
// Numerical COM-centered 6DOF dynamics; inertial r,v and body omega/inertia.
// Quaternion: scalar-first Hamilton, active body -> inertial. SI-compatible;
// no gravity or external wrench is implicit. A step holds supplied wrench fixed.
#include "lie_groups.hpp"
#include <algorithm>
#include <numbers>

namespace mathphysics::rigid_body {
using Vec=Vector3<double>;
using Mat=Matrix<double,3,3>;
struct State {Vec position{},velocity{};std::array<double,4> quaternion{1,0,0,0};Vec omega_body{};};
struct Wrench {Vec force_body{},torque_body{},force_inertial{};};
struct Derivative {Vec position,velocity;std::array<double,4> quaternion;Vec omega_body;};
struct StepOptions {double tolerance=1e-12;std::size_t max_iterations=20;};
struct StepResult {State state;std::size_t iterations;double residual;};
namespace detail {
inline void finite(const Vec& v){for(auto x:v)if(!std::isfinite(x))throw std::invalid_argument("nonfinite dynamics vector");}
inline double inf_norm(const Vec& v){return std::max({std::abs(v[0]),std::abs(v[1]),std::abs(v[2])});}
inline Vec add(const Vec& a,const Vec& b){return {a[0]+b[0],a[1]+b[1],a[2]+b[2]};}
inline Vec sub(const Vec& a,const Vec& b){return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};}
inline Vec scale(const Vec& a,double x){return {x*a[0],x*a[1],x*a[2]};}
inline std::array<double,4> normalized(std::array<double,4> q){
  const double n=std::hypot(std::hypot(q[0],q[1]),std::hypot(q[2],q[3]));
  if(!std::isfinite(n)||n==0)throw std::invalid_argument("invalid internal quaternion");
  for(auto& x:q)x/=n;
  return q;
}
inline State validated(State s){
  finite(s.position);finite(s.velocity);finite(s.omega_body);
  (void)Su2<double>::from_unit_quaternion(s.quaternion,1e-10);
  s.quaternion=normalized(s.quaternion);return s;
}
inline void validate_wrench(const Wrench& w){finite(w.force_body);finite(w.torque_body);finite(w.force_inertial);}
// Pivoted solve for a general (Newton Jacobian) matrix, not assumed SPD.
inline Vec solve3(Mat a,Vec b){
  for(std::size_t k=0;k<3;++k){
    std::size_t p=k;for(std::size_t i=k+1;i<3;++i)if(std::abs(a(i,k))>std::abs(a(p,k)))p=i;
    if(!std::isfinite(a(p,k))||a(p,k)==0)throw std::runtime_error("singular Newton Jacobian");
    if(p!=k){for(std::size_t j=k;j<3;++j)std::swap(a(p,j),a(k,j));std::swap(b[p],b[k]);}
    for(std::size_t i=k+1;i<3;++i){const double f=a(i,k)/a(k,k);for(std::size_t j=k+1;j<3;++j)a(i,j)-=f*a(k,j);b[i]-=f*b[k];}
  }
  Vec x{};for(int i=2;i>=0;--i){double v=b[std::size_t(i)];for(std::size_t j=std::size_t(i)+1;j<3;++j)v-=a(std::size_t(i),j)*x[j];x[std::size_t(i)]=v/a(std::size_t(i),std::size_t(i));}
  finite(x);return x;
}
inline double largest_symmetric_eigenvalue(const Mat& a){
  const double p1=a(0,1)*a(0,1)+a(0,2)*a(0,2)+a(1,2)*a(1,2);
  if(p1==0)return std::max({a(0,0),a(1,1),a(2,2)});
  const double q=a.trace()/3;auto b=a-Mat::identity().scaled(q);
  const double p=std::sqrt((b(0,0)*b(0,0)+b(1,1)*b(1,1)+b(2,2)*b(2,2)+2*p1)/6);
  if(p==0)return q;
  b=b.scaled(1/p);const double angle=std::acos(std::clamp(b.determinant()/2,-1.0,1.0))/3;
  return q+2*p*std::cos(angle);
}
} // namespace detail

class Model {
public:
  // J is the full inertia about the COM, in body axes, not necessarily principal.
  // Tiny symmetry roundoff (1e-12 relative) is averaged; nonphysical/SPD failures reject.
  Model(double mass,Mat inertia):mass_(mass),inertia_(std::move(inertia)){
    if(!std::isfinite(mass_)||mass_<=0)throw std::invalid_argument("mass must be positive finite");
    for(const auto& row:inertia_.rows())for(auto x:row){if(!std::isfinite(x))throw std::invalid_argument("nonfinite inertia");scale_=std::max(scale_,std::abs(x));}
    if(scale_==0)throw std::invalid_argument("inertia must be SPD");
    for(std::size_t i=0;i<3;++i)for(std::size_t j=i+1;j<3;++j){
      if(std::abs(inertia_(i,j)/scale_-inertia_(j,i)/scale_)>1e-12)throw std::invalid_argument("inertia must be symmetric");
      const double mean=inertia_(i,j)/2+inertia_(j,i)/2;inertia_(i,j)=mean;inertia_(j,i)=mean;
    }
    Mat a;for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j)a(i,j)=inertia_(i,j)/scale_;
    for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<=i;++j){
      double x=a(i,j);for(std::size_t k=0;k<j;++k)x-=chol_(i,k)*chol_(j,k);
      if(i==j){if(!std::isfinite(x)||x<=0)throw std::invalid_argument("inertia must be SPD");chol_(i,j)=std::sqrt(x);}else chol_(i,j)=x/chol_(j,j);
    }
    if(2*detail::largest_symmetric_eigenvalue(a)>a.trace()+1e-12)throw std::invalid_argument("inertia violates physical triangle condition");
    for(std::size_t j=0;j<3;++j){Vec e{};e[j]=1;const auto col=solve_inertia(e);for(std::size_t i=0;i<3;++i)inverse_(i,j)=col[i];}
  }
  double mass()const{return mass_;}
  const Mat& inertia()const{return inertia_;}
  Vec solve_inertia(const Vec& b)const{
    detail::finite(b);Vec y{},x{};
    for(std::size_t i=0;i<3;++i){double v=b[i]/scale_;for(std::size_t j=0;j<i;++j)v-=chol_(i,j)*y[j];y[i]=v/chol_(i,i);}
    for(int i=2;i>=0;--i){double v=y[std::size_t(i)];for(std::size_t j=std::size_t(i)+1;j<3;++j)v-=chol_(j,std::size_t(i))*x[j];x[std::size_t(i)]=v/chol_(std::size_t(i),std::size_t(i));}
    detail::finite(x);return x;
  }
  Vec momentum_body(const State& state)const{const auto s=detail::validated(state);const auto h=inertia_.apply(s.omega_body);detail::finite(h);return h;}
  Vec momentum_inertial(const State& state)const{
    const auto s=detail::validated(state);return Su2<double>::from_unit_quaternion(s.quaternion).to_so3().apply(momentum_body(s));
  }
  double kinetic_energy(const State& state)const{
    const auto s=detail::validated(state);const double e=0.5*(mass_*dot(s.velocity,s.velocity)+dot(s.omega_body,momentum_body(s)));
    if(!std::isfinite(e))throw std::overflow_error("nonfinite energy");
    return e;
  }
  Derivative rhs(const State& state,const Wrench& wrench={})const{
    const auto s=detail::validated(state);detail::validate_wrench(wrench);
    const auto r=Su2<double>::from_unit_quaternion(s.quaternion).to_so3();
    auto acceleration=detail::add(r.apply(wrench.force_body),wrench.force_inertial);for(auto& x:acceleration)x/=mass_;
    auto dq=hamilton(s.quaternion,std::array<double,4>{0,s.omega_body[0],s.omega_body[1],s.omega_body[2]});for(auto& x:dq)x*=0.5;
    const auto angular=solve_inertia(detail::sub(wrench.torque_body,cross(s.omega_body,momentum_body(s))));
    detail::finite(acceleration);for(auto x:dq)if(!std::isfinite(x))throw std::overflow_error("quaternion RHS overflow");
    return {s.velocity,acceleration,dq,angular};
  }
  // Implicit midpoint in h=J omega, Cayley attitude increment, midpoint COM step.
  // No hidden dt cap; a nonconvergent/singular/overflowing step fails explicitly.
  StepResult step(const State& state,double dt,const Wrench& wrench={},StepOptions options={})const{
    const auto s=detail::validated(state);detail::validate_wrench(wrench);
    if(!std::isfinite(dt)||dt<=0)throw std::invalid_argument("dt must be positive finite");
    if(!std::isfinite(options.tolerance)||options.tolerance<=0||options.tolerance>1e-4||options.max_iterations==0||options.max_iterations>100)throw std::invalid_argument("invalid Newton options");
    const auto h0=momentum_body(s);
    auto h1=detail::add(h0,detail::scale(detail::sub(wrench.torque_body,cross(s.omega_body,h0)),dt));
    const double threshold=options.tolerance*(1+detail::inf_norm(h0)+dt*detail::inf_norm(wrench.torque_body));
    if(!std::isfinite(threshold))throw std::overflow_error("Newton threshold overflow");
    Vec wm{};double residual=0;std::size_t iteration=0;
    for(;;){
      detail::finite(h1);const auto hm=detail::scale(detail::add(h0,h1),0.5);wm=solve_inertia(hm);
      const auto f=detail::sub(detail::sub(h1,h0),detail::scale(detail::sub(wrench.torque_body,cross(wm,hm)),dt));
      detail::finite(f);residual=detail::inf_norm(f);if(residual<=threshold)break;
      if(iteration>=options.max_iterations)throw std::runtime_error("implicit midpoint Newton did not converge");
      const auto jac=Mat::identity()+(hat(wm)-hat(hm)*inverse_).scaled(dt/2);
      h1=detail::sub(h1,detail::solve3(jac,f));++iteration;
    }
    auto increment=detail::normalized({1,dt*wm[0]/2,dt*wm[1]/2,dt*wm[2]/2});
    const double half_scalar=std::sqrt((1+increment[0])/2);
    const std::array<double,4> half{half_scalar,increment[1]/(2*half_scalar),increment[2]/(2*half_scalar),increment[3]/(2*half_scalar)};
    const auto qmid=detail::normalized(hamilton(s.quaternion,half));
    State out=s;out.quaternion=detail::normalized(hamilton(s.quaternion,increment));out.omega_body=solve_inertia(h1);
    const auto force=detail::add(Su2<double>::from_unit_quaternion(qmid).to_so3().apply(wrench.force_body),wrench.force_inertial);
    auto dv=force;for(auto& x:dv)x=(x/mass_)*dt;out.velocity=detail::add(s.velocity,dv);
    out.position=detail::add(s.position,detail::scale(detail::add(s.velocity,out.velocity),dt/2));
    out=detail::validated(out);return {out,iteration,residual};
  }
private:
  double mass_,scale_=0;Mat inertia_,chol_,inverse_;
};
} // namespace mathphysics::rigid_body
