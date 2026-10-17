#include "symbolic.hpp"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
using namespace mathphysics;
std::string quoted(const std::string& s){
  std::string out="\"";
  for(char c:s){if(c=='"'||c=='\\')out+='\\';out+=c;}
  return out+'"';
}
template<class T>std::string text(const T& v){std::ostringstream s;s<<v;return s.str();}
template<class T>std::string element(const T& v){return text(v);}
template<class T>std::string element(const Complex<T>& v){std::ostringstream os;os<<'('<<v.re<<")+I*("<<v.im<<')';return os.str();}
template<class T,std::size_t R,std::size_t C>
std::string matrix(const Matrix<T,R,C>& m){
  std::ostringstream out;out<<'[';
  for(std::size_t i=0;i<R;++i){if(i)out<<',';out<<'[';for(std::size_t j=0;j<C;++j){if(j)out<<',';out<<quoted(element(m(i,j)));}out<<']';}
  return out.str()+']';
}
template<class T>std::string numeric_matrix(const Matrix<T,3,3>& m){
  std::ostringstream out;out<<std::setprecision(17)<<'[';
  for(std::size_t i=0;i<3;++i){if(i)out<<',';out<<'[';for(std::size_t j=0;j<3;++j){if(j)out<<',';out<<m(i,j);}out<<']';}return out.str()+']';
}
void symbolic(){
  std::cout<<"{\"status\":\"COMPUTED\",\"epsilon3\":[";
  bool first=true;for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int k=0;k<3;++k){if(!first)std::cout<<',';first=false;std::cout<<epsilon3(i,j,k);}
  std::cout<<"],\"E\":[";
  for(std::size_t i=0;i<3;++i){if(i)std::cout<<',';std::cout<<matrix(So3Algebra<Exact>::basis(i).matrix());}
  std::cout<<"],\"sigma\":[";
  for(std::size_t i=0;i<3;++i){if(i)std::cout<<',';std::cout<<matrix(pauli<Exact>(i));}
  std::cout<<"],\"T\":[";
  for(std::size_t i=0;i<3;++i){if(i)std::cout<<',';std::cout<<matrix(generator<Exact>(i,Su2Basis::HalfPauliOverI));}
  std::cout<<"],\"F\":[";
  for(std::size_t i=0;i<3;++i){if(i)std::cout<<',';std::cout<<matrix(generator<Exact>(i,Su2Basis::PauliOverI));}
  const Vector3<Polynomial> v{Polynomial::var("v0"),Polynomial::var("v1"),Polynomial::var("v2")};
  const Vector3<Exact> axis{Exact{},Exact{},Exact{1}};
  const auto full=su2_exponential(axis,Su2Basis::PauliOverI),half=su2_exponential(axis,Su2Basis::HalfPauliOverI);
  const auto u=SymbolicSu2::stereographic({Polynomial::var("t0"),Polynomial::var("t1"),Polynomial::var("t2")},Chart::Plus);
  const auto m=u.matrix();const auto r=u.so3_matrix();
  std::cout<<"],\"hat\":"<<matrix(hat(v))<<",\"r_exp_z\":"<<matrix(so3_exponential(axis))
    <<",\"u_exp_full_z\":"<<matrix(full.matrix)<<",\"u_exp_half_z\":"<<matrix(half.matrix)
    <<",\"full_argument_scale\":"<<quoted(text(full.argument_scale))<<",\"half_argument_scale\":"<<quoted(text(half.argument_scale))
    <<",\"chart_u\":{\"numerator\":"<<matrix(m.numerator())<<",\"denominator\":"<<quoted(text(m.denominator()))<<'}'
    <<",\"chart_r\":{\"numerator\":"<<matrix(r.numerator())<<",\"denominator\":"<<quoted(text(r.denominator()))<<"}}\n";
}
void fixtures(){
  std::cout<<std::setprecision(17)<<"{\"status\":\"COMPUTED\",\"scalar\":\"f64\",\"cases\":[";
  for(int i=0;i<25;++i){
    const double x=i;
    const Vector3<double> v=i==0?Vector3<double>{}:Vector3<double>{0.13*(x-12),-0.07*(x+2),0.03*(i*i%17+1)};
    const auto r=So3Algebra<double>(v).exp();const auto u=Su2Algebra<double>(v).exp();
    if(i)std::cout<<',';
    std::cout<<"{\"id\":"<<i<<",\"v\":["<<v[0]<<','<<v[1]<<','<<v[2]<<"],\"rotation\":"<<numeric_matrix(r.matrix())<<",\"quaternion\":[";
    for(std::size_t j=0;j<4;++j){if(j)std::cout<<',';std::cout<<u.quaternion()[j];}
    std::cout<<"]}";
  }std::cout<<"]}\n";
}
int main(int argc,char** argv){
  try{if(argc==2&&std::string(argv[1])=="symbolic")symbolic();else if(argc==2&&std::string(argv[1])=="fixtures")fixtures();else throw std::invalid_argument("commands: symbolic | fixtures");return 0;}
  catch(const std::exception&){std::cout<<"{\"status\":\"INCONCLUSIVE\"}\n";return 2;}
}
