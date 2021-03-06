#pragma once

#include <initializer_list>
#include <boost/fusion/container/generation/make_vector.hpp>
#include <boost/fusion/include/make_vector.hpp>

#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Utils.h" 
#include "UHH2/core/include/AnalysisModule.h"


template<typename P>
inline void sort_by_eta(std::vector<P> & particles){
  std::sort(particles.begin(), particles.end(), [](const P & p1, const P & p2){return fabs(p1.eta()) > fabs(p2.eta());});
}

template<typename P>
inline void sort_by_eta(std::vector<P*> & particles){
  std::sort(particles.begin(), particles.end(), [](const P* p1, const P* p2){return fabs(p1->eta()) > fabs(p2->eta());});
}

inline std::vector<std::unique_ptr<Selection>> make_uvec(std::unique_ptr<Selection> a, std::unique_ptr<Selection> b){
  std::vector<std::unique_ptr<Selection>> my_vec;
  my_vec.push_back(move(a)); my_vec.push_back(move(b));
  return my_vec;
}

inline std::vector<std::unique_ptr<Selection>> make_uvec(std::unique_ptr<Selection> a, std::unique_ptr<Selection> b, std::unique_ptr<Selection> c){
  std::vector<std::unique_ptr<Selection>> my_vec;
  my_vec.push_back(move(a)); my_vec.push_back(move(b)); my_vec.push_back(move(c));
  return my_vec;
}

inline std::vector<std::unique_ptr<Selection>> make_uvec(std::unique_ptr<Selection> a, std::unique_ptr<Selection> b, std::unique_ptr<Selection> c, std::unique_ptr<Selection> d){
  std::vector<std::unique_ptr<Selection>> my_vec;
  my_vec.push_back(move(a)); my_vec.push_back(move(b)); my_vec.push_back(move(c)); my_vec.push_back(move(d));
  return my_vec;
}
inline std::vector<std::unique_ptr<Selection>> make_uvec(std::unique_ptr<Selection> a, std::unique_ptr<Selection> b, std::unique_ptr<Selection> c, std::unique_ptr<Selection> d, std::unique_ptr<Selection> e){
  std::vector<std::unique_ptr<Selection>> my_vec;
  my_vec.push_back(move(a)); my_vec.push_back(move(b)); my_vec.push_back(move(c)); my_vec.push_back(move(d)); my_vec.push_back(move(e));
  return my_vec;
}
inline std::vector<std::unique_ptr<Selection>> make_uvec(std::unique_ptr<Selection> a, std::unique_ptr<Selection> b, std::unique_ptr<Selection> c, std::unique_ptr<Selection> d, std::unique_ptr<Selection> e, std::unique_ptr<Selection> f){
  std::vector<std::unique_ptr<Selection>> my_vec;
  my_vec.push_back(move(a)); my_vec.push_back(move(b)); my_vec.push_back(move(c)); my_vec.push_back(move(d)); my_vec.push_back(move(e)); my_vec.push_back(move(f));
  return my_vec;
}
inline std::vector<std::unique_ptr<Selection>> make_uvec(std::unique_ptr<Selection> a, std::unique_ptr<Selection> b, std::unique_ptr<Selection> c, std::unique_ptr<Selection> d, std::unique_ptr<Selection> e, std::unique_ptr<Selection> f, std::unique_ptr<Selection> g){
  std::vector<std::unique_ptr<Selection>> my_vec;
  my_vec.push_back(move(a)); my_vec.push_back(move(b)); my_vec.push_back(move(c)); my_vec.push_back(move(d)); my_vec.push_back(move(e)); my_vec.push_back(move(f)); my_vec.push_back(move(g));
  return my_vec;
}

inline std::vector<std::unique_ptr<AnalysisModule>> make_uvec(std::unique_ptr<AnalysisModule> a, std::unique_ptr<AnalysisModule> b, std::unique_ptr<AnalysisModule> c, std::unique_ptr<AnalysisModule> d){
  std::vector<std::unique_ptr<AnalysisModule>> my_vec;
  my_vec.push_back(move(a)); my_vec.push_back(move(b)); my_vec.push_back(move(c)); my_vec.push_back(move(d));
  return my_vec;
}


template<typename T>
int count_objects(const vector<T> & objects, const Event & event, const boost::optional<std::function<bool (const T &, const Event & )>> & object_id){
  int n_objects = objects.size();
  if(object_id){
    n_objects = 0;
    for(const auto & obj : objects){
      if((*object_id)(obj, event)) ++n_objects;
    }
  }  
  return n_objects;
} 
 
/*
//template<typename T>
vector<unique_ptr<Selection>> make_uvec(std::initializer_list<unique_ptr<Selection>> list){
  vector<unique_ptr<Selection>> my_vec;
  for (unsigned int i=0; i<list.size();i++)
    my_vec.push_back(move(list.at(i)));
  return my_vec;
}
*/


/*
template <class... Params>
void f(Params... params) {
    std::array<int, sizeof...(params)> list = {params...};
}

 std::initializer_list<int>::iterator it;  // same as: const int* it
    for ( it=args.begin(); it!=args.end(); ++it)


*/
/*
template<typename... Args>
  std::vector<Args> make_vec(Args... args){
  std::vector<Args> my_vec = {args...};
  return my_vec;
}

*/
