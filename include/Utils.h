#pragma once

#include "UHH2/core/include/Selection.h"

#include <initializer_list>
#include "UHH2/core/include/Event.h"

template<typename P>
inline void sort_by_eta(std::vector<P> & particles){
    std::sort(particles.begin(), particles.end(), [](const P & p1, const P & p2){return p1.eta() > p2.eta();});
}

template<typename P>
inline void sort_by_eta(std::vector<P*> & particles){
    std::sort(particles.begin(), particles.end(), [](const P* p1, const P* p2){return p1->eta() > p2->eta();});
}



//template<typename T>
vector<unique_ptr<Selection>> make_uvec(std::initializer_list<unique_ptr<Selection>> list){
  vector<unique_ptr<Selection>> my_vec;
  for (unsigned int i=0; i<list.size();i++)
    my_vec.push_back(move(list.at(i)));
  return my_vec;
}



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
