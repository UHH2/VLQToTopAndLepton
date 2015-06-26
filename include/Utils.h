#pragma once

#include "UHH2/core/include/Event.h"

template<typename P>
inline void sort_by_eta(std::vector<P> & particles){
    std::sort(particles.begin(), particles.end(), [](const P & p1, const P & p2){return p1.eta() > p2.eta();});
}

template<typename P>
inline void sort_by_eta(std::vector<P*> & particles){
    std::sort(particles.begin(), particles.end(), [](const P* p1, const P* p2){return p1->eta() > p2->eta();});
}
