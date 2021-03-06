#pragma once

#include "UHH2/core/include/Hists.h"

/**  \brief Example class for booking and filling histograms
 * 
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */

class VLQToTopAndLeptonHists: public uhh2::Hists {
 public:
    // use the same constructor arguments as Hists for forwarding:
    VLQToTopAndLeptonHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~VLQToTopAndLeptonHists();
 private:
    TH1F* N_jets, *eta_jet1, *eta_jet2, *eta_jet3, *eta_jet4;
    TH2F* pt_eta_jet1, *pt_eta_jet2, *pt_eta_jet3, *pt_eta_jet4;
    TH1F* N_mu, *pt_mu, *eta_mu, *reliso_mu;
    TH1F* N_pv;
};
