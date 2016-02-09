#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

template<typename T>
class JetPtAndMultFixerWeight: public uhh2::AnalysisModule {
public:
    explicit JetPtAndMultFixerWeight(uhh2::Context & ctx,
                            string const & h_in,
                            float offset, float gradient,
                            string const & h_weight = "weight_jetpt",
                            bool apply_event_weight = false) :
        h_in_(ctx.get_handle<std::vector<T>>(h_in)),
        offset_(offset), gradient_(gradient),
        cov_p0_p0_(0.), cov_p0_p1_(0.), cov_p1_p1_(0.),
        h_weight_(ctx.declare_event_output<float>(h_weight)),
        h_weight_up_(ctx.declare_event_output<float>(h_weight+"_up")),
        h_weight_down_(ctx.declare_event_output<float>(h_weight+"_down")),
        apply_event_weight_(apply_event_weight) {
            auto dataset_type = ctx.get("dataset_type");
            is_mc = dataset_type == "MC";
            if (!is_mc) {
                cout << "Warning: JetPtAndMultFixerWeight will not have an effect on "
                <<" this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
                return;
            }
        }

    explicit JetPtAndMultFixerWeight(uhh2::Context & ctx,
                            string const & h_in,
                            float offset, float gradient,
                            float cov_p0_p0, float cov_p0_p1, float cov_p1_p1,
                            string const & h_weight = "weight_jetpt",
                            bool apply_event_weight = false) :
        h_in_(ctx.get_handle<std::vector<T>>(h_in)),
        offset_(offset), gradient_(gradient),
        cov_p0_p0_(cov_p0_p0), cov_p0_p1_(cov_p0_p1), cov_p1_p1_(cov_p1_p1),
        h_weight_(ctx.declare_event_output<float>(h_weight)),
        h_weight_up_(ctx.declare_event_output<float>(h_weight+"_up")),
        h_weight_down_(ctx.declare_event_output<float>(h_weight+"_down")),
        apply_event_weight_(apply_event_weight) {
            auto dataset_type = ctx.get("dataset_type");
            is_mc = dataset_type == "MC";
            if (!is_mc) {
                cout << "Warning: JetPtAndMultFixerWeight will not have an effect on "
                <<" this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
                return;
            }
        }

    virtual bool process(uhh2::Event & event) override {
        if (!event.is_valid(h_in_))
            return false;
        auto const & coll = event.get(h_in_);
        float weight = 1.0f;
        float weight_up = 1.0f;
        float weight_down = 1.0f;
        if (is_mc) {
	  weight = (0.5-coll.size()*0.1) + offset_;
	  /*
            for (auto const & part : coll) {
                float part_pt = part.pt();
                float sf = offset_ + part_pt * gradient_;
                float sf_err = std::sqrt(cov_p0_p0_ + 2 * part_pt * cov_p0_p1_ + part_pt * part_pt * cov_p1_p1_);

                weight *= std::min(1.0f, sf);
                weight_up *= std::min(1.0f, sf + sf_err);
                weight_down *= std::min(1.0f, sf - sf_err);
            }
	  */
	    //if(coll.size()<5)std::cout<<"event weight factor "<<weight<<" # Ak4Jets "<<coll.size()<<std::endl; 

	    if (apply_event_weight_) {
	      event.weight *= weight;
	    }


	}
        event.set(h_weight_, weight);
        event.set(h_weight_up_, weight_up);
        event.set(h_weight_down_, weight_down);
        // event.weight *= weight;
        return true;
    }

private:
    uhh2::Event::Handle<std::vector<T>> h_in_;
    float offset_, gradient_;
    float cov_p0_p0_, cov_p0_p1_, cov_p1_p1_;
    uhh2::Event::Handle<float> h_weight_, h_weight_up_, h_weight_down_;
    bool is_mc, apply_event_weight_;
    

};
