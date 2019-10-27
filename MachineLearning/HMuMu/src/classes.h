#include <vector>
#include <TVector3.h>
namespace MachineLearning_HMuMu {
    struct dictionary {
        std::vector<TVector3> v_v3;
        std::vector<char> v_char;
        std::vector<std::vector<int> > vv_int;
        std::vector<std::vector<TVector3> > vv_v3;
    };
}
