#include<iostream>
#include<string>
#include<vector>
#include<cmath>
#include<limits>
#include<chrono>
#include<algorithm>
#include<numeric>
#include<fstream>

using Clock = std::chrono::high_resolution_clock;

namespace Const {
    constexpr double Joul = 4184;//Joul-Calorie conversion
    constexpr double Avogadro = 6.022E+23;
    constexpr double Boltzmann = 0.0138065;
    constexpr double pNnm = 1.0e21;
}
// =========================================Experiment conditions============================================
namespace Condition {
    constexpr double Temperature = 298.15;
    constexpr double kT = Temperature * Const::Boltzmann;
    constexpr int ArmLength = 2200;//total length of the 2 dsDNA arms, unit is base-pair.
    constexpr double PillarStiffness = 0.07406;//spring constant of the pillar/optical trap/micro-needle/etc that is used for stretching.
    constexpr double SaltConc = 100; //salt concentration in mM, this the the salt concentraion or the unzip experiment;
}
// ========================================DNA mechanical parameters=========================================
//These are pseudo constants..they are salt dependent and temperature dependent
namespace DNAParams{
    constexpr double LPDS = 51.97;//dsDNA persistence length
    constexpr double KDS = 1318;//dsDNA elastic modulus
    constexpr double L0DS = 0.338;//dsDNA contour length per bp
    constexpr double LPSS = 0.765;//ssDNA persistence length
    constexpr double KSS = 470;//ssDNA elastic modulus  
    constexpr double L0SS = 0.554;//ssDNA contour length per nt
}

#define VALID_MAX_FORCE 10000000.0
#define VALID_MIN_FORCE 0.2

// ====================================basepair energy measured by Huguet et al==============================
//ref: Huguet, Josep M., et al. (2010) PNAS
namespace BPEnergy {
    constexpr int LUTsize = 4;
    constexpr double LUTdH[LUTsize][LUTsize] = {
        {-7.28, -4.63, -5.21, -5.80},//aa, at, ag, ac
        {-8.31, -7.28, -8.96, -8.16},//ta, tt, tg, tc
        {-8.16, -5.80, -8.57, -10.1},//ga, gt, gg, gc
        {-8.96, -5.21, -9.66, -8.57} //ca, ct, cg, cc
    };
    constexpr double LUTdS[LUTsize][LUTsize] = {
        {-20.28, -11.62, -12.89, -14.46},//aa, at, ag, ac
        {-25.06, -20.28, -24.48, -22.46},//ta, tt, tg, tc
        {-22.46, -14.46, -22.30, -25.96},//ga, gt, gg, gc
        {-24.48, -12.89, -24.43, -22.30} //ca, ct, cg, cc
    };
    constexpr double LUTm[LUTsize][LUTsize] = {
        {0.145, 0.117, 0.070, 0.099},//aa, at, ag, ac
        {0.091, 0.145, 0.091, 0.155},//ta, tt, tg, tc
        {0.155, 0.099, 0.063, 0.079},//ga, gt, gg, gc
        {0.091, 0.070, 0.132, 0.063} //ca, ct, cg, cc
    };

    //Constexpr class
    class constexpr_map_class {
        public:
        constexpr int operator[] (char key) const {
                return bp2idx (key);
        }

        private:
        constexpr int bp2idx(char base) const {
            switch (base) {
            case 'a':
            case 'A':
                return 0;

            case 't':
            case 'T':
                return 1;

            case 'g':
            case 'G':
                return 2;

            case 'c':
            case 'C':
                return 3;

            default:
                return -1;
            }
        };
    } constexpr bp2idx_map; // a tip: can either put constexpr here, or before "class"

    static_assert (bp2idx_map['a'] == 0, "Error.");
    static_assert (bp2idx_map['A'] == 0, "Error.");
    static_assert (bp2idx_map['t'] == 1, "Error.");
    static_assert (bp2idx_map['T'] == 1, "Error.");
    static_assert (bp2idx_map['g'] == 2, "Error.");
    static_assert (bp2idx_map['G'] == 2, "Error.");
    static_assert (bp2idx_map['c'] == 3, "Error.");
    static_assert (bp2idx_map['C'] == 3, "Error.");
    //static_assert (bp2idx_map['C'] == 4, "Error.");//will cause error
}


//Calculate DNA sequence's energy
namespace DNAsequence {
    const double EffSaltConc = log(Condition::SaltConc * 0.001) / 298.0;//298.0 K is where the energy was measured in Huguet paper, it is not equal to Condition::temperature
    const double factor_pNnm = Const::pNnm * Const::Joul / Const::Avogadro / Condition::kT;//convert to 'pN-nm' unit;
    double lookup_bp_energy(char bp1, char bp2) {
        int idx1 = BPEnergy::bp2idx_map[bp1];
        int idx2 = BPEnergy::bp2idx_map[bp2];
        double energy = BPEnergy::LUTdH[idx1][idx2] - (BPEnergy::LUTm[idx1][idx2] * 
                        EffSaltConc + BPEnergy::LUTdS[idx1][idx2] * 0.001 ) * Condition::Temperature;
        return - energy * factor_pNnm;//convert to 'pN-nm' unit;
    }
    std::vector<double> calculate_sequence_energy(const std::string & sequence) {
        //from the DNA sequence, calculate the energy at every j_unzipped.

        std::vector<double> sequenceEnergy;
        if (sequence.size() < 2) {
            std::cerr << "Error: Sequence length must be greater than or equal to 2" << std::endl;
            return sequenceEnergy;
        }

        double accum = 0.0;
        *std::back_inserter(sequenceEnergy) = 0.0;
        std::transform(sequence.cbegin() + 1, sequence.cend(), std::back_inserter(sequenceEnergy), 
        [&accum](const char& bp2) {
            accum += lookup_bp_energy(*(&bp2 - 1), bp2);
            return accum;
            }
        );//equivalent to std::partial_sum()

        return sequenceEnergy;
    }
}

// these will become dynamic inputs in the real application
const char sequence[] = "TGGTTTACTCCTATACCGAGAAAAAACGTATTCGTAAGGATTTTGGTAAACGTCCACAAGTTCTGGATGTACCTTATCTCCTTTCTATCCAGCTTGACTCGTTTCAGAAATTTATCGAGCAAGATCCTGAAGGGCAGTATGGTCTGGAAGCTGCTTTCCGTTCCGTATTCCCGATTCAGAGCTACAGCGGTAATTCCGAGCTGCAATACGTCAGCTACCGCCTTGGCGAACCGGTGTTTGACGTCCAGGAATGTCAAATCCGTGGCGTGACCTATTCCGCACCGCTGCGCGTTAAACTGCGTCTGGTGATCTATGAGCGCGAAGCGCCGGAAGGCACCGTAAAAGACATTAAAGAACAAGAAGTCTACATGGGCGAAATTCCGCTCATGACAGACAACGGTACCTTTGTTATCAACGGTACTGAGCGTGTTATCGTTTCCCAGCTGCACCGTAGTCCGGGCGTCTTCTTTGACTCCGACAAAGGTAAAACCCACTCTTCGGGTAAAGTGCTGTATAACGCGCGTATCATCCCTTACCGTGGTTCCTGGCTGGACTTCGAATTCGATCCGAAGGACAACCTGTTCGTACGTATCGACCGTCGCCGTAAACTGCCTGCGACCATCATTCTGCGCGCCCTGAACTACACCACAGAGCAGATCCTCGACCTGTTCTTTGAAAAAGTTATCTTTGAAATCCGTGATAACAAGCTGCAGATGGAACTGGTGCCGGAACGCCTGCGTGGTGAAACCGCATCTTTTGACATCGAAGCTAACGGTAAAGTGTACGTAGAAAAAGGCCGCCGTATCACTGCGCGCCACATTCGCCAGCTGGAAAAAGACGACGTCAAACTGATCGAAGTCCCGGTTGAGTACATCGCAGGTAAAGTGGTTGCTAAAGACTATATTGATGAGTCTACCGGCGAGCTGATCTGCGCAGCGAACATGGAGCTGAGCCTGGATCTGCTGGCTAAGCTGAGCCAGTCTGGTCACAAGCGTATCGAAACGCTGTTCACCAACGATCTGGATCACGGCCCATATATCTCTGAAACCTTACGTGTCGACCCAACTAACGACCGTCTGAGCGCACTGGTAGAAATCTACCGCATGATGCGCCCTGGCGAGCCGCCGACTCGTGAAGCAGCTGAAAGCCTGTTCGAGAACCTGTTCTTCTCCGAAGACCGTTATGACTTGTCTGCGGTTGGTCGTATGAAGTTCAACCGTTCTCTGCTGCGCGAAGAAATCGAAGGTTCCGGTATCCTGAGCAAAGACGACATCATTGATGTTATGAAAAAGCTCATCGATATCCGTAACGGTAAAGGCGAAGTCGATGATATCGACCACCTCGGCAACCGTCGTATCCGTTCCGTTGGCGAAATGGCGGAAAACCAGTTCCGCGTTGGCCTGGTACGTGTAGAGCGTGCGGTGAAAGAGCGTCTGTCTCTGGGCGATCTGGATACCCTGATGCCACAGGATATGATCAACGCCAAGCCGATTTCCGCAGCAGTGAAAGAGTTCTTCGGTTCCAGCCAGCTGTCTCAGTTTATGGACCAGAACAACCCGCTGTCTGAGATTACGCACAAACGTCGTATCTCCGCACTCGGCCCAGGCGGTCTGACCCGTGAACGTGCAGGCTTCGAAGTTCGAGACGTACACCCGACTCACTACGGTCGCGTATGTCCAATCGAAACCCCTGAAGGTCCGAACATCGGTCTGATCAACTCTCTGTCCGTGTACGCACAGACTAACGAATACGGCTTCCTTGAGACTCCGTATCGTAAAGTGACCGACGGTGTTGTAACTGACGAAATTCACTACCTGTCTGCTATCGAAGAAGGCAACTACGTTATCGCCCAGGCGAACTCCAACTTGGATGAAGAAGGCCACTTCGTAGAAGACCTGGTAACTTGCCGTAGCAAAGGCGAATCCAGCTTGTTCAGCCGCGACCAGGTTGACTACATGGACGTATCCACCCAGCAGGTGGTATCCGTCGGTGCGTCCCTGATCCCGTTCCTGGAACACGATGACGCCAACCGTGCATTGATGGGTGCGAACATGCAACGTCAGGCCGTTCCGACTCTGCGCGCTGATAAGCCGCTGGTTGGTACTGGTATGGAACGTGCTGTTGCCGTTGACTCCGGTGTAACTGCGGTAGCTAAACGTGGTGGTGTCGTTCAGTACGTGGATGCTTCCCGTATCGTTATCAAAGTTAACGAAGACGAGATGTATCCGGGTGAAGCAGGTATCGACATCTACAACCTGACCAAATACACCCGTTCTAACCAGAACACCTGTATCAACCAGATGCCGTGTGTGTCTCTGGGTGAACCGGTTGAACGTGGCGACGTGCTGGCAGACGGTCCGTCCACCGACCTCGGTGAACTGGCGCTTGGTCAGAACATGCGCGTAGCGTTCATGCCGTGGAATGGTTACAACTTCGAAGACTCCATCCTCGTATCCGAGCGTGTTGTTCAGGAAGACCGTTTCACCACCATCCACATTCAGGAACTGGCGTGTGTGTCCCGTGACACCAAGCTGGGTCCGGAAGAGATCACCGCTGACATCCCGAACGTGGGTGAAGCTGCGCTCTCCAAACTGGATGAATCCGGTATCGTTTACATTGGTGCGGAAGTGACCGGTGGCGACATTCTGGTTGGTAAGGTAACGCCGAAAGGTGAAACTCAGCTGACCCCAGAAGAAAAACTGCTGCGTGCGATCTTCGGTGAGAAAGCCTCTGACGTTAAAGACTCTTCTCTGCGCGTACCAAACGGTGTATCCGGTACGGTTATCGACGTTCAGGTCTTTACTCGCGATGGCGTAGAAAAAGACAAACGTGCGCTGGAAATCGAAGAAATGCAGCTCAAACAGGCGAAGAAAGACCTGTCTGAAGAACTGCAGATCCTCGAAGCGGGTCTGTTCAGCCGTATCCGTGCTGTGCTGGTAGCCGGTGGCGTTGAAGCTGAGAAGCTCGACAAACTGCCGCGCGATCGCTGGCTGGAGCTGGGCCTGACAGACGAAGAGAAACAAAATCAGCTGGAACAGCTGGCTGAGCAGTATGACGAACTGAAACACGAGTTCGAGAAGAAACTCGAAGCGAAACGCCGCAAAATCACCCAGGGCGACGATCTGGCACCGGGCGTGCTGAAGATTGTTAAGGTATATCTGGCGGTTAAACGCCGTATCCAGCCTGGTGACAAGATGGCAGGTCGTCACGGTAACAAGGGTGTAATTTCTAAGATCAACCCGATCGAAGATATGCCTTACGATGAAAACGGTACGCCGGTAGACATCGTACTGAACCCGCTGGGCGTACCGTCTCGTATGAACATCGGTCAGATCCTCGAAACCCACCTGGGTATGGCTGCGAAAGGTATCGGCGACAAGATCAACGCCATGCTGAAACAGCAGCAAGAAGTCGCGAAACTGCGCGAATTCATCCAGCGTGCGTACGATCTGGGCGCTGACGTTCGTCAGAAAGTTGACCTGAGTACCTTCAGCGATGAAGAAGTTATGCGTCTGGCTGAAAACCTGCGCAAAGGTATGCCAATCGCAACGCCGGTGTTCGACGGTGCGAAAGAAGCAGAAATTAAAGAGCTGCTGAAACTTGGCGACCTGCCGACTTCCGGTCAGATCCGCCTGTACGATGGTCGCACTGGTGAACAGTTCGAGCGTCCGGTAACCGTTGGTTACATGTACATGCTGAAACTGAACCACCTGGTCGACGACAAGATGCACGCGCGTTCCACCGGTTCTTACAGCCTGGTTACTCAGCAGCCGCTGGGTGGTAAGGCACAGTTCGGTGGTCAGCGTTTCGGGGAGATGGAAGTGTGGGCGCTGGAAGCATACGGCGCAGCATACACCCTGCAGGAAATGCTCACCGTTAAGTCTGATGACGTGAACGGTCGTACCAAGATGTATAAAAACATCGTGGACGGCAACCATCAGATGGAGCCGGGCATGCCAGAATCCTTCAACGTATTGTTGAAAGAGATTCGTTCGCTGGGTATCAACATCGAACTGGAAGACGAGTAA";
const std::vector<double> seq_energy = DNAsequence::calculate_sequence_energy(sequence);

namespace MyMath{
    __inline__ __device__ double Langevin(double x) {
        return 1.0/tanh(x) - 1.0/x;
    }
    __inline__ __device__ double Langevin_integ(double x) {
        return log(sinh(x)/x);
    }
}

// ==========================================DNA mechanical models===========================================
namespace DNAModel {
    //================================WLC/FJC model: parameter definitions==================================
    //phi = x/L (L = contour length)
    //alpha = fA/kT (A = persistence length)
    //k0_eff=k0A/kT (K0 = elastic modulus)

    //====================================Marko-Siggia 1995 WLC========================================
    __inline__ __device__ double phi2alpha_MS(double phi){ 
        return phi + 0.25 / ((1.0 - phi) * (1.0 - phi)) - 0.25; 
    }


    //==============================WLC high force, Odijk 1995 macromolecules===========================
    __inline__ __device__ double alpha2phi_Odijk95(double alpha, double k0_eff) { 
        // if (alpha < 0.25) {
        //     return 0.0;// undefined at alpha == 0, just give it a small value
        //     //I can do this because this is force, energy must be calculated correctly!!
        // }
        return 1.0 - 0.5 / sqrt(alpha) + alpha / k0_eff; 
    }
    __inline__ __device__ double integ_phidalpha_Odijk95(double alpha, double k0_eff) { 
        return alpha - sqrt(alpha) + 0.5 * alpha * alpha / k0_eff; 
    }
    __inline__ __device__ double integ_alphadphi_Odijk95(double alpha, double k0_eff) {
        return alpha * alpha2phi_Odijk95(alpha, k0_eff) - integ_phidalpha_Odijk95(alpha, k0_eff);
    }
    // ================================MODIFIED VERSION OF FJC, Smith 1995 macromolecules================================
    //Modified version specific for ssDNA force region, and keeps accuracy
    //For ssDNA, alpha = (force * lp_ss / kT) = force /5.4, a force range of (0.1 ~ 60) is alpha < 12
    //My homemade Langevin_integ function should be accurate enough in this region.
    __inline__ __device__ double alpha2phi_Smith95_m(double alpha, double k0_eff) {//"m" means modified
        return MyMath::Langevin(2.0 * alpha) + alpha / k0_eff;
    }
    __inline__ __device__ double integ_phidalpha_Smith95_m(double alpha, double k0_eff) { 
        return 0.5 * MyMath::Langevin_integ(2.0 * alpha) + 0.5 * alpha * alpha / k0_eff;
    }
    __inline__ __device__ double integ_alphadphi_Smith95_m(double alpha, double k0_eff) {
        //integ actually starts from 1, but it's OK since it is for partition function calculation
        return alpha * alpha2phi_Smith95_m(alpha, k0_eff) - integ_phidalpha_Smith95_m(alpha, k0_eff);
    }
}
__inline__ __device__ double lz_ds (double force) {//dsDNA's length per base
    return Condition::ArmLength * DNAParams::L0DS * 
            DNAModel::alpha2phi_Odijk95(force * DNAParams::LPDS / Condition::kT, DNAParams::KDS * DNAParams::LPDS / Condition::kT);
}
__inline__ __device__ double lz_ss (double force, int j) {//ssDNA's length per base
    return 2.0 * j * DNAParams::L0SS * 
            DNAModel::alpha2phi_Smith95_m(force * DNAParams::LPSS / Condition::kT, DNAParams::KSS * DNAParams::LPSS / Condition::kT);
}
__inline__ __device__ double le_ds (double force) {//function version of dsDNA's energy per bp:
    return Condition::ArmLength * Condition::kT * DNAParams::L0DS * 
            DNAModel::integ_alphadphi_Odijk95(force * DNAParams::LPDS / Condition::kT, DNAParams::KDS * DNAParams::LPDS / Condition::kT) / DNAParams::LPDS;
}
__inline__ __device__ double le_ss (double force, int j) {//function version of ssDNA's energy per bp:
    return 2.0 * j * Condition::kT * DNAParams::L0SS * 
            DNAModel::integ_alphadphi_Smith95_m(force * DNAParams::LPSS / Condition::kT, DNAParams::KSS * DNAParams::LPSS / Condition::kT) / DNAParams::LPSS;
}
__inline__ __device__ double delta_ext(double force, double j, double ext) {//func used to find force so the system total extension = ext
    return ext - force/Condition::PillarStiffness - lz_ds (force) - lz_ss (force, j);//increasing function with force
}

__inline__ __device__ double find_force(int j, double ext) {
    double f1 = ValidRange::ValidMinForce;
    double f2 = ValidRange::VALID_MAX_FORCE;

    double y1 = delta_ext(f1, j, ext);
    double y2 = delta_ext(f2, j, ext);

    double fm = 0.0;
    double ym = 0.0;
    for (int cnt = 0; cnt < 50; ++cnt) {
        fm = (f1 + f2) * 0.5;
        ym = delta_ext(fm, j, ext);
        f1 = (y1 < y2 && ym < 0 || y1 > y2 && ym > 0) ? fm : f1;
        f2 = (y1 < y2 && ym > 0 || y1 > y2 && ym < 0) ? fm : f2;

        y1 = (y1 < y2 && ym < 0 || y1 > y2 && ym > 0) ? ym : y1;
        y2 = (y1 < y2 && ym > 0 || y1 > y2 && ym < 0) ? ym : y2;
    }
    return fm;
}

struct dp {// a data point
    int extension_total = 0;//in nm;
    double extension_DNA = 0.0;//in nm;
    double force_average = 0.0;//in pN
    double force_SD = 0.0;//in pN
    double junzipped_average = 0.0;//#bp unzipped
    double junzipped_SD = 0.0;//#bp unzipped
};

constexpr double energy_threshold = 50.0;//don't calculate exp(-p/kT) and set probability to 0;


__global__ dp calculate_array(int extension, const std::vector<double> & seq_energy) {
        
    std::vector<double> temp_e(seq_energy.size(), 0.0);
    std::vector<double> temp_f(seq_energy.size(), 0.0);

    double min_e = 1.0e100;
    for (int j = 0; j < seq_energy.size(); ++j) {
        double f = find_force(j, extension);
        double energy = 0.0;
        if (f >= VALID_MAX_FORCE) {
            energy = std::nan("1");//make this a large number, meaning that do not use the value
        } else if (f <= VALID_MIN_FORCE) {
            energy = 0; // these DNA mechanical models are invalid at such low forces, just set energy to zero
        } else {
            energy = (0.5 * f * f / Condition::PillarStiffness + le_ds(f) + le_ss(f, j))/Condition::kT;
        }
        
        temp_f.at(j) = f;
        temp_e.at(j) = energy + seq_energy[j];

        if (min_e > temp_e.at(j)) {
            min_e = temp_e.at(j);
        }
    }

    double prob = 0;
    double Fprob = 0;
    double FFprob = 0;
    double Jprob = 0;
    double JJprob = 0;
    double p,f;
    long long mean_iter = 0;
    int max_iter = 0;
    int min_iter = 10000000;
    for (int j = 0; j < seq_energy.size(); ++j) {

        p = temp_e.at(j) - min_e;
        p = p > energy_threshold ? 0.0 : std::exp(-p);
        f = temp_f.at(j);

        prob += p;
        Fprob += f * p;
        FFprob += f * f * p;
        Jprob += j * p;
        JJprob += j * j * p;
    }

    dp point;
    point.extension_total = extension;
    point.force_average = Fprob/prob;
    point.extension_DNA = extension - point.force_average/Condition::PillarStiffness;
    point.force_SD = std::sqrt(FFprob/prob -  (Fprob/prob) * (Fprob/prob));
    point.junzipped_average = Jprob/prob;
    point.junzipped_SD = std::sqrt(JJprob/prob -  (Jprob/prob) * (Jprob/prob));
    return point;
}

int main(int argc, char* argv[]) {
    constexpr int ext_start = 0;//in nm
    constexpr int ext_end = 1500;//in nm
    constexpr int ext_step = 1;//in nm

    std::vector<dp> results;

    auto start_time = Clock::now();

    for (int ext = ext_start; ext <= static_cast<int>(1.2 * seq_energy.size()); ext += ext_step) {
        results.push_back(calculate_array(ext, seq_energy));
        if (ext % 100 == 0) {
            std::cout << "Calculated extension: " << ext << " nm" s<< std::endl;
        }
    }

    auto end_time = Clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << "Computation Time: " << elapsed.count() << " seconds" << std::endl;

    //Output results
    if (argc > 1) { // argc > 1 == verbose mode
        std::cout << "Extension_total(nm),Extension_DNA(nm),Force_average(pN),Force_SD(pN),Junzipped_average(bp),Junzipped_SD(bp)" << std::endl;
        for (const auto& point : results) {
            std::cout << point.extension_total << ","
                    << point.extension_DNA << ","
                    << point.force_average << ","
                    << point.force_SD << ","
                    << point.junzipped_average << ","
                    << point.junzipped_SD << std::endl;
        }
    }
    // output results to a CSV file
    std::ofstream outfile("unzipping_results.csv");
    outfile << "Extension_total(nm),Extension_DNA(nm),Force_average(pN),Force_SD(pN),Junzipped_average(bp),Junzipped_SD(bp)" << std::endl;
    for (const auto& point : results) {
        outfile << point.extension_total << ","
                << point.extension_DNA << ","
                << point.force_average << ","
                << point.force_SD << ","
                << point.junzipped_average << ","
                << point.junzipped_SD << std::endl;
    }
    outfile.close();

    return 0;
}
