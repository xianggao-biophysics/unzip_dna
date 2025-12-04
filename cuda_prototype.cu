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

namespace Constants {
    // 1. Physical constants
    constexpr float Joule = 4184;//Joule-Calorie conversion
    constexpr float Avogadro = 6.022E+23;
    constexpr float Boltzmann = 0.0138065;
    constexpr float pNnm = 1.0e21;

    // 2. Experiment conditions:
    constexpr float Temperature = 298.15;
    constexpr float kT = Temperature * Boltzmann;
    constexpr float SaltConc = 100; //salt concentration in mM, this the the salt concentraion or the unzip experiment;

    // 3. Mechanical parameters 
    //DNA's parameters are pseudo constants -- they are salt dependent and temperature dependent
    __constant__ int ArmLength = 2200;//total length of the 2 dsDNA arms, unit is base-pair.
    __constant__ float TrapStiffness = 0.07406;//spring constant of the pillar/optical trap/micro-needle/etc that is used for stretching.
    
    __constant__ float LPDS = 51.97;//dsDNA persistence length
    __constant__ float KDS = 1318;//dsDNA elastic modulus
    __constant__ float L0DS = 0.338;//dsDNA contour length per bp
    __constant__ float LPSS = 0.765;//ssDNA persistence length
    __constant__ float KSS = 470;//ssDNA elastic modulus  
    __constant__ float L0SS = 0.554;//ssDNA contour length per nt
}

#define VALID_MAX_FORCE 10000000.0
#define VALID_MIN_FORCE 0.2

// convertion from sequence to energy can be done on CPU side, since it is only done once
// ====================================basepair energy measured by Huguet et al==============================
//ref: Huguet, Josep M., et al. (2010) PNAS
namespace BPEnergy {
    constexpr int LUTsize = 4;
    constexpr float LUTdH[LUTsize][LUTsize] = {
        {-7.28, -4.63, -5.21, -5.80},//aa, at, ag, ac
        {-8.31, -7.28, -8.96, -8.16},//ta, tt, tg, tc
        {-8.16, -5.80, -8.57, -10.1},//ga, gt, gg, gc
        {-8.96, -5.21, -9.66, -8.57} //ca, ct, cg, cc
    };
    constexpr float LUTdS[LUTsize][LUTsize] = {
        {-20.28, -11.62, -12.89, -14.46},//aa, at, ag, ac
        {-25.06, -20.28, -24.48, -22.46},//ta, tt, tg, tc
        {-22.46, -14.46, -22.30, -25.96},//ga, gt, gg, gc
        {-24.48, -12.89, -24.43, -22.30} //ca, ct, cg, cc
    };
    constexpr float LUTm[LUTsize][LUTsize] = {
        {0.145, 0.117, 0.070, 0.099},//aa, at, ag, ac
        {0.091, 0.145, 0.091, 0.155},//ta, tt, tg, tc
        {0.155, 0.099, 0.063, 0.079},//ga, gt, gg, gc
        {0.091, 0.070, 0.132, 0.063} //ca, ct, cg, cc
    };

    int bp2idx(char base) {
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
    
    float EffSaltConc = log(Constants::SaltConc * 0.001) / 298.0;//298.0 K is where the energy was measured in Huguet paper, it is not equal to Condition::temperature
    float factor_pNnm = Constants::pNnm * Constants::Joule / Constants::Avogadro / Constants::kT;//convert to 'pN-nm' unit;
    float lookup_bp_energy(char bp1, char bp2) {
        int idx1 = bp2idx(bp1);
        int idx2 = bp2idx(bp2);
        float energy = LUTdH[idx1][idx2] - (LUTm[idx1][idx2] * 
                        EffSaltConc + LUTdS[idx1][idx2] * 0.001 ) * Constants::Temperature;
        return - energy * factor_pNnm;//convert to 'pN-nm' unit;
    }
    std::vector<float> calculate_sequence_energy(const std::string & sequence) {
        //from the DNA sequence, calculate the energy at every j_unzipped.

        std::vector<float> sequenceEnergy;
        if (sequence.size() < 2) {
            std::cerr << "Error: Sequence length must be greater than or equal to 2" << std::endl;
            return sequenceEnergy;
        }

        float accum = 0.0;
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
constexpr char sequence[] = "TGGTTTACTCCTATACCGAGAAAAAACGTATTCGTAAGGATTTTGGTAAACGTCCACAAGTTCTGGATGTACCTTATCTCCTTTCTATCCAGCTTGACTCGTTTCAGAAATTTATCGAGCAAGATCCTGAAGGGCAGTATGGTCTGGAAGCTGCTTTCCGTTCCGTATTCCCGATTCAGAGCTACAGCGGTAATTCCGAGCTGCAATACGTCAGCTACCGCCTTGGCGAACCGGTGTTTGACGTCCAGGAATGTCAAATCCGTGGCGTGACCTATTCCGCACCGCTGCGCGTTAAACTGCGTCTGGTGATCTATGAGCGCGAAGCGCCGGAAGGCACCGTAAAAGACATTAAAGAACAAGAAGTCTACATGGGCGAAATTCCGCTCATGACAGACAACGGTACCTTTGTTATCAACGGTACTGAGCGTGTTATCGTTTCCCAGCTGCACCGTAGTCCGGGCGTCTTCTTTGACTCCGACAAAGGTAAAACCCACTCTTCGGGTAAAGTGCTGTATAACGCGCGTATCATCCCTTACCGTGGTTCCTGGCTGGACTTCGAATTCGATCCGAAGGACAACCTGTTCGTACGTATCGACCGTCGCCGTAAACTGCCTGCGACCATCATTCTGCGCGCCCTGAACTACACCACAGAGCAGATCCTCGACCTGTTCTTTGAAAAAGTTATCTTTGAAATCCGTGATAACAAGCTGCAGATGGAACTGGTGCCGGAACGCCTGCGTGGTGAAACCGCATCTTTTGACATCGAAGCTAACGGTAAAGTGTACGTAGAAAAAGGCCGCCGTATCACTGCGCGCCACATTCGCCAGCTGGAAAAAGACGACGTCAAACTGATCGAAGTCCCGGTTGAGTACATCGCAGGTAAAGTGGTTGCTAAAGACTATATTGATGAGTCTACCGGCGAGCTGATCTGCGCAGCGAACATGGAGCTGAGCCTGGATCTGCTGGCTAAGCTGAGCCAGTCTGGTCACAAGCGTATCGAAACGCTGTTCACCAACGATCTGGATCACGGCCCATATATCTCTGAAACCTTACGTGTCGACCCAACTAACGACCGTCTGAGCGCACTGGTAGAAATCTACCGCATGATGCGCCCTGGCGAGCCGCCGACTCGTGAAGCAGCTGAAAGCCTGTTCGAGAACCTGTTCTTCTCCGAAGACCGTTATGACTTGTCTGCGGTTGGTCGTATGAAGTTCAACCGTTCTCTGCTGCGCGAAGAAATCGAAGGTTCCGGTATCCTGAGCAAAGACGACATCATTGATGTTATGAAAAAGCTCATCGATATCCGTAACGGTAAAGGCGAAGTCGATGATATCGACCACCTCGGCAACCGTCGTATCCGTTCCGTTGGCGAAATGGCGGAAAACCAGTTCCGCGTTGGCCTGGTACGTGTAGAGCGTGCGGTGAAAGAGCGTCTGTCTCTGGGCGATCTGGATACCCTGATGCCACAGGATATGATCAACGCCAAGCCGATTTCCGCAGCAGTGAAAGAGTTCTTCGGTTCCAGCCAGCTGTCTCAGTTTATGGACCAGAACAACCCGCTGTCTGAGATTACGCACAAACGTCGTATCTCCGCACTCGGCCCAGGCGGTCTGACCCGTGAACGTGCAGGCTTCGAAGTTCGAGACGTACACCCGACTCACTACGGTCGCGTATGTCCAATCGAAACCCCTGAAGGTCCGAACATCGGTCTGATCAACTCTCTGTCCGTGTACGCACAGACTAACGAATACGGCTTCCTTGAGACTCCGTATCGTAAAGTGACCGACGGTGTTGTAACTGACGAAATTCACTACCTGTCTGCTATCGAAGAAGGCAACTACGTTATCGCCCAGGCGAACTCCAACTTGGATGAAGAAGGCCACTTCGTAGAAGACCTGGTAACTTGCCGTAGCAAAGGCGAATCCAGCTTGTTCAGCCGCGACCAGGTTGACTACATGGACGTATCCACCCAGCAGGTGGTATCCGTCGGTGCGTCCCTGATCCCGTTCCTGGAACACGATGACGCCAACCGTGCATTGATGGGTGCGAACATGCAACGTCAGGCCGTTCCGACTCTGCGCGCTGATAAGCCGCTGGTTGGTACTGGTATGGAACGTGCTGTTGCCGTTGACTCCGGTGTAACTGCGGTAGCTAAACGTGGTGGTGTCGTTCAGTACGTGGATGCTTCCCGTATCGTTATCAAAGTTAACGAAGACGAGATGTATCCGGGTGAAGCAGGTATCGACATCTACAACCTGACCAAATACACCCGTTCTAACCAGAACACCTGTATCAACCAGATGCCGTGTGTGTCTCTGGGTGAACCGGTTGAACGTGGCGACGTGCTGGCAGACGGTCCGTCCACCGACCTCGGTGAACTGGCGCTTGGTCAGAACATGCGCGTAGCGTTCATGCCGTGGAATGGTTACAACTTCGAAGACTCCATCCTCGTATCCGAGCGTGTTGTTCAGGAAGACCGTTTCACCACCATCCACATTCAGGAACTGGCGTGTGTGTCCCGTGACACCAAGCTGGGTCCGGAAGAGATCACCGCTGACATCCCGAACGTGGGTGAAGCTGCGCTCTCCAAACTGGATGAATCCGGTATCGTTTACATTGGTGCGGAAGTGACCGGTGGCGACATTCTGGTTGGTAAGGTAACGCCGAAAGGTGAAACTCAGCTGACCCCAGAAGAAAAACTGCTGCGTGCGATCTTCGGTGAGAAAGCCTCTGACGTTAAAGACTCTTCTCTGCGCGTACCAAACGGTGTATCCGGTACGGTTATCGACGTTCAGGTCTTTACTCGCGATGGCGTAGAAAAAGACAAACGTGCGCTGGAAATCGAAGAAATGCAGCTCAAACAGGCGAAGAAAGACCTGTCTGAAGAACTGCAGATCCTCGAAGCGGGTCTGTTCAGCCGTATCCGTGCTGTGCTGGTAGCCGGTGGCGTTGAAGCTGAGAAGCTCGACAAACTGCCGCGCGATCGCTGGCTGGAGCTGGGCCTGACAGACGAAGAGAAACAAAATCAGCTGGAACAGCTGGCTGAGCAGTATGACGAACTGAAACACGAGTTCGAGAAGAAACTCGAAGCGAAACGCCGCAAAATCACCCAGGGCGACGATCTGGCACCGGGCGTGCTGAAGATTGTTAAGGTATATCTGGCGGTTAAACGCCGTATCCAGCCTGGTGACAAGATGGCAGGTCGTCACGGTAACAAGGGTGTAATTTCTAAGATCAACCCGATCGAAGATATGCCTTACGATGAAAACGGTACGCCGGTAGACATCGTACTGAACCCGCTGGGCGTACCGTCTCGTATGAACATCGGTCAGATCCTCGAAACCCACCTGGGTATGGCTGCGAAAGGTATCGGCGACAAGATCAACGCCATGCTGAAACAGCAGCAAGAAGTCGCGAAACTGCGCGAATTCATCCAGCGTGCGTACGATCTGGGCGCTGACGTTCGTCAGAAAGTTGACCTGAGTACCTTCAGCGATGAAGAAGTTATGCGTCTGGCTGAAAACCTGCGCAAAGGTATGCCAATCGCAACGCCGGTGTTCGACGGTGCGAAAGAAGCAGAAATTAAAGAGCTGCTGAAACTTGGCGACCTGCCGACTTCCGGTCAGATCCGCCTGTACGATGGTCGCACTGGTGAACAGTTCGAGCGTCCGGTAACCGTTGGTTACATGTACATGCTGAAACTGAACCACCTGGTCGACGACAAGATGCACGCGCGTTCCACCGGTTCTTACAGCCTGGTTACTCAGCAGCCGCTGGGTGGTAAGGCACAGTTCGGTGGTCAGCGTTTCGGGGAGATGGAAGTGTGGGCGCTGGAAGCATACGGCGCAGCATACACCCTGCAGGAAATGCTCACCGTTAAGTCTGATGACGTGAACGGTCGTACCAAGATGTATAAAAACATCGTGGACGGCAACCATCAGATGGAGCCGGGCATGCCAGAATCCTTCAACGTATTGTTGAAAGAGATTCGTTCGCTGGGTATCAACATCGAACTGGAAGACGAGTAA";
const std::vector<float> seq_energy = BPEnergy::calculate_sequence_energy(sequence);

inline __device__ float Langevin(float x) {
    return 1.0/tanh(x) - 1.0/x;
}
inline __device__ float Langevin_integ(float x) {
    return log(sinh(x)/x);
}

// ==========================================DNA mechanical models===========================================
namespace DNAModel {
    //================================WLC/FJC model: parameter definitions==================================
    //phi = x/L (L = contour length)
    //alpha = fA/kT (A = persistence length)
    //k0_eff=k0A/kT (K0 = elastic modulus)

    //====================================Marko-Siggia 1995 WLC========================================
    inline __device__ float phi2alpha_MS(float phi){ 
        return phi + 0.25 / ((1.0 - phi) * (1.0 - phi)) - 0.25; 
    }


    //==============================WLC high force, Odijk 1995 macromolecules===========================
    inline __device__ float alpha2phi_Odijk95(float alpha, float k0_eff) { 
        // if (alpha < 0.25) {
        //     return 0.0;// undefined at alpha == 0, just give it a small value
        //     //I can do this because this is force, energy must be calculated correctly!!
        // }
        return 1.0 - 0.5 / sqrt(alpha) + alpha / k0_eff; 
    }
    inline __device__ float integ_phidalpha_Odijk95(float alpha, float k0_eff) { 
        return alpha - sqrt(alpha) + 0.5 * alpha * alpha / k0_eff; 
    }
    inline __device__ float integ_alphadphi_Odijk95(float alpha, float k0_eff) {
        return alpha * alpha2phi_Odijk95(alpha, k0_eff) - integ_phidalpha_Odijk95(alpha, k0_eff);
    }
    // ================================MODIFIED VERSION OF FJC, Smith 1995 macromolecules================================
    //Modified version specific for ssDNA force region, and keeps accuracy
    //For ssDNA, alpha = (force * lp_ss / kT) = force /5.4, a force range of (0.1 ~ 60) is alpha < 12
    //My homemade Langevin_integ function should be accurate enough in this region.
    inline __device__ float alpha2phi_Smith95_m(float alpha, float k0_eff) {//"m" means modified
        return Langevin(2.0 * alpha) + alpha / k0_eff;
    }
    inline __device__ float integ_phidalpha_Smith95_m(float alpha, float k0_eff) { 
        return 0.5 * Langevin_integ(2.0 * alpha) + 0.5 * alpha * alpha / k0_eff;
    }
    inline __device__ float integ_alphadphi_Smith95_m(float alpha, float k0_eff) {
        //integ actually starts from 1, but it's OK since it is for partition function calculation
        return alpha * alpha2phi_Smith95_m(alpha, k0_eff) - integ_phidalpha_Smith95_m(alpha, k0_eff);
    }
}
inline __device__ float lz_ds (float force) {//dsDNA's length per base
    return Condition::ArmLength * DNAParams::L0DS * 
            DNAModel::alpha2phi_Odijk95(force * DNAParams::LPDS / Condition::kT, DNAParams::KDS * DNAParams::LPDS / Condition::kT);
}
inline __device__ float lz_ss (float force, int j) {//ssDNA's length per base
    return 2.0 * j * DNAParams::L0SS * 
            DNAModel::alpha2phi_Smith95_m(force * DNAParams::LPSS / Condition::kT, DNAParams::KSS * DNAParams::LPSS / Condition::kT);
}
inline __device__ float le_ds (float force) {//function version of dsDNA's energy per bp:
    return Condition::ArmLength * Condition::kT * DNAParams::L0DS * 
            DNAModel::integ_alphadphi_Odijk95(force * DNAParams::LPDS / Condition::kT, DNAParams::KDS * DNAParams::LPDS / Condition::kT) / DNAParams::LPDS;
}
inline __device__ float le_ss (float force, int j) {//function version of ssDNA's energy per bp:
    return 2.0 * j * Condition::kT * DNAParams::L0SS * 
            DNAModel::integ_alphadphi_Smith95_m(force * DNAParams::LPSS / Condition::kT, DNAParams::KSS * DNAParams::LPSS / Condition::kT) / DNAParams::LPSS;
}
inline __device__ float delta_ext(float force, float j, float ext) {//func used to find force so the system total extension = ext
    return ext - force/Condition::TrapStiffness - lz_ds (force) - lz_ss (force, j);//increasing function with force
}

inline __device__ float find_force(int j, float ext) {
    float f1 = ValidRange::ValidMinForce;
    float f2 = ValidRange::VALID_MAX_FORCE;

    float y1 = delta_ext(f1, j, ext);
    float y2 = delta_ext(f2, j, ext);

    float fm = 0.0;
    float ym = 0.0;
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
    float extension_DNA = 0.0;//in nm;
    float force_average = 0.0;//in pN
    float force_SD = 0.0;//in pN
    float junzipped_average = 0.0;//#bp unzipped
    float junzipped_SD = 0.0;//#bp unzipped
};

constexpr float energy_threshold = 50.0;//don't calculate exp(-p/kT) and set probability to 0;


__global__ dp calculate_array(int extension, const std::vector<float> & seq_energy) {
        
    std::vector<float> temp_e(seq_energy.size(), 0.0);
    std::vector<float> temp_f(seq_energy.size(), 0.0);

    float min_e = 1.0e100;
    for (int j = 0; j < seq_energy.size(); ++j) {
        float f = find_force(j, extension);
        float energy = 0.0;
        if (f >= VALID_MAX_FORCE) {
            energy = std::nan("1");//make this a large number, meaning that do not use the value
        } else if (f <= VALID_MIN_FORCE) {
            energy = 0; // these DNA mechanical models are invalid at such low forces, just set energy to zero
        } else {
            energy = (0.5 * f * f / Condition::TrapStiffness + le_ds(f) + le_ss(f, j))/Condition::kT;
        }
        
        temp_f.at(j) = f;
        temp_e.at(j) = energy + seq_energy[j];

        if (min_e > temp_e.at(j)) {
            min_e = temp_e.at(j);
        }
    }

    float prob = 0;
    float Fprob = 0;
    float FFprob = 0;
    float Jprob = 0;
    float JJprob = 0;
    float p,f;
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
    point.extension_DNA = extension - point.force_average/Condition::TrapStiffness;
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
    std::chrono::duration<float> elapsed = end_time - start_time;
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
