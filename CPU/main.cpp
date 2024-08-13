#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <cmath>
#include <chrono>
#include <string>
#include <cmath>

#include "ThreadPool.h"//https://github.com/progschj/ThreadPool
#include "lookup_table_generator.h"//pure header to generate look up data to speed up calculation

//==========================================what is DNA unzipping experiment==============================================
//the DNA structure used in an unzipping experiment is something like this:
//
//            Force
//             ↑
//             ||
//             ||
//               =================
//             ||
//             ||
//             ↓
//            Force
//
//This is called a "Y structure". A Y structure has 2 dsDNA arms and a dsDNA trunk.
//when the two arms are stretched as shown, the trunk will be unzipped (ie, paired nucleotides will be separated).
//the force during this "unzip" precess is measured against the end-to-end distance of the arms. 
//The sequence of the trunk determines the "force vs extension" profile.
//AND can be calculated theoretically.

//for more information, see:

//[1] Essevaz-Roulet, Baptiste, Ulrich Bockelmann, and Francois Heslot. (1997) PNAS.
//[2] Bockelmann, Ulrich, et al. (2002) Biophysical journal.
//[3] Huguet, Josep M., et al. (2010) PNAS.


// Define an alias for the high-resolution clock
using Clock = std::chrono::high_resolution_clock;

struct dp {// a data point
    int extension_total = 0;//in nm;
    double extension_DNA = 0.0;//in nm;
    double force_average = 0.0;//in pN
    double force_SD = 0.0;//in pN
    double junzipped_average = 0.0;//#bp unzipped
    double junzipped_SD = 0.0;//#bp unzipped
};

//Calculate DNA sequence's energy
namespace DNAsequence {
    const double EffSaltConc = MyMath::Ln(Condition::SaltConc * 0.001) / 298.0;//298.0 K is where the energy was measured in Huguet paper, it is not equal to Condition::temperature
    const double factor_pNnm = Const::pNnm * Const::Joul / Const::Avogadro / Condition::kT;//convert to 'pN-nm' unit;
}
constexpr double energy_threshold = 50.0;//don't calculate exp(-e/kT) and set probability to 0;

//=============================================util functions==========================================
//single read. todo: read the fasta directly, multiple line read.
std::string readtxt_firstline(const std::string & path) {
    std::ifstream file(path);

    if (!file.is_open()) {
        std::cerr << "Error opening the file." << std::endl;
        return "0";
    }

    std::string line;
    if (getline(file, line)) {
        std::cout << "Trunk sequence '" + path + "' read successfully." << std::endl;
    } else {
        std::cerr << "Failed to read a line from the file." << std::endl;
    }

    file.close();

    return line;
}
//generate a file path for output result from the input file path.
std::string creat_path_out(const std::string & path_in)
{
    size_t lastSlashIndex = path_in.rfind('/');

    std::string parentPath;
    if (lastSlashIndex != std::string::npos) {
        parentPath = path_in.substr(0, lastSlashIndex) + "/";
    }

    std::string fullname = path_in.substr(lastSlashIndex + 1, std::string::npos);

    return parentPath + "out_" + fullname.substr(0, fullname.rfind('.')) + ".csv";;
}
namespace DNAsequence {
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
dp calculate_array(int extension, const std::vector<double> & seq_energy) {
        
    std::vector<double> temp_e(seq_energy.size(), 0.0);
    std::vector<double> temp_f(seq_energy.size(), 0.0);

    double min_e = 1.0e100;
    for (int j = 0; j < seq_energy.size(); ++j) {
        temp_f.at(j) = Force(j,extension);
        temp_e.at(j) = Energy(j,extension);
        // temp_f.at(j) = lookup(j,extension,Lut_force);
        // temp_e.at(j) = lookup(j,extension,Lut_energy);
        
        temp_e.at(j) += seq_energy[j];

        if (min_e > temp_e.at(j)) {
            min_e = temp_e.at(j);
        }
    }

    double prob = 0;
    double Fprob = 0;
    double FFprob = 0;
    double Jprob = 0;
    double JJprob = 0;
    double temp1,temp2,temp3;
    for (int j = 0; j < seq_energy.size(); ++j) {

        temp1 = temp_e.at(j) - min_e;
        temp1 = temp1 > energy_threshold ? 0.0 : std::exp(-temp1);
        temp2 = temp_f.at(j);

        prob += temp1;

        temp3 = temp2 * temp1;
        Fprob += temp3;
        FFprob += temp3 * temp2;

        temp3 = j * temp1;
        Jprob += temp3;
        JJprob += j * temp3;

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

//==================================================main===========================================
int main(int argc, char * argv[]) {

    auto start = Clock::now();

    if (argc <2) {
        std::cerr << "Please provide a filename (argv[1])";
        return 1;
    }
    
    const std::string sequence = readtxt_firstline(argv[1]);
    const std::vector<double> seq_energy = DNAsequence::calculate_sequence_energy(sequence);
    //test, the "41740.760955375" is from my python code
    //std::cout << "Sequence energy differs from python program by "<< *(seq_energy.cend() - 1) - 10140.0933068047 << std::endl;// no difference, good
    
    int numThreads = std::thread::hardware_concurrency();
    std::cout << "Number of threads: " << numThreads << std::endl;
    ThreadPool pool(numThreads);

    std::vector< std::future< dp > > results;

    for(int extension = 1; extension < static_cast<int>(1.2 * seq_energy.size()); ++extension) {
        results.emplace_back(
            pool.enqueue([extension, &seq_energy]{
                return calculate_array(extension, seq_energy);
            })
        );
    }

    std::vector<dp> result_array;
    for(auto && result: results)
        result_array.emplace_back(result.get());

    //prepare file for output
    const std::string path_out = argc >= 3 ? argv[2] : creat_path_out(argv[1]);

    std::ofstream fout(path_out);
    if (!fout) {
        std::cerr << "Error opening file. Cannot save result" << std::endl;
        return 1;
    }

    fout << "total extension (nm),DNA extension (nm),average force (pN),sd force (pN),average bp unzipped,sd bp unzipped" << std::endl;

    char delimit = ',';
    std::for_each(result_array.cbegin(), result_array.cend(), [&fout, delimit](const dp & point) {
        fout << point.extension_total << delimit;
        fout << point.extension_DNA << delimit;
        fout << point.force_average << delimit;
        fout << point.force_SD << delimit;
        fout << point.junzipped_average << delimit;
        fout << point.junzipped_SD << std::endl;;
    });

    // Close the file
    fout.close();
    std::cout << "Result has been written to '" + path_out << "'." <<std::endl;
    
    auto end = Clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "Execution time: " << static_cast<double>(elapsedTime) / 1'000'000.0 << " s." << std::endl;
    
    return 0;
}