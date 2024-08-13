#pragma once

#include <limits>
#include <array>

//======================================parameters for constexpr LUT=========================================
//define resolution and size of the look-up tables
//j = unzipped #bp
//ext = an index to indicate the total length of the system
#ifndef J_SIZE
#define J_SIZE 200
#endif

#ifndef J_RESELUTION
#define J_RESELUTION 8000/J_SIZE
#endif

#ifndef EXT_SIZE
#define EXT_SIZE 200
#endif

#ifndef EXT_RESELUTION
#define EXT_RESELUTION 8000/EXT_SIZE
#endif

//This is the header that metaprograms the constexpr lookup table for unzipping curve calculation, it can speed up the calculation by roughly 100 times!
//
//
//
//
//
//
//

// =================================================Constants================================================
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
// ====================================valid force range for the model=======================================
//to increase speed and to make the models constexpr feasible, I simplified the models
//the modified models are only precious in a certain range
//energy/force calculation above the range may not be accurate
namespace ValidRange{
    constexpr double ValidMaxForce = 10000000.0;
    constexpr double ValidMinForce = 0.2;
}
// =====================anonymous namespace for functions used only once in this header====================
namespace {
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
    // =====================================constexpr version of math functions==================================
    namespace MyMath{
        constexpr double Inf = std::numeric_limits<double>::infinity(); 
        constexpr double VeryLargeNumber = 1.0e20;
        constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
        constexpr double Pi = 3.1415926535897932384626433832795;
        constexpr size_t precision_ln = 30; 

        //Newton-Raphson method for root calculation
        constexpr double Sqrt(double x) {
            if (x < 0) {
                return NaN;
            }
            if (x == 0) {
                return 0.0;
            }

            double prev = 1.0;
            double curr = x;
            while (curr != prev) {
                prev = curr;
                curr = 0.5 * (prev + x / prev);
            }

            return curr;
        }   
        constexpr double Cbrt(double x) {
            if (x == 0) {
                return 0.0;
            }
            double prev = 1.0;
            double curr = x;
            while (curr != prev) {
                prev = curr;
                curr = 0.3333333333333333333333333333333333333 * (2.0 * prev + x / (prev * prev));
            }
            return curr;
        }
        constexpr double Pow(double x, int n) {
            double res = 1;
            for (int i = 0; i < n; ++i){
                res *= x;
            }
            return res;
        }
        constexpr double Square(double x){
            return x * x;
        }
        constexpr double Cubic(double x) {
            return x * x * x;
        }
        constexpr double Abs(double x) {
            return x > 0.0 ? x : -x;
        }
        constexpr double Ln(double x) {
            if (x < 0.0) {
                return 1.0 / 0.0;//err;
            }
            
            if (x < 1.0) {
                return -Ln(1.0 / x);
            }
            if (x > 3.0 ) {;
                return 1.0 + Ln(x * 0.36787944117144232159552377016146);
            }
            //calculate x if 1.0 <= x <= 4.0
            double y = (x - 1) / (x + 1);
            const double ysq = y * y;
            double ypow = y; //y^+1
            double sum = ypow;
            for (size_t i = 1; i <= precision_ln; ++i) {
                ypow *= ysq;
                sum += ypow/(2.0 * i + 1);
            }
            return 2.0 * sum;
        }
        constexpr double Coth(double x) {
            if ( x > 7.9608220) {
                // I compared Coth(x) and the np.cosh(x)/np.sinh(x), and 1, and figured out this value.
                //above this value, 1.0 is more close to Coth(x).
                //this value depends on how many terms are used of course.
                return 1.0;
            }
            double product = 1.0;//x**0/0!
            double sum1 = product;
            double sum2 = 0.0f;
            int i = 1;
            for (; i <= 10; ++i){
                product /= 2*i - 1;//x**1/1!, ..
                product *= x;
                sum2 += product;

                product /= 2*i; //x**2/2!, ...
                product *= x;
                sum1 += product;
            }
            product *= x;
            product /= 2*i + 1;//x**1/1!, ..
            sum2 += product;
            
            return sum1/sum2;
        }
        // constexpr double Tanh_coarse(double x) {
        //     // the polyn formula is copied from 
        //     // https://math.stackexchange.com/questions/107292/rapid-approximation-of-tanhx
        //     // the precision is so bad
        //     return 
        //     (-.67436811832e-5+(.2468149110712040+(.583691066395175e-1+.3357335044280075e-1*x)*x)*x)/
        //     (.2464845986383725+(.609347197060491e-1+(.1086202599228572+.2874707922475963e-1*x)*x)*x);
        // }
        constexpr double Tanh(double x) {
            return 1.0 / Coth(x);
        }
        constexpr double Langevin(double x) {
            if ( x > 7.9608220) {
                // I compared Coth(x) and the np.cosh(x)/np.sinh(x), and 1, and figured out this value.
                //above this value, 1.0 is more close to Coth(x).
                //this value depends on how many terms are used of course.
                return 1.0 - 1.0 / x;
            }
            double product = x*x;//x**0/0!
            double factorial = 6;//2!
            double sum1 = 1.0/2.0 - 1.0 / 6.0;
            double sum2 = 1.0;
            for (int i = 2; i <= 10; ++i){
                sum2 += product /factorial;
                factorial *= 2.0 * i;
                sum1 += product * (1.0/factorial - 1.0 /factorial/(2.0 * i + 1.0));
                factorial *= (2.0 * i + 1.0);
                product *= x*x;
            }
            return x*sum1/sum2;
        }
        constexpr double Langevin_integ(double x) {
            // = ln(sinh(x)/x)
            double sum = 0.0;
            double factor = 1.0;
            double product = 1.0;
            for (double i = 1.0; i < 20.0; ++i) {
                sum += factor;
                factor *= (x * x /(i * 2.0) / (i * 2.0 + 1.0));
            }
            //return sum;
            return Ln(sum);
        }

        //some tests
        constexpr double Coth_test_res = Coth(0.9);//1.3960672530300118350929600819912
        static_assert(Coth_test_res > 1.3960672530);//1.3960672525087865
        static_assert(Coth_test_res < 1.3960672531);

        constexpr double Sqrt_test_res = Sqrt(4.0);
        static_assert(Sqrt_test_res > 1.9999999999);
        static_assert(Sqrt_test_res < 2.0000000001);

        constexpr double Cbrt_test_res = Cbrt(8.0);
        static_assert(Cbrt_test_res > 1.9999999999);
        static_assert(Cbrt_test_res < 2.0000000001);

        constexpr double Ln_test_res = Ln(45);
        static_assert(Ln_test_res > 3.80666);
        static_assert(Ln_test_res < 3.80667);
    }
    // ==========================================DNA mechanical models===========================================
    namespace DNAModel {
        //================================WLC/FJC model: parameter definitions==================================
        //phi = x/L (L = contour length)
        //alpha = fA/kT (A = persistence length)
        //k0_eff=k0A/kT (K0 = elastic modulus)

        //====================================Marko-Siggia 1995 WLC========================================
        constexpr double phi2alpha_MS(double phi){ 
            return phi + 0.25 / MyMath::Square(1.0 - phi) - 0.25; 
        }

        //For now I haven't implemented constexpr cos() and acos(), so I cannot make this function constexpr.
        // double alpha2phi_MS(double alpha) {//My ppt: xxxxxxxx
        //     alpha = alpha - 0.75;
        //     double p = -Square(alpha) / 3.0;
        //     double q = -2.0 * Cubic(alpha) / 27.0 + 0.25;
        //     //Cadano's formula is only correct when fA/kT-0.75 < 3.0/np.cbrt(16.0)
        //     if (alpha < 3.0/Cbrt(16.0)) {
        //         double m = Sqrt(Square(q)/4.0 + Cubic(p)/27.0);
        //         double n= -q/2.0;
        //         return 1.0 + alpha/3.0 + Cbrt(n + m) + Cbrt(n - m);
        //     } else { //use trigonometric solution
        //         return 1.0 + alpha / 3.0 + 2.0 * Sqrt(-p / 3.0) * cos(acos(1.5 * q * Sqrt(-3.0 / p) / p) / 3.0 - 2.0 * Pi * 2.0 / 3.0);
        //     }
        // }

        //==========================extensible MS, or Modified-MS model, Wang et al. 1997===============================
        //NO additional function to calculate x from f, because if we can substitute the x/L-f/K by:
        //phi = x/L-f/K and alpha = fA/kT
        //Then we can calculate:
        //constexpr double alpha2phi_MMS(){}
        //need to finish alpha2phi_MS first
        //todo: phi2alpha_MMS() will definitely need extra work.

        //==============================WLC high force, Odijk 1995 macromolecules===========================
        constexpr double alpha2phi_Odijk95(double alpha, double k0_eff) { 
            // if (alpha < 0.25) {
            //     return 0.0;// undefined at alpha == 0, just give it a small value
            //     //I can do this because this is force, energy must be calculated correctly!!
            // }
            return 1.0 - 0.5 / MyMath::Sqrt(alpha) + alpha / k0_eff; 
        }
        constexpr double integ_phidalpha_Odijk95(double alpha, double k0_eff) { 
            return alpha - MyMath::Sqrt(alpha) + 0.5 * alpha * alpha / k0_eff; 
        }
        constexpr double integ_alphadphi_Odijk95(double alpha, double k0_eff) {
            return alpha * alpha2phi_Odijk95(alpha, k0_eff) - integ_phidalpha_Odijk95(alpha, k0_eff);
        }
        // ================================MODIFIED VERSION OF FJC, Smith 1995 macromolecules================================
        //Modified version specific for ssDNA force region, and keeps accuracy
        //For ssDNA, alpha = (force * lp_ss / kT) = force /5.4, a force range of (0.1 ~ 60) is alpha < 12
        //My homemade Langevin_integ function should be accurate enough in this region.
        constexpr double alpha2phi_Smith95_m(double alpha, double k0_eff) {//"m" means modified
            return MyMath::Langevin(2.0 * alpha) + alpha / k0_eff;
        }
        constexpr double integ_phidalpha_Smith95_m(double alpha, double k0_eff) { 
            return 0.5 * MyMath::Langevin_integ(2.0 * alpha) + 0.5 * alpha * alpha / k0_eff;
        }
        constexpr double integ_alphadphi_Smith95_m(double alpha, double k0_eff) {
            //integ actually starts from 1, but it's OK since it is for partition function calculation
            return alpha * alpha2phi_Smith95_m(alpha, k0_eff) - integ_phidalpha_Smith95_m(alpha, k0_eff);
        }

        //================================FJC, Smith 1995 macromolecules================================
        //the resolution is not very good, and the calculation is too slow
        // constexpr double alpha2phi_Smith95(double alpha, double k0_eff) {
        //     return (MyMath::Coth(2.0 * alpha) - 0.5 / alpha) * (1.0 + alpha / k0_eff);
        // }
        // constexpr double integ_phidalpha_Smith95(double alpha, double k0_eff) { 
        //     const int n = 100;//integration from alpha/n to alpha!!!, not from 0 so there is a minor error
        //     const double delta = alpha/n;
        //     double sum = 0.0;
        //     double a = 0.0;
        //     for (int i = 1; i < n + 1; ++i) {//means i can = n.
        //         a = delta * i;
        //         sum += alpha2phi_Smith95(a, k0_eff);
        //     }
        //     sum -= 0.5 * alpha2phi_Smith95(delta, k0_eff);
        //     sum -= 0.5 * alpha2phi_Smith95(alpha, k0_eff);
        //     return sum * delta;
        // }
        // constexpr double integ_alphadphi_Smith95(double alpha, double k0_eff) {
        //     //integ actually starts from 1, but it's OK since it is for partition function calculation
        //     return alpha * alpha2phi_Smith95(alpha, k0_eff) - integ_phidalpha_Smith95(alpha, k0_eff);
        // }

        // ================================High force OF FJC, Smith 1995 macromolecules================================
    // The langevian function is ill formed (at least for the computer) when alpha -> 0 so avoid it=========
    // These are high force version of Smith's extensible FJC, the coth(..)  term equals to 1 when force is large
    // For ssDNA, since we also calculate dE/dJ, the low force region energy contribution thus becomes important, 
    // we cannot use this high force approximation
    // constexpr double alpha2phi_Smith95_hf(double alpha, double k0_eff) {
    //     return (1.0 - 0.5 / alpha) * (1.0 + alpha / k0_eff);
    // }
    // constexpr double integ_phidalpha_Smith95_hf(double alpha, double k0_eff) { 
    //     //remove the error of integration from zero to ValidMinAlpha and replace it with Langevin_0_to_MinAlpha
    //     return (0.5 * alpha + k0_eff - 0.5)  * alpha / k0_eff - 0.5 * MyMath::Ln(alpha);
    // }
    // constexpr double integ_alphadphi_Smith95_hf(double alpha, double k0_eff) {
    //     return alpha * alpha2phi_Smith95_hf(alpha, k0_eff) - integ_phidalpha_Smith95_hf(alpha, k0_eff);
    // }
    }

    //use a pre-calucated constexpr LUT will speed up the calculation by at least 1 order!!
    constexpr int j_size = J_SIZE;//size of the j index dimension of the LUT. basically the max length of trunk
    constexpr int ext_size = EXT_SIZE;//size of the extension dimension of LUT, the larger the number, the more precise the result.
    constexpr int ext_resolution = J_RESELUTION;//resolution of the total extension dimention: actual total extension is index * resolution
    constexpr int j_resolution = EXT_RESELUTION;//resolution of the j-index dimention: actual j is index * resolution

    using lut_type = std::array<std::array<float,ext_size>,j_size>;

    //=================================utility functions to create constexpr lut=================================
    constexpr double lz_ds (double force) {//dsDNA's length per base
        return Condition::ArmLength * DNAParams::L0DS * 
                DNAModel::alpha2phi_Odijk95(force * DNAParams::LPDS / Condition::kT, DNAParams::KDS * DNAParams::LPDS / Condition::kT);
    }
    constexpr double lz_ss (double force, int j) {//ssDNA's length per base
        return 2.0 * j * DNAParams::L0SS * 
                DNAModel::alpha2phi_Smith95_m(force * DNAParams::LPSS / Condition::kT, DNAParams::KSS * DNAParams::LPSS / Condition::kT);
    }
    constexpr double le_ds (double force) {//function version of dsDNA's energy per bp:
        return Condition::ArmLength * Condition::kT * DNAParams::L0DS * 
                DNAModel::integ_alphadphi_Odijk95(force * DNAParams::LPDS / Condition::kT, DNAParams::KDS * DNAParams::LPDS / Condition::kT) / DNAParams::LPDS;
    }
    constexpr double le_ss (double force, int j) {//function version of ssDNA's energy per bp:
        return 2.0 * j * Condition::kT * DNAParams::L0SS * 
                DNAModel::integ_alphadphi_Smith95_m(force * DNAParams::LPSS / Condition::kT, DNAParams::KSS * DNAParams::LPSS / Condition::kT) / DNAParams::LPSS;
    }
    constexpr double delta_ext(double force, double j, double ext) {//func used to find force so the system total extension = ext
        return ext - force/Condition::PillarStiffness - lz_ds (force) - lz_ss (force, j);//increasing function with force
    }

    constexpr double tor_binary_search = 1.0e-3; // 0.001 nm
    constexpr double find_force(int j, double ext) {//j == length of unzipped trunk, ext total extension
        //simple binary search to get force so calc_z_diff(force) == 0
        //force function must be monotonic
        double f1 = ValidRange::ValidMinForce;
        double f2 = ValidRange::ValidMaxForce;

        double y1 = delta_ext(f1, j, ext);
        double y2 = delta_ext(f2, j, ext);

        if (y1 * y2 >= 0) {
            if (y1 < 0){
                return ValidRange::ValidMaxForce + 1.0;//force is too large
            } else {
                return ValidRange::ValidMinForce - 1.0;//force is too small
            }
        }

        int cnt = 0;//in case the root is not found
        while (++cnt < 10000) {
            
            double fm = (f1 + f2) * 0.5;
            double ym = delta_ext(fm, j, ext);
            
            if (MyMath::Abs(ym) <= tor_binary_search) {
                return fm;
            }

            //(&& has higher precedence)
            if (y1 < y2 && ym > 0 || y1 > y2 && ym < 0) {
                f2 = fm;
                y2 = ym;
            } else if (y1 < y2 && ym < 0 || y1 > y2 && ym > 0) {
                f1 = fm;
                y1 = ym;
            } else {
                return MyMath::VeryLargeNumber;//means some weird error
            }
        }
        return MyMath::VeryLargeNumber;//meaning that the root is not found
    }

    // ===================calculate force look-up table=========================
    constexpr lut_type Lut_force = []{

        lut_type arr;

        for (int j = 0; j < j_size; ++j) {
            for (int k = 0; k < ext_size; ++k) {
                arr[j][k] = find_force(j * j_resolution, k * ext_resolution);
            }
        }
        return arr;
    }();
    // ====================calculate energy look-up table=======================
    constexpr lut_type Lut_energy = []{

        lut_type arr;

        double f = 0.0;
        for (int j = 0; j < j_size; ++j) {
            for (int k = 0; k < ext_size; ++k) {
                f = Lut_force[j][k];
                if (f >= ValidRange::ValidMaxForce || f <= ValidRange::ValidMinForce) {
                    arr[j][k] = MyMath::VeryLargeNumber;//make this a large number, meaning that do not use the value
                } else {
                    arr[j][k] = (0.5 * f * f / Condition::PillarStiffness + le_ds(f) + le_ss(f, j * j_resolution))/Condition::kT;
                }
            }
        }

        return arr;
    }();
}
// =====================construct constexpr lut class objs==================
constexpr class lookup_class {
public:
    constexpr lookup_class(lut_type lut_in) : lut(lut_in) {};
    constexpr double operator() (double j0, double extension) const {
        return lookup (j0, extension);
    };
    constexpr lut_type  get_lut() const {
        return lut;
    }

private:
    const lut_type lut;
    constexpr double lookup(double j0, double extension) const {

        if (j0 < 0 || j0 >= j_size * j_resolution || extension < 0 || extension >= ext_size * ext_resolution) {
            //error
            return -1.0;
        }

        double j = j0 / static_cast<double>(j_resolution);
        double k = extension / static_cast<double>(ext_resolution);
        
        //because j and k are positive.
        //static_cast<int> has the same result as std::floor() from cmath lib
        //Legend says that cast is 3 times faster
        double j1 = static_cast<int>(j);
        double k1 = static_cast<int>(k);

        double j2 = j1 + 1;
        double k2 = k1 + 1;

        //https://en.wikipedia.org/wiki/Bilinear_interpolation
        return lut[j1][k1] * (j2 - j) * (k2 - k) + 
                lut[j1][k2] * (j2 - j) * (k - k1) + 
                lut[j2][k1] * (j - j1) * (k2 - k) + 
                lut[j2][k2] * (j - j1) * (k - k1); 
        
    }
} Force {Lut_force}, Energy {Lut_energy};
// constexpr lookup_class Force {Lut_force};
// constexpr lookup_class Energy {Lut_energy};