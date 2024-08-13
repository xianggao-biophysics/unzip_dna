# DNA Unzipping Curve Calculator (CPU version)  

*At the current stage after my optimization, this program is 1,000-2,000 times faster than my proof-of-concept python code.*  

## How to build the project  

The build command should work on Windows and Linux:  

```bash
>g++ -std=c++20 -fconstexpr-ops-limit=100000000000 main.cpp utils.cpp -DJ_SIZE=200 -DEXT_SIZE=200 -o main.exe
```
The lookup table size is controlled by macro **J_SIZE**. The larger the number, the more accurate the result. Increasing this number can significantly increase the compiling time.  
  
For Windows OS, run **build_project.bat** to build with *additional options and testing*. For Linux user, ChatGPT should be able to translate the .bat file to a shell script. If the .bat file is ran successfully, you should see something like this:  
![image](tutorial.png)  
The compile time is quite long if the look-up table size is big!  

## Do some tests using the example data

A test sequence is provided (NEB_H5alpha_Accessory_colonization_factor_AcfD). You can run the executable on this example like this:  

```bash
>main.exe NEB_H5alpha_Accessory_colonization_factor_AcfD.txt out.csv
```

The first argument is input file name. The second argument determines the output file name and is optional.  
Run "plot.py" to plot the test sequence's result.  

## Goal of this program  

My goal is to make the unzipping curve calculation fast (so I can calculate the unzipping curves of thousands of genes in an acceptable time). However, there is no better method other than brute-force partition function calculation for now. The only thing I can do is to make each loop faster. I decided to calculate something ahead of time and save it in the program as constexpr variables.

It took some thinking to move majority of the calculation from run-time to compile time. After several attempts, I **"constexpred"** most of the calculation overhead. Two look-up tables (LUTs) are created to hold these data. These LUTs are saved in **constexpr std::arrays** so I have to use c++20 (or above). The drawback is that the compile time is very long, thousands of times longer than a straightforward c++ program.  
  
On Aug/15/2023, I implemented multithreading. The execution speed increased by another factor of 10-20.  

## DNA unzipping theory  

![image](theory_vs_experiment_v1.png)

**Figure above shows DNA unzipping experiment on a 4.4 kb DNA**. Single-molecule measurement (blue) and the theoretically prediced curve (black) agree well.  
  
Further reading on DNA unzipping experiment and theory:  

[1] Essevaz-Roulet, Baptiste, Ulrich Bockelmann, and Francois Heslot. (1997) PNAS  
[2] Bockelmann, Ulrich, et al. (2002) Biophysical journal  
[3] Huguet, Josep M., et al. (2010) PNAS  


# DNA Unzipping Curve Calculator (GPU version)  

This is my first time using CUDA. It runs pretty fast and can finish the whole *E.coli*  genome in about a minute on my GeForce RTX 3070 Laptop.  

## Instructions to Use (with Example Data)  

1. The example genome is NEB H5alpha, the genbank file is in "examples" folder.  
2. Run **parse_h5alpha.py** to parse the *200304_NEB_H5alpha_genbank.txt* and *200304_NEB_H5alpha.txt* to individual files. These files should be saved in folder "parsed"  
3. Build the project on Windows. This .bat file will automatically run a test using the sequences in "parse" folder  

```bash
> build.bat
```
To batch-process sequence files in your own folder, use this command  

```bash
> main path/to/your/own/folder
```

Result (I only showed the first 23 genes):  

![image](examples/result.svg)

Prediction vs Experiment:

![image](reference/theory_vs_experiment.png)

**DNA unzipping experiment vs theory**. The prediction (**${\color{black}-}$**) aligns well with the single-molecule unzipping curve (**${\color{red}-}$**).  

## Why CUDA  

Single-molecule approaches such as optical tweezers, can unzip a single dsDNA at single-molecular level.

![image](reference/sm_DNA_unzipping_exp_schematics.png)  

The theoretical prediction of an unzipping trace (*unzipping force vs total extension*) can be obtained from partition functin of the system. The partition function $Z$ at a total extension of the system $z$ is

$$Z(z) = \sum_{j=0}^{j_{max}}e^{-G(j,z)/kT}$$

where $G(j,z)$ is the total energy of the system. Force $F(z)$ can be obtained from partition function as follows:  

$$F(z) = -kT\frac{\partial }{\partial z}\mathrm{ln}Z$$

To calculate the force-extension curve, we need to obtain the energy $G(j,z)$ at every j and z. However, $G$ has a **non-analytical complex form**, while the scales of $j$ and $z$ are large (usually in a range of 1,000-10,000). Moreover, we need to calculate the unzipping traces for every gene in the genome. Even *E. Coli* has thousands of genes. Therefore, Calculation of $G$ for all $j$ and $z$ for all DNA sequences is better to be on a GPU.  

*(Besides using GPU, I also further sped up the calculation using a trick I developed in my previous project [unzipDNA_CPU](https://github.com/Taomihog/unzipDNA_CPU).)*

## Further Reading  

[1] Essevaz-Roulet, Baptiste, Ulrich Bockelmann, and Francois Heslot. (1997) PNAS  
[2] Bockelmann, Ulrich, et al. (2002) Biophysical journal  
[3] Huguet, Josep M., et al. (2010) PNAS  
