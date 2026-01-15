# athenak-tutorial-staticNS

## Install AthenaK

Follow the [wiki](https://github.com/IAS-Astrophysics/athenak/wiki/) page of AthenaK to download it. To build and install:
```
mkdir build
cd build
cmake -D Kokkos_ARCH_NATIVE=ON -D PROBLEM=dyngr_tov ../
make -j
```
The executable will be created at `build/src/athenak`. The option `-D PROBLEM=dyngr_tov` is to specify the problem generator for the static NS. Problem generator files are located under `src/pgen/`.

To compile with MPI, add the option: `-D Athena_ENABLE_MPI=ON`.

## Run simulation

First, create an individual directory under the top AthenaK directory for the specific run:
```
mkdir runs
cd runs
mkdir staticNS
cd staticNS
```
Copy the input file given by AthenaK to the current directory:
```
cp ../../inputs/dyngr/mag_tov.athinput ./
```
Then run the executable:
```
../../build/src/athena -i mag_tov.athinput time/nlim=100
```
Output files are saved under `runs/staticNS/`.

To run with MPI:
```
mpirun -np 4 ../../build/src/athena -i mag_tov.athinput
```

## Use tabulated EOS

We use [PyCompOSE](https://github.com/computationalrelativity/PyCompOSE/tree/master) to generate tabulated EoS (with the `.athtab` extension). Download the files listed in the header of the scripts in [PyCompOSE](https://github.com/computationalrelativity/PyCompOSE/tree/master), including all files from https://compose.obspm.fr, and `.h5` and `.pizza` files from Zenodo. To use the generated EoS tables, set the input argument `dyn_eos = compose`, and add `table = PATH_TO_THE_TABLE` (3-D) in the `<mhd>` and `table = PATH_TO_THE_TABLE_SLICE` (at 0 temperature) in the `<problem>` blocks.

## Plot
Create a directory for plotting scripts:
```
mkdir plots
cd plots
```
### Plot a slice
Use `plot_image.py` from https://github.com/jfields7/plot-tools. Since the generated binary files are 3-D, we need to specify a slice location to make a 2-D plot.
```
image = Image("../bin/tov.mhd_w_bcc.00001.bin", extent, "dens", slice_loc=['z', 0.0])
```
