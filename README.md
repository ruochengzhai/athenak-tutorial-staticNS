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

## Run the simulation

First, create an individual directory under the top AthenaK directory for the specific run:
```
mkdir runs
cd runs
mkdir staticNS
cd staticNS
```
Then run the executable:
```
../../build/src/athena -i ../../inputs/dyngr/mag_tov.athinput time/nlim=100
```
Output files are saved under `runs/staticNS/`.

