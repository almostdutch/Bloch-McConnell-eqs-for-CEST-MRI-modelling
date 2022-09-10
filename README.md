# Bloch McConnell equations simulator

![](https://github.com/almostdutch/Bloch-McConnell-eqs-for-CEST-MRI-modelling/blob/master/cest_movie.gif)

Program to model CEST-MRI effects with an arbitrary number of pools using a fully-scalable system of differential equations


For the explanation of different CEST-MRI sequences, see Khlebnikov, V. et al. Comparison of pulsed three-dimensional CEST acquisition schemes at 7 tesla: 
steady state versus pseudosteady state. Magn. Reson. Med. 77, 2280–2287 (2017).

For the full metabolic analysis of CEST-MRI spectra, see Khlebnikov V. et al. Analysis of chemical exchange saturation transfer contributions from brain metabolites 
to the Z-spectra at various field strengths and pH. Sci. Reports. 9, 1098 (2019). 


cest_master_ss.m is a master script containing sequence parameters for a steady-state CEST-MRI sequences
cest_master_ps.m is a master script containing sequence parameters for a pseudosteady-state CEST-MRI sequences

Folder **pulses** contains RF profiles of preparation (MT) and excitation pulses\
Folder **routines** contains the core routines of the simulator

**Possibilities:** the simulator is implemented as a generalized fully scalable system of differential equations for an arbitrary number of CEST-MRI pools, 
which means that absolutely any (!) sequence parameters can be editted and tuned at the user descretion. In addition to simulating the CEST-MRI spectrum, 
the simulator can analyze the spectrum in terms of the underlying basis functions that constitute the spectrum. Furthermore, it's possible to monitor 
the full evolution of Mz component over time, which can be important, for example, to check if a steady-state is reached.

Example of simulating a steady-state CEST-MRI sequence:
![](https://github.com/almostdutch/Bloch-McConnell-eqs-for-CEST-MRI-modelling/blob/master/ss-CEST-MRI-sequences.jpg)

Example of simulating a pseudosteady-state CEST-MRI sequence:
![](https://github.com/almostdutch/Bloch-McConnell-eqs-for-CEST-MRI-modelling/blob/master/ps-CEST-MRI-sequences.jpg)





