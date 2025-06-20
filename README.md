# Efficient-and-faithful-reconstruction-of-dynamical-attractors-using-Differentiators
Code used to generate the figures in Main Text and Supplementary Material of the paper "Efficient and faithful reconstruction of dynamical attractors using Differentiators" by U. Sutulovic, D. Proverbio, R. Katz, and G. Giordano.
Supplementary material for the paper can be found at the Zenodo link https://zenodo.org/records/15706888 .

# License

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

The full text of the GNU General Public License can be found in the file "LICENSE.txt".

# Usage

- `Example_Lorenz.m` : Computes an estimate of the dynamical attractor of the Lorenz system using the Differentiator. Set 'include_grassberger' to 1 to also include attractor reconstruction by means of the Schreiber-Grassberger method. Select the noise type by setting the variable 'noise_type' and select if it corrupts the measurement as additive or multiplicative noise with the variable 'additive_noise'. Before the section 'Lorenz system simulation', it is also possible to set the Differentiator parameters and the simulation options.

- `Example_Hindmarsh_Rose.m` : Computes an estimate of the dynamical attractor of the Hindmarsh-Rose system using the Differentiator. Select the dynamical regime with the variable 'regime'. Select the noise type by setting the variable 'noise_type' and select if it corrupts the measurement as additive or multiplicative noise with the variable 'additive_noise'. Before the section 'Hindmarsh-Rose system simulation', it is also possible to set the Differentiator parameters and the simulation options.

- `Example_Epileptor.m` : Computes an estimate of the dynamical attractor of the Epileptor system using the Differentiator. Select the noise type by setting the variable 'noise_type' and select if it corrupts the measurement as additive or multiplicative noise with the variable 'additive_noise'. Before the section 'Epileptor model simulation', it is also possible to set the Differentiator parameters and the simulation options.

- `systems` : This folder contains the functions to simulate the systems considered in the examples.

- `util` : This folder contains utility functions to retrieve the paper results, including the Differentiator and the Schreiber-Grassberger method.
