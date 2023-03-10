# 111-2-MS_Homework #1

> 學期：111-2 

> 課程名稱：分子模擬（Molecular Simulation）

> 授課教師：Shiang-Tai Lin (林祥泰)\
> Email: stlin@ntu.edu.tw


## Homework #1: Compute system energy from force field

- Due: March 9 (Thursday)

### 1. Consider a water molecule with the coordinates of each atom (in unit of Angstroms) shown in the figure bellow

H (0, 0, 1)
O (0, 0, 0)
H (0, 0.9, -0.35)

Suppose the atomic charge (in unit of electrons) on oxygen atom (red) is -0.82 and on hydrogen atom is 0.41, and the forcefield parameters to be used are as follows: 

- the force constant for OH bond is Kb=500 kcal/mol Å2 and the equilibrium bond length is 1 Å
- the force constant for HOH angle bend is K0=120 kcal/mol and the equilibrium angle is 109.470
- the Lennard-Jones-12-6 parameters for oxygen are (R0=3.5532 Å, D0=0.1848 kcal/mol) and for hydrogen are (R0=0.9 Å, D0=0.01 kcal/mol)

(a) Use the BOND class learned in class and provide the positions of the three atoms of water to determine the length and energy of the two OH bonds.

(b) Use the ANGLE class to determine the included angle and the energy of the angle between H-O-H.

(c) Report the energy components (bond, cosine harmonic angle, torsion, inversion, van der Waals, and coulomb) and the total energy of this water molecule.

(d) Base on this forcefield, in what situation will the total energy of water become zero? Comment on the meaning of the total energy you find in (c).

### 2. A cluster of two water molecules is shown in the figure below. Suppose the atomic charge, force field parameters for bond and angle terms are the same as those in problem 1. Modify the code so that it determines the total energy of the water cluster.

H (0.8, 0, 3.95)
H (-0.8, 0, 3.95) O (0, 0, 3.36)

H (0, 0, 1)
O (0, 0, 0)
H (0, 0.9, -0.4)


(a) Report the internal energy of each of the water molecules

(b) Report the total energy of the water cluster

(c) Report the binding energy of the two water molecules (i.e., the difference of energies from (a) and (b)).
