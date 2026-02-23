# CASTEP-ADP tool
```.adp``` input file should have filename ```<seed>.adp``` (where the ```.md``` and ```.cell``` files are of the form ```<seed>.md``` and ```<seed>.cell```), this should contain paramters detailed below

to run the tool in command line: ```python3 castep-adp.py <seed>```

there is an optional ```--dryrun``` flag: ```python3 castep-adp.py <seed> --dryrun``` which will just read and check ```<seed>.adp``` for inputs and write out the header to the ```<seed>.out``` file

## input file parameters
```calculate_adp```: calculate the atomic displacement paramters, default value = True (cannot be false if ```write_jmol``` is set to True)

```calculate_ke```: calculate the kinetic energy tensor from averaged velocities from a .md file, default value = False

```equilibration_timesteps```: number of equilibration timesteps for the MD calculation, default value = 0

```r_equilibrium```: specifies the method used to calculate $r_{eq}$ for calculating $U$, default value = finite 
    
- ```finite```: calculates a finite temperature $r_{eq}$ from the MD calculation
- ```zero```: uses a zero temperature $r_{eq}$ (requires a .cell file with zero temperature positions)

```write_adp```: writes 3 orthogonal vectors for the atomic displacement parameter, default value = False

```write_ke```: writes the kinetic energy tensor $T_{ij} = \frac{m}{2}\langle v_iv_j\rangle$, default value = True if ```calculate_ke``` is True, otherwise False

```write_ke_vectors```: writes the eigenvectors from the kinetic energy tensor, default value = False

```write_jmol```: writes a jmol script (.spt) file to display the atomic displacement parameter ellipsoids in jmol, default value = False

```write_r_eq```: writes the equilibrium position, calculated by method specified in ```r_equilibrium```, default value = True

```write_uij```: writes the covariance matrix $U_{ij}=\langle \Delta x_i \Delta x_j\rangle$ for the atomic displacement parameter, default value = True if ```calculate_adp``` is True, otherwise False 