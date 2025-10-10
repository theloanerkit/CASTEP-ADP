# CASTEP-ADP tool
```.adp``` input file should have filename ```<seed>.adp```

## input file parameters
```calculate_adp```: calculate the atomic displacement paramters, default value = True (cannot be false if ```write_jmol``` is set to True)

```calculate_ke```: calculate the kinetic energy tensor from averaged velocities from a .md file, default value = False (written out in electron volts)

```equilibration_timesteps```: number of equilibration timesteps for the MD calculation, default value = 0

```r_equilibrium```: specifies the method used to calculate $r_{eq}$ for calculating $U$, default value = finite 
    
- ```finite```: calculates a finite temperature $r_{eq}$ from the MD calculation
- ```zero```: uses a zero temperature $r_{eq}$ (requires a .cell file with zero temperature positions)


```write_jmol```: writes a jmol script (.spt) file to display the atomic displacement parameter ellipsoids in jmol, default value = False