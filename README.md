# CASTEP-ADP tool


## input file parameters
```equilibration_timesteps```: number of equilibration timesteps for the MD calculation, default value = 0

```r_equilibrium```: specifies the method used to calculate $r_{eq}$ for calculating $U$, default value = finite 
    
- ```finite```: calculates a finite temperature $r_{eq}$ from the MD calculation
- ```zero```: uses a zero temperature $r_{eq}$ (requires a .cell file with zero temperature positions)


```write_jmol```: writes a jmol script (.spt) file to display the atomic displacement parameter ellipsoids in jmol, default value = False