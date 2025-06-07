# Commands
- `gfortran -o q2df_generalization precision.f90 q2df_generalization.f90`
- `./q2df_generalization`
- For debugging: `gfortran -g -fcheck=all -Wall -fbacktrace -o q2df_generalization precision.f90 q2df_generalization.f90`

## Other:
- Note you cannot do huge(1.0_WP) - 1.0_WP, have to use multiplication
- 1473 failures with robust (up from 22)