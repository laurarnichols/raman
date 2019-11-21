title: Todo
author: Laura Nichols
date: 11/21/2019

* `../src/ramanIntensity.f90:103:` Change this to be separate variables
* `../src/ramanIntensity.f90:124:` Make this an input variable rather than hardcode
* `../src/ramanIntensity.f90:136:` Change `input.txt` to be read from input file
* `../src/ramanIntensity.f90:137:` Change output to go to command line like QE
* `../src/ramanIntensity.f90:145:` Change this to `read(12,*)`
* `../src/ramanIntensity.f90:167:` Figure out why `interval` and `count2` are allocatable
* `../src/ramanIntensity.f90:197:` Figure out why set this here as it is overwritten below
* `../src/ramanIntensity.f90:287:` Take this out of the loop
* `../src/ramanIntensity.f90:301:` Take this out of the loop
* `../src/ramanIntensity.f90:311:` Figure out where `zfactor2` comes from
* `../src/ramanIntensity.f90:313:` Since multiple them, figure out why have `exp(0.5*tmp)` as it cancels
