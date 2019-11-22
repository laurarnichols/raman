title: Todo
author: Laura Nichols
date: 11/22/2019

* `../src/ramanIntensity.f90:107:` Change this to be separate variables
* `../src/ramanIntensity.f90:134:` Make this an input variable rather than hardcode
* `../src/ramanIntensity.f90:146:` Change `input.txt` to be read from input file
* `../src/ramanIntensity.f90:147:` Change output to go to command line like QE
* `../src/ramanIntensity.f90:155:` Change this to `read(12,*)`
* `../src/ramanIntensity.f90:177:` Figure out why `interval` and `count2` are allocatable
* `../src/ramanIntensity.f90:297:` Take this out of the loop
* `../src/ramanIntensity.f90:311:` Take this out of the loop
* `../src/ramanIntensity.f90:321:` Figure out where `zfactor2` comes from
* `../src/ramanIntensity.f90:323:` Since multiple them, figure out why have `exp(0.5*tmp)` as it cancels
* `../src/ramanIntensity.f90:359:` Add detailed derivation of this in a separate page
