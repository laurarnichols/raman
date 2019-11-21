title: Todo
author: Laura Nichols
date: 11/21/2019

* `../src/ramanIntensity.f90:106:` Change this to be separate variables
* `../src/ramanIntensity.f90:130:` Make this an input variable rather than hardcode
* `../src/ramanIntensity.f90:142:` Change `input.txt` to be read from input file
* `../src/ramanIntensity.f90:143:` Change output to go to command line like QE
* `../src/ramanIntensity.f90:151:` Change this to `read(12,*)`
* `../src/ramanIntensity.f90:173:` Figure out why `interval` and `count2` are allocatable
* `../src/ramanIntensity.f90:292:` Take this out of the loop
* `../src/ramanIntensity.f90:306:` Take this out of the loop
* `../src/ramanIntensity.f90:316:` Figure out where `zfactor2` comes from
* `../src/ramanIntensity.f90:318:` Since multiple them, figure out why have `exp(0.5*tmp)` as it cancels
* `../src/ramanIntensity.f90:354:` Add detailed derivation of this in a separate page
