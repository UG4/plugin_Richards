# Richards plugin #

**Richards** is a **UG4-Plugin** providing tools for the implementation of unsaturated flow.
Copyright 2019-2023 UG4 Developers, University Frankfurt

### Available Models ###

This provides saturation and (relative) permeabilities for the following models:

* van Genuchten-Mualem

$$
	Seff:= (1/1+(alpha * psi )^n)^m,
 	S    = thetaR + (thetaS-thetaR)*Seff
 
$$ 
* Exponential
$$ 
	Seff = exp(-p.alpha * (pc/p0)
	S    = thetaR + (thetaS-thetaR)*Seff
 $$
 
* ...

Export is via UserDataNumber. JSON constructors and factories are provided.

### Installation ###
Please install/clone this repository through UG4's package manager
[ughub](https://github.com/UG4/ughub).

### Dependencies ###

Depends on:
* [ugcore](https://github.com/UG4/ugcore).
* [AutodiffForUG4](https://github.com/UG4/external_AutodiffForUG4).
* [JSONForUG4](https://github.com/UG4/external_JSON).

Requires a C++-17 compiler.

cmake -DUSE_JSON=ON -DUSE_AUTODIFF=ON
