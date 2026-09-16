# Richards plugin #

**Richards** is a **UG4-Plugin** providing tools for the implementation of unsaturated flow.

## License
Copyright 2019-2023 UG4 Developers, University Frankfurt.

The plugin is available under [EUPL v1.2](https://joinup.ec.europa.eu/sites/default/files/custom-page/attachment/2020-03/EUPL-1.2%20EN.txt).

Details on UG4 can be found [here](https://github.com/UG4/ugcore).

### Available Models ###

This provides saturation and (relative) permeabilities for the following models:

* van Genuchten-Mualem
* Exponential
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
```
cmake -DUSE_JSON=ON -DUSE_AUTODIFF=ON
```
