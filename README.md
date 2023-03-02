# MLD&MTD  
<a href="https://github.com/romeroqe/mld-mtd"><img src="https://shields.io/github/v/release/romeroqe/mld-mtd" alt="Release"></a>
<a href="http://creativecommons.org/licenses/by/4.0/"><img src="https://shields.io/github/license/romeroqe/mld-mtd" alt="License"></a>
<a href="https://zenodo.org/badge/latestdoi/524132263"><img src="https://zenodo.org/badge/524132263.svg" alt="DOI"></a>

This is a methodology to locate the minimum and maximum depth of the thermocline, its thickness, and its strength by fitting the sigmoid function to the temperature profiles in the global ocean. This methodology can be applied to both global and local studies in those areas of the ocean where the water column can be divided into three layers according to its thermal structure.

## Installation
To use this methodology, simply download the `mldmtd.py` file.

## Demo
To locate the minimum and maximum depth of the thermocline of an Argo float profile::

```python
from mldmtd import getProfileDataFromArgoNc, thermocline

df = getProfileDataFromArgoNc('.../aoml/4900432/profiles/D4900432_106.nc')
pres_mtd, temp_mtd, pres_mld, temp_mld, r2, N2T, pres_pred, temp_pred = thermocline(df)
```

## License
  
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
