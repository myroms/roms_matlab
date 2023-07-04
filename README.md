# **ROMS Matlab Processing Scripts**

![roms-jedi logo](https://www.myroms.org/trac/roms_src_600px.png)

The Regional Ocean Modeling System (**ROMS**) framework is intended for users
interested in ocean modeling. Please check https://github.com/myroms for
instructions on registering into the **ROMS** community and downloading its source
code and test repositories.

This repository contains useful Matlab scripts that can be used for configuring
**ROMS** applications and pre-and post-processing input and output data. Use the
following command to download the **ROMS** Matlab repository:
```
git clone https://github.com/myroms/roms_matlab.git
```
| Tools   |  Description |
| ------------- | -------------  |
| **4dvar**         | **ROMS** 4D-Var data assimilation observations processing |
| **bathymetry** | **ROMS** bathymetry extraction and processing |
| **bin** | **CSH** and **BASH** scripts |
| **boundary** | **ROMS** lateral boundaries conditions processing  |
| **coastlines** | **ROMS** application coastline processing |
| **colormaps** | Color palettes for plotting |
| **coupling** | Coupling melding weights to combine **DATA** and **ESM** components |
| **forcing** | **ROMS** atmospheric fields processing |
| **grid** | **ROMS** application grid processing |
| **initial** | **ROMS** initial conditions processing |
| **ioda** | **ROMS-JEDI** data assimilation observations in **IODA** format |
| **landmask** | **ROMS** application grid land/sea mask processing |
| **mex** | Deprecated NetCDF interface to Matlab |
| **netcdf** | Usefull NetCDF processing scripts |
| **m_map** | Mapping package for Matlab developed at **UBC** |
| **seagrid** | Deprecated **ROMS** grid generation tool |
| **seawater** | CSIRO Seawater library functions |
| **tidal_ellipse** | Tidal parameter conversion |
| **t_tide** | **ROMS** tidal forcing processing |
| **utility** | Miscellaneous **ROMS** processing and plotting scripts |

Notice that other **Matlab** and **Python** tools are available in the ocean
community and can be used to process **ROMS** data. It is up to the user to
select the ones more appropriate for their applications. However, some **grid**
and **utility** sub-directories scripts are the official version of **ROMS** for
computing numerical kernel configurations and fields.

The **doxygen** version of **ROMS** is available at:
```
https://www.myroms.org/doxygen
```
The **WikiROMS** documentation and tutorials portal is available at:
```
https://www.myroms.org/wiki
```
