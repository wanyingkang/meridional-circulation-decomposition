# meridional-circulation-decomposition
This code computes the meridional streamfunction and decomposes it into components driven by eddy momentum/heat transport, adiabatic heating and surface friction using Kuo-Eliasson Equation (see Kuo 1956); and it can be used to reproduce the results in Kang et al. 2019 (necessary data included).

The data input should be the time-averaged 3D output of a GCM with regular grid. The following variables are required: 

V: meridional air motion averaged over a period of time

U: zonal air motion averaged over a period of time

T: temperature averaged over a period of time

OMEGA: vertical motion in pressure coordinate averaged over a period of time

VU: V*U for each time step averaged over a period of time

VT: V*T for each time step averaged over a period of time

OMEGAT: OMEGA*T for each time step averaged over a period of time

OMEGAU: OMEGA*U for each time step averaged over a period of time

PS: surface pressure averaged over a period of time

TREF: the reference temperature profile. This only applies to the Held-Suarez setup where temperature is simply relaxed toward a given profile Tref, and is used to compuate the parameterized heating rate. For more realistic GCM, this heat rate should be replaced by individual adiabatic heat source, including latent heating, radiative cooling etc. 



Kuo, H.L., 1956. Forced and free meridional circulations in the atmosphere. Journal of Atmospheric Sciences 13, 561–568.

Kang, W., M. Cai and E. Tziperman, 2019, Tropical and Extratropical General Circulation with a Meridional Reversed Temperature Gradient as Expected in a High Obliquity Planet, Icarus
