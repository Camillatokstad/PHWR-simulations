from cProfile import label
import openmc 
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd 
import openmc.deplete

sp = openmc.StatePoint('statepoint.100.h5')
tally = sp.tallies[1]

flux = tally.get_slice(scores=['flux'])
prompt = tally.get_slice(scores=['prompt-nu-fission'])
flux.std_dev.shape = (100, 100)
flux.mean.shape = (100, 100)
prompt.std_dev.shape = (100, 100)
prompt.mean.shape = (100, 100)
# Plotting the results
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
#im1 = ax1.imshow(flux.mean, extent=[-166.6, 166.6, -166.6, 166.6], origin='lower')
im1 = ax1.imshow(flux.mean, extent = [-0.7,0.7,-0.7,0.7], origin= "lower")
ax1.set_title("Neutron Flux Distribution")
#plt.colorbar(im1, ax=ax1)

#im2 = ax2.imshow(prompt.mean, extent=[-166.6, 166.6, -166.6, 166.6], origin='lower')
im2 = ax2.imshow(prompt.mean, extent = [-0.7,0.7,-0.7,0.7], origin= "lower")
ax2.set_title("Fission Site Heatmap")
#plt.colorbar(im2, ax=ax2)



plt.show()

print(sp.tallies)
erange = np.logspace(np.log10(10e-5),np.log10(20e6),501)
tally2 = sp.tallies[2]
flx = tally2.mean.ravel()
plt.loglog(erange[:-1], flx) 
plt.grid()
plt.xlabel("Energy eV") 
plt.ylabel("Flux [n/cm-src]") 
plt.title("Neutron energy spectrum") 
plt.show() #plt.savefig() for heplab

"""
# Compute resonance escape probability using tally arithmetic
therm_abs_rate = sp.get_tally(name='therm. abs. rate')

fuel_therm_abs_rate = sp.get_tally(name='fuel therm. abs. rate')
therm_util = fuel_therm_abs_rate / therm_abs_rate

values = therm_util.get_pandas_dataframe()

print(therm_util)
print(values)
"""
#k-value over time
results = openmc.deplete.ResultsList.from_hdf5("./depletion_results.h5")
time, k = results.get_eigenvalue()
time /= (24 * 60 * 60)  # convert back to days from seconds

plt.errorbar(time, k[:, 0], yerr=k[:, 1])
plt.xlabel("Time [d]")
plt.ylabel("$k_{eff}\pm \sigma$")
plt.show()


#concentration of elements
_time, u5 = results.get_atoms("1", "U235")
_time, pu239 = results.get_atoms("1", "Pu239")
_time, u8 = results.get_atoms("1", "U238")
_time, pu241 = results.get_atoms("1", "Pu241")

plt.plot(time, u5, label="U235")
plt.plot(time,pu239, label="Pu239")
plt.xlabel("Time [d]")
plt.ylabel("Number of atoms")
plt.legend()
plt.show()


_time, u5_fission = results.get_reaction_rate("1", "U235", "fission")
plt.plot(time, u5_fission, label="U235 fission events")
plt.xlabel("Time [d]")
plt.ylabel("Fission reactions / s")
plt.legend()
plt.show()

plt.plot(time,pu239/u5, label="Pu239/U235")
plt.xlabel("Time [d]")
plt.ylabel("Pu239/U235")
plt.legend()
plt.show()

"""

prompt_nu_fission = tally.get_values(scores=["prompt-nu-fission"]).sum()
delayed_nu_fission = tally.get_values(scores=["delayed-nu-fission"]).sum()

print(f"Total production of neutrons due to prompt fission: {prompt_nu_fission:.5e}")
print(f"Total production of delayed neutrons due to fission: {delayed_nu_fission:.5e}")

total_neutron_production = prompt_nu_fission + delayed_nu_fission
print(f"Total production of neutrons due to fission: {total_neutron_production:.5e}")

delayed_neutron_fraction = (delayed_nu_fission/total_neutron_production)
print(f"Delayed neutron fraction: {delayed_neutron_fraction:.5e}")
"""

prompt_nu_fission = tally.get_values(scores=["prompt-nu-fission"]).sum()
delayed_nu_fission = tally.get_values(scores=["delayed-nu-fission"]).sum()
total_neutron_production = prompt_nu_fission + delayed_nu_fission
print(f"Total production of neutrons due to fission: {total_neutron_production:.5e}")

E_p = 3.044*10**-11
eta = 2.7
P_0 = 1200*10**6
V = 20172040.76

N=((0.95*P_0*eta)/E_p)/V
print(N)

test =prompt_nu_fission/V
print(test) 