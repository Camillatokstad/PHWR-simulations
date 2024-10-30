import openmc
import numpy as np
import matplotlib.pyplot as plt

#Define materials
urox= openmc.Material(name= "Fuel UO2")
urox.add_nuclide("U235", 0.05)
urox.add_nuclide("U238", 0.95)
urox.add_nuclide("O16", 2)
urox.set_density("g/cm3", 10.98)

helium = openmc.Material(name = "Helium")
helium.add_element("He", 1)
helium.set_density("g/cm3", 0.00000178)

light_water = openmc.Material(name="Water")
light_water.add_element("H",2)
light_water.add_nuclide("O16",1)
light_water.set_density("g/cm3", 1)

heavy_water = openmc.Material(name="Heavy Water")
heavy_water.add_nuclide("H2", 2)
heavy_water.add_nuclide("O16", 1)
heavy_water.set_density("g/cm3", 1.11)

graphite = openmc.Material(name = "Graphite")
graphite.add_nuclide("C0",1)
graphite.set_density("g/cm3",1.86)

zircaloy4 = openmc.Material(name = "Zircaloy")
zircaloy4.add_element("Zr", 0.985)
zircaloy4.add_element("Sn", 0.015)
zircaloy4.set_density("g/cm3", 6.56)

#Define geometry
fuel_pin = -openmc.ZCylinder(r = 0.40)
helium_gap = -openmc.ZCylinder(r=0.42)
cladding_region = -openmc.ZCylinder(r = 0.46)

z_min = openmc.ZPlane(-1, boundary_type = "reflective")
z_max = openmc.ZPlane(1, boundary_type = "reflective")

fuel_pin = fuel_pin & +z_min & -z_max 
helium_gap = helium_gap & +z_min & -z_max
cladding_region = cladding_region & +z_min & -z_max

box = -openmc.model.RectangularPrism(width=1.4, height=1.4, boundary_type = "reflective")
box = box & +z_min & -z_max

helium_gap = helium_gap &~ fuel_pin
cladding_region = cladding_region &~ helium_gap &~ fuel_pin
box = box &~cladding_region &~ helium_gap &~ fuel_pin 

#Create cells
fuel_cell = openmc.Cell(fill=urox, region = fuel_pin)
helium_cell = openmc.Cell(fill=helium, region=helium_gap)
cladding_cell = openmc.Cell(fill=zircaloy4, region=cladding_region)
lightwater_cell = openmc.Cell(fill=light_water, region=box)
heavy_water_cell = openmc.Cell(fill=heavy_water, region=box)
graphite_cell = openmc.Cell(fill=graphite, region=box)


pin_universe = openmc.Universe(name = "Pincell")
#pin_universe.add_cells([fuel_cell,helium_cell,cladding_cell,lightwater_cell])
pin_universe.add_cells([fuel_cell,helium_cell,cladding_cell,heavy_water_cell])
#pin_universe.add_cells([fuel_cell,helium_cell,cladding_cell,graphite_cell])
geometry = openmc.Geometry()
geometry.root_universe = pin_universe 
geometry.export_to_xml()
pin_universe.plot(basis = "xy", pixels = [800,800], color_by = "material",
          colors = {urox: "green", helium: "brown",
zircaloy4: "gray", heavy_water : "blue"}) 
plt.show()
#plt.savefig(fname = "My_Pincell")

#Run simulations
settings = openmc.Settings()
settings.batches = 100 #How many iterations of neutron generations
settings.inactive = 10 #How many generations used to determine sample spots 
settings.particles = 10000 #Number of neutrons in each generation
materials = openmc.Materials([heavy_water, urox, helium, zircaloy4])
openmc.Materials.cross_sections = "endfb-vii.1-hdf5/cross_sections.xml"
materials.export_to_xml()
#The following section tells OpenMC where to sample particle sites.
bounds = [-0.7,-0.7,-1,0.7,0.7,1]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:],only_fissionable=True)
settings.source = openmc.IndependentSource(space=uniform_dist)
settings.export_to_xml()




mesh = openmc.RegularMesh() #Making a mesh instance
mesh.dimension = [100,100] #Setting the bins of the mesh to 100 x 100 bins 
mesh.lower_left= [-0.7,-0.7] #Setting lower left coordinates 
mesh.upper_right= [0.7,0.7] #Setting upper right coordinates
mesh_filter = openmc.MeshFilter(mesh) #Making a meshfilter


tallies_file = openmc.Tallies()

tally = openmc.Tally(name = "Tally")
tally.filters = [mesh_filter]
tally.scores = ["flux", 'prompt-nu-fission', "delayed-nu-fission", "nu-fission"] #Designating the objects to be tallied 
tallies_file.append(tally) #Appending the tally to our tallies file

energy_bins = np.logspace(np.log10(10e-5),np.log10(20e6),501)  # Example energy bins
energy_filter = openmc.EnergyFilter(energy_bins)

tally2 = openmc.Tally(name = "Energy") #Making a new tally only looking at flux 
tally2.filters = [energy_filter] #Designating the energy filter 
tally2.scores = ["flux"]
tallies_file.append(tally2)


# Resonance Escape Probability tallies
therm_abs_rate = openmc.Tally(name='therm. abs. rate')
therm_abs_rate.scores = ['absorption']
therm_abs_rate.filters = [openmc.EnergyFilter([0., 0.625])]
tallies_file.append(therm_abs_rate)

# Thermal Flux Utilization tallies
fuel_therm_abs_rate = openmc.Tally(name='fuel therm. abs. rate')
fuel_therm_abs_rate.scores = ['absorption']
fuel_therm_abs_rate.filters = [openmc.EnergyFilter([0., 0.625]),
                               openmc.CellFilter([fuel_cell])]
tallies_file.append(fuel_therm_abs_rate)
tallies_file.export_to_xml()


if __name__ == "__main__":
        openmc.run()









