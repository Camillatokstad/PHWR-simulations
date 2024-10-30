
import openmc.deplete
import openmc
import numpy as np
import matplotlib.pyplot as plt


pin_volume = np.pi*0.4**2*300
assembly_volume = pin_volume*17*17-(9*pin_volume)
core_volume = 2*(6*assembly_volume)+2*(8*assembly_volume)+2*(10*assembly_volume)+2*(12*assembly_volume)+6*(14*assembly_volume)
print(core_volume)

pin_volume = np.pi*0.7**2*300
assembly_volume = pin_volume*17*17-(9*pin_volume)
whole_core = 2*(6*assembly_volume)+2*(8*assembly_volume)+2*(10*assembly_volume)+2*(12*assembly_volume)+6*(14*assembly_volume)
print(whole_core)
#Define materials
urox= openmc.Material(name= "Fuel UO2")
urox.add_nuclide("U235", 0.05)
urox.add_nuclide("U238", 0.95)
urox.add_nuclide("O16", 2)
urox.set_density("g/cm3", 10.98)
urox.volume = core_volume
urox.temperature = 1200

print(urox)

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
heavy_water.temperature = 600

graphite = openmc.Material(name = "Graphite")
graphite.add_nuclide("C0",1)
graphite.set_density("g/cm3",1.86)

zircaloy4 = openmc.Material(name = "Zircaloy")
zircaloy4.add_element("Zr", 0.985)
zircaloy4.add_element("Sn", 0.015)
zircaloy4.set_density("g/cm3", 6.56)

boron_carbide = openmc.Material(name="Boron Carbide B4C") 
boron_carbide.add_element("B", 1)  # Boron
boron_carbide.add_nuclide("C0", 1)  # Carbon
boron_carbide.set_density("g/cm3", 2.46)

steel = openmc.Material(name = "Steel")
steel.add_nuclide("Fe56",0.98)
steel.add_nuclide("C0",0.02)
steel.set_density("g/cm3", 7.85)

lead = openmc.Material(name= "Lead")
lead.add_nuclide("Pb207", 1)
lead.set_density("g/cm3", 11.33)

#Define geometry
fuel_pin = -openmc.ZCylinder(r = 0.40)
helium_gap = -openmc.ZCylinder(r=0.42)
cladding_region = -openmc.ZCylinder(r = 0.45)

z_min = openmc.ZPlane(-1, boundary_type = "transmission")
z_max = openmc.ZPlane(1, boundary_type = "transmission")

fuel_pin = fuel_pin & +z_min & -z_max 
helium_gap = helium_gap & +z_min & -z_max
cladding_region = cladding_region & +z_min & -z_max

box = -openmc.model.RectangularPrism(width=1.4, height=1.4, boundary_type = "transmission")
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
pin_universe.plot(basis = "xz", pixels = [800,800], color_by = "material",
          colors = {urox: "green", helium: "brown",
zircaloy4: "gray", heavy_water : "blue"}) 
#plt.show()
#plt.savefig(fname = "My_Pincell")

#Run simulations
settings = openmc.Settings()
settings.batches = 100 #How many iterations of neutron generations
settings.inactive = 10 #How many generations used to determine sample spots 
settings.particles = 10000 #Number of neutrons in each generation
materials = openmc.Materials([heavy_water, urox, helium, zircaloy4, boron_carbide,steel, lead])
openmc.Materials.cross_sections = "endfb-vii.1-hdf5/cross_sections.xml"
materials.export_to_xml()
#The following section tells OpenMC where to sample particle sites.
bounds = [-166.6,-166.6,-150,166.6,166.6,150]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:],only_fissionable=True)
settings.source = openmc.IndependentSource(space=uniform_dist)
settings.export_to_xml()





box2 = -openmc.model.RectangularPrism(width=1.4, height=1.4, boundary_type = "transmission")
box2 = box2 & +z_min & -z_max

water_cell2 = openmc.Cell(fill=heavy_water, region=box2)

water_cell_control = openmc.Cell(fill=heavy_water, region=fuel_pin)
water_cell_control.region=openmc.model.rectangular_prism(width=1.4,height=1.4)&~fuel_pin
water_universe = openmc.Universe(name="Water", cells=[water_cell2])

control_rod_cell = openmc.Cell(fill=boron_carbide, region=fuel_pin)
control_universe = openmc.Universe(name="Control Rod Pincell", cells=[control_rod_cell, water_cell_control])

#Constructing geometry
f = pin_universe
c = control_universe
w = water_universe
#Dimensions of the pincells
pinx = 1.4
piny = 1.4
pinz = 2
assembly_x = 17
assembly_y = 17
assembly_z = 150
#Uranium-fuel pin

universii = np.zeros((assembly_z,assembly_x,assembly_y), dtype = openmc.Universe) #Making array for universes
for z in range(0,assembly_z):
    for i in range (0,assembly_x):
        for j in range(0,assembly_y):
            if i in [2,8,14] and j in [2,8,14] and z in np.arange(0):
                universii[z][i][j]=c
            elif i in [2,8,14] and j in [2,8,14] and z in np.arange(0,150):   
                universii[z][i][j] = w
            else:
                universii[z][i][j] = f

f_assembly = openmc.RectLattice(name = "fuel_assembly")
f_assembly.lower_left = (-assembly_x*abs(pinx)/2,-assembly_y*abs(piny)/2,-assembly_z)
f_assembly.universes = universii
f_assembly.pitch = (1.4,1.4,2)
#Making outer universe

water_cell = openmc.Cell(fill = heavy_water)
outer_uni = openmc.Universe(cells = [water_cell])
f_assembly.outer = outer_uni
#Putting the lattice into a simulation
outedge = "transmission" #Infinite reactor or singular assembly
outcap= openmc.ZPlane(assembly_z, boundary_type = outedge)
outbot = openmc.ZPlane(-assembly_z, boundary_type = outedge)
outcell = openmc.model.rectangular_prism(width = 23.8, height = 23.8, boundary_type = outedge) 
outcell = outcell &+outbot &-outcap
f_cell = openmc.Cell(fill = f_assembly, region = outcell) #making cell for assembly 
cell_list = [f_cell]
assembly_uni = openmc.Universe(cells = cell_list) #Putting the assembly into a universe 
geometry = openmc.Geometry()
geometry.root_universe = assembly_uni
#geometry.export_to_xml()

assembly_uni.plot(basis = "xz", pixels = [800,800], color_by = "material",
          colors = {urox: "green", helium: "brown",
zircaloy4: "gray", heavy_water : "blue",  boron_carbide: "red"}) 






dim = 23.8
box_waterass=-openmc.model.RectangularPrism(width=dim, height=dim)& +outbot & -outcap
waterasscell = openmc.Cell(region=box_waterass, fill=heavy_water)
water_u = openmc.Universe(cells=[waterasscell])

a = assembly_uni
e = water_u

assemblies_x = 14
assemblies_y = 14

core = np.zeros((assemblies_x, assemblies_y), dtype = openmc.Universe)

#Top left corner
core[:,:] = a
core[0,:4] = e
core[1,:3] = e
core[2,:2] = e
core[3,:1] = e
#Top right corner
core[0,-4:] = e
core[1,-3:] = e
core[2,-2:] = e
core[3,-1:] = e
#Bottom left
core[-1,:4] = e
core[-2,:3] = e
core[-3,:2] = e
core[-4,:1] = e
#bottom right
core[-1,-4:] = e
core[-2,-3:] = e
core[-3,-2:] = e
core[-4,-1:] = e


core_lattice = openmc.RectLattice(name = "Core lattice ")

core_lattice.lower_left = (-assembly_x*abs(pinx)/2* len(core[1]),-assembly_y*abs(piny)/2* len(core[1]))
core_lattice.universes = core
core_lattice.pitch = (pinx*assembly_x,piny*assembly_y)
core_lattice.outer = outer_uni


coreedge = "transmission"
zmaxcore = openmc.ZPlane(assembly_z, boundary_type = coreedge)
zmincore = openmc.ZPlane(-assembly_z, boundary_type = coreedge)

corecontainment = -openmc.ZCylinder(r = 200, boundary_type = coreedge)
corecontainment = corecontainment &-zmaxcore & +zmincore

core_steel_shield = - openmc.ZCylinder(r = 210, boundary_type = "transmission")
core_steel_z_max = openmc.ZPlane(assembly_z+ 10, boundary_type = "transmission")
core_steel_z_min = openmc.ZPlane(-assembly_z- 10, boundary_type = "transmission")
core_steel_shield = core_steel_shield &+ core_steel_z_min &- core_steel_z_max
core_steel_shield = core_steel_shield &~ corecontainment

core_lead_shield = -openmc.ZCylinder(r = 220, boundary_type = "vacuum")
core_lead_z_max = openmc.ZPlane(assembly_z + 20, boundary_type = "vacuum")
core_lead_z_min = openmc.ZPlane(-assembly_z - 20, boundary_type = "vacuum")
core_lead_shield = core_lead_shield & + core_lead_z_min & - core_lead_z_max
core_lead_shield = core_lead_shield &~ corecontainment &~ core_steel_shield

steel_containment = openmc.Cell(name = "Containment cell")
steel_containment.fill = steel
steel_containment.region = core_steel_shield

lead_shield = openmc.Cell(name = "Lead shielding")
lead_shield.fill = lead
lead_shield.region = core_lead_shield

corecell = openmc.Cell(name = "Core cell")
corecell.fill = core_lattice
corecell.region = corecontainment

core_universe = openmc.Universe(name = "Core_uni", cells = [corecell, steel_containment, lead_shield])

geometry = openmc.Geometry()
geometry.root_universe = core_universe
geometry.export_to_xml()


core_universe.plot(basis = "xy", pixels = [2000,2000], color_by = "material",
          colors = {urox: "green", helium: "brown",
zircaloy4: "gray", heavy_water : "blue",  boron_carbide: "red", lead: "black", steel: "pink"}) 

#plt.show()




mesh = openmc.RegularMesh() #Making a mesh instance
mesh.dimension = [400,400] #Setting the bins of the mesh to 100 x 100 bins 
mesh.lower_left= [-166.6,-166.6] #Setting lower left coordinates 
mesh.upper_right= [166.6,166.6] #Setting upper right coordinates

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

"""
#Designating the model to run depletion on
reactor_model = openmc.Model()
reactor_model.geometry = geometry
reactor_model.materials = materials
reactor_model.settings = settings
reactor_model.tallies = tallies_file 
chain = openmc.deplete.Chain.from_xml("./chain_endfb71_pwr.xml") #Designating the depletion chain being used for our depletion calculation
operator = openmc.deplete.CoupledOperator(reactor_model, "./chain_endfb71_pwr.xml") 
power = 1200e6 #Designating power settings in watts
time_steps = [30 * 24 * 60 * 60] * 6
#Number of timesteps and their length in this cae 6 steps of one month length in seconds.
integrator = openmc.deplete.PredictorIntegrator(operator, time_steps, power) #Running the depletion calculation 



if __name__ == "__main__":
        #openmc.run()
        integrator.integrate()
"""
tallies_file.export_to_xml()

if __name__ == "__main__":
        openmc.run()
