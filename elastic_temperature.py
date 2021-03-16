
from argparse import ArgumentParser
from lammps import PyLammps

def potential(lmp, args):
    """ set up potential and minimization """
    ff_string = ' '
    ff_string = ff_string.join(args.elements) # merge all element string to one string
    lmp.kim_interactions(ff_string)

    # Setup neighbor style
    lmp.neighbor(1.0, "nsq")
    lmp.neigh_modify("once no every 1 delay 0 check yes")

    # reset the current timestep to zero
    lmp.reset_timestep(0)

    # Setup MD
    lmp.timestep(lmp.variables["timestep"].value)
    lmp.fix(4, "all", "nve")
    if lmp.variables["thermostat"].value == 1:
        lmp.fix(5, "all", "langevin", lmp.variables["temp"].value, lmp.variables["temp"].value, lmp.variables["tdamp"].value, lmp.variables["seed"].value)

    lmp.fix("avp", "all", "ave/time", lmp.variables["nevery"].value, lmp.variables["nrepeat"].value, lmp.variables["nfreq"].value, "c_thermo_press", "mode vector")

    # Setup output
    lmp.thermo(args.nthermo)
    lmp.thermo_style("custom step temp pe press f_avp[1] f_avp[2] f_avp[3] f_avp[4] f_avp[5] f_avp[6]")
    lmp.thermo_modify("norm no")

    return

def displace(lmp, args, idir):
    """computes the response to a small strain """

    if idir == 1:
        lmp.variable("len0 equal {}".format(lmp.variables["lx0"].value))
    elif idir == 2 or idir == 6:
        lmp.variable("len0 equal {}".format(lmp.variables["ly0"].value))
    else:
        lmp.variable("len0 equal {}".format(lmp.variables["lz0"].value))

    # Reset box and simulation parameters
    lmp.clear()
    lmp.box("tilt large")
    lmp.kim_init(args.kim_model, "metal", "unit_conversion_mode")
    lmp.read_restart("restart.equil")
    lmp.change_box("all triclinic")
    potential(lmp, args)

    # Negative deformation
    lmp.variable("delta equal -${up}*${len0}")
    lmp.variable("deltaxy equal -${up}*xy")
    lmp.variable("deltaxz equal -${up}*xz")
    lmp.variable("deltayz equal -${up}*yz")

    if idir == 1:
        lmp.change_box("all x delta 0 ${delta} xy delta ${deltaxy} xz delta ${deltaxz} remap units box")
    elif idir == 2:
        lmp.change_box("all y delta 0 ${delta} yz delta ${deltayz} remap units box")
    elif idir == 3:
        lmp.change_box("all z delta 0 ${delta} remap units box")
    elif idir == 4:
        lmp.change_box("all yz delta ${delta} remap units box")
    elif idir == 5:
        lmp.change_box("all xz delta ${delta} remap units box")
    else:
        lmp.change_box("all xy delta ${delta} remap units box")

    # Run MD
    lmp.run(lmp.variables["nequil"].value)
    lmp.run(lmp.variables["nrun"].value)

    # Obtain new stress tensor
    lmp.variable("pxx1 equal f_avp[1]")
    lmp.variable("pyy1 equal f_avp[2]")
    lmp.variable("pzz1 equal f_avp[3]")
    lmp.variable("pxy1 equal f_avp[4]")
    lmp.variable("pxz1 equal f_avp[5]")
    lmp.variable("pyz1 equal f_avp[6]")

    # Compute elastic constant from pressure tensor
    C1neg = lmp.variables["d1"].value
    C2neg = lmp.variables["d2"].value
    C3neg = lmp.variables["d3"].value
    C4neg = lmp.variables["d4"].value
    C5neg = lmp.variables["d5"].value
    C6neg = lmp.variables["d6"].value

    # Reset box and simulation parameters
    lmp.clear()
    lmp.box("tilt large")
    lmp.kim_init(args.kim_model, "metal", "unit_conversion_mode")
    lmp.read_restart("restart.equil")
    lmp.change_box("all triclinic")
    potential(lmp, args)

    # Positive deformation
    lmp.variable("delta equal ${up}*${len0}")
    lmp.variable("deltaxy equal ${up}*xy")
    lmp.variable("deltaxz equal ${up}*xz")
    lmp.variable("deltayz equal ${up}*yz")

    if idir == 1:
        lmp.change_box("all x delta 0 ${delta} xy delta ${deltaxy} xz delta ${deltaxz} remap units box")
    elif idir == 2:
        lmp.change_box("all y delta 0 ${delta} yz delta ${deltayz} remap units box")
    elif idir == 3:
        lmp.change_box("all z delta 0 ${delta} remap units box")
    elif idir == 4:
        lmp.change_box("all yz delta ${delta} remap units box")
    elif idir == 5:
        lmp.change_box("all xz delta ${delta} remap units box")
    else:
        lmp.change_box("all xy delta ${delta} remap units box")

    # Run MD
    lmp.run(lmp.variables["nequil"].value)
    lmp.run(lmp.variables["nrun"].value)

    # Obtain new stress tensor
    lmp.variable("pxx1 equal f_avp[1]")
    lmp.variable("pyy1 equal f_avp[2]")
    lmp.variable("pzz1 equal f_avp[3]")
    lmp.variable("pxy1 equal f_avp[4]")
    lmp.variable("pxz1 equal f_avp[5]")
    lmp.variable("pyz1 equal f_avp[6]")

    # Compute elasic constant from pressure tensor
    C1pos = lmp.variables["d1"].value
    C2pos = lmp.variables["d2"].value
    C3pos = lmp.variables["d3"].value
    C4pos = lmp.variables["d4"].value
    C5pos = lmp.variables["d5"].value
    C6pos = lmp.variables["d6"].value

    # Combine posiive and negative
    lmp.variable("C1{} equal {}".format(idir, 0.5*(C1neg+C1pos)))
    lmp.variable("C2{} equal {}".format(idir, 0.5*(C2neg+C2pos)))
    lmp.variable("C3{} equal {}".format(idir, 0.5*(C3neg+C3pos)))
    lmp.variable("C4{} equal {}".format(idir, 0.5*(C4neg+C4pos)))
    lmp.variable("C5{} equal {}".format(idir, 0.5*(C5neg+C5pos)))
    lmp.variable("C6{} equal {}".format(idir, 0.5*(C6neg+C6pos)))

    return

def elastic():
    """ Compute elastic constant tensor for a crystal at finite temperature

     Written by Aidan Thompson (Sandia, athomps@sandia.gov)

     Global constants:
       up = the deformation magnitude (strain units)

     To run this on a different system, it should only be necessary to
     modify the files init.mod and potential.mod. In order to calculate
     the elastic constants correctly, care must be taken to specify
     the correct units in init.mod (units, cfac and cunits). It is also
     important to verify that the MD sampling of stress components
     is generating accurate statistical averages.
     One indication of this is that the elastic constants are insensitive
     to the choice of the variable ${up} in init.mod. Another is to  check for finite size effects. """

    parser = ArgumentParser(description='A python script to compute elastic properties of bulk materials')

    parser.add_argument("input_data_file", help="The full path & name of the lammps data file.")
    parser.add_argument("kim_model", help="the KIM ID of the interatomic model archived in OpenKIM")
    parser.add_argument("elements", nargs='+', default=['Au'], help="a list of N chemical species, which defines a mapping between atom types in LAMMPS to the available species in the OpenKIM model")
    parser.add_argument("adiabatic", default=False, help="determines if the computations is in NVE or NVT ensemble")
    args = parser.parse_args()

    L = PyLammps()

    L.units("metal")

    # Define the finite deformation size. Try several values to verify that results do not depend on it.
    L.variable("up equal 2.0e-2")

    # Define MD parameters
    L.variable("nevery equal 10")                  # sampling interval
    L.variable("nrepeat equal 10")                 # number of samples
    L.variable("nfreq equal ${nevery}*${nrepeat}") # length of one average
    L.variable("nthermo equal ${nfreq}")           # interval for thermo output
    L.variable("nequil equal 10*${nthermo}")       # length of equilibration run
    L.variable("nrun equal 3*${nthermo}")          # length of equilibrated run
    L.variable("temp equal 2000.0")                # temperature of initial sample
    L.variable("timestep equal 0.001")             # timestep
    L.variable("tdamp equal 0.01")                 # time constant for thermostat
    L.variable("seed equal 123457")                # seed for thermostat

    # metal units, elastic constants in GPa
    cfac = 1.0e-4

    L.boundary("p", "p", "p") # periodic boundary conditions in all three directions
    L.box("tilt large") # to avoid termination if the final simulation box has a high tilt factor

    # use the OpenKIM model to set the energy interactions
    L.kim_init(args.kim_model, "metal", "unit_conversion_mode")

    L.read_data(args.input_data_file)

    potential(L, args)

    # Need to set mass to something, just to satisfy LAMMPS
    mass_dictionary = {'H': 1.00797, 'He': 4.00260, 'Li': 6.941, 'Be': 9.01218, 'B': 10.81, 'C': 12.011, 'N': 14.0067, 'O': 15.9994, 'F': 18.998403, 'Ne': 20.179, 'Na': 22.98977, 'Mg': 24.305, 'Al': 26.98154, 'Si': 28.0855, 'P': 30.97376, 'S': 32.06, 'Cl': 35.453, 'K': 39.0983, 'Ar': 39.948, 'Ca': 40.08, 'Sc': 44.9559, 'Ti': 47.90, 'V': 50.9415, 'Cr': 51.996, 'Mn': 54.9380, 'Fe': 55.847, 'Ni': 58.70, 'Co': 58.9332, 'Cu': 63.546, 'Zn': 65.38, 'Ga': 69.72, 'Ge': 72.59, 'As': 74.9216, 'Se': 78.96, 'Br': 79.904, 'Kr': 83.80, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.9059, 'Zr': 91.22, 'Nb': 92.9064, 'Mo': 95.94, 'Tc': (98), 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.4, 'Ag': 107.868, 'Cd': 112.41, 'In': 114.82, 'Sn': 118.69, 'Sb': 121.75, 'I': 126.9045, 'Te': 127.60, 'Xe': 131.30, 'Cs': 132.9054, 'Ba': 137.33, 'La': 138.9055, 'Ce': 140.12, 'Pr': 140.9077, 'Nd': 144.24, 'Pm': (145), 'Sm': 150.4, 'Eu': 151.96, 'Gd': 157.25, 'Tb': 158.9254, 'Dy': 162.50, 'Ho': 164.9304, 'Er': 167.26, 'Tm': 168.9342, 'Yb': 173.04, 'Lu': 174.967, 'Hf': 178.49, 'Ta': 180.9479, 'W': 183.85, 'Re': 186.207, 'Os': 190.2, 'Ir': 192.22, 'Pt': 195.09, 'Au': 196.9665, 'Hg': 200.59, 'Tl': 204.37, 'Pb': 207.2, 'Bi': 208.9804, 'Po': (209), 'At': (210), 'Rn': (222), 'Fr': (223), 'Ra': 226.0254, 'Ac': 227.0278, 'Pa': 231.0359, 'Th': 232.0381, 'Np': 237.0482, 'U': 238.029}
    for itype in range(1, len(args.elements)+1):
        L.mass(itype, mass_dictionary.get(args.elements[itype-1], 1.0e-20))

    # Compute initial state at zero pressure
    L.run(L.variables["nequil"].value)

    if not args.adiabatic:
        L.variable("thermostat equal 0")
    else:
        L.variable("thermostat equal 1")

    print L.variables("thermostat").eval

    L.run(L.variables["nrun"].value)

    L.variable("lx0 equal {}".format(L.eval("lx")))
    L.variable("ly0 equal {}".format(L.eval("ly")))
    L.variable("lz0 equal {}".format(L.eval("lz")))

    # These formulas define the derivatives w.r.t. strain components
    L.variable("d1 equal -(v_pxx1-{})/(v_delta/v_len0)*{}".format(L.eval("f_avp[1]"), cfac))
    L.variable("d2 equal -(v_pyy1-{})/(v_delta/v_len0)*{}".format(L.eval("f_avp[2]"), cfac))
    L.variable("d3 equal -(v_pzz1-{})/(v_delta/v_len0)*{}".format(L.eval("f_avp[3]"), cfac))
    L.variable("d4 equal -(v_pyz1-{})/(v_delta/v_len0)*{}".format(L.eval("f_avp[4]"), cfac))
    L.variable("d5 equal -(v_pxz1-{})/(v_delta/v_len0)*{}".format(L.eval("f_avp[5]"), cfac))
    L.variable("d6 equal -(v_pxy1-{})/(v_delta/v_len0)*{}".format(L.eval("f_avp[6]"), cfac))

    # Write restart
    L.write_restart("restart.equil")

    for idir in range(1, 7):
        displace(L, args, idir)

    postprocess_and_output(L)
    return

def postprocess_and_output(lmp):
    """Compute the moduli and print everything to screen """

    # Output final values
    C11all = lmp.variables["C11"].value
    C22all = lmp.variables["C22"].value
    C33all = lmp.variables["C33"].value

    C12all = 0.5*(lmp.variables["C12"].value + lmp.variables["C21"].value)
    C13all = 0.5*(lmp.variables["C13"].value + lmp.variables["C31"].value)
    C23all = 0.5*(lmp.variables["C23"].value + lmp.variables["C32"].value)

    C44all = lmp.variables["C44"].value
    C55all = lmp.variables["C55"].value
    C66all = lmp.variables["C66"].value

    C14all = 0.5*(lmp.variables["C14"].value + lmp.variables["C41"].value)
    C15all = 0.5*(lmp.variables["C15"].value + lmp.variables["C51"].value)
    C16all = 0.5*(lmp.variables["C16"].value + lmp.variables["C61"].value)

    C24all = 0.5*(lmp.variables["C24"].value + lmp.variables["C42"].value)
    C25all = 0.5*(lmp.variables["C25"].value + lmp.variables["C52"].value)
    C26all = 0.5*(lmp.variables["C26"].value + lmp.variables["C62"].value)

    C34all = 0.5*(lmp.variables["C34"].value + lmp.variables["C43"].value)
    C35all = 0.5*(lmp.variables["C35"].value + lmp.variables["C53"].value)
    C36all = 0.5*(lmp.variables["C36"].value + lmp.variables["C63"].value)

    C45all = 0.5*(lmp.variables["C45"].value + lmp.variables["C54"].value)
    C46all = 0.5*(lmp.variables["C46"].value + lmp.variables["C64"].value)
    C56all = 0.5*(lmp.variables["C56"].value + lmp.variables["C65"].value)

    # Average moduli for cubic crystals
    C11cubic = (C11all + C22all + C33all)/3.0
    C12cubic = (C12all + C13all + C23all)/3.0
    C44cubic = (C44all + C55all + C66all)/3.0

    bulkmodulus = (C11cubic + 2*C12cubic)/3.0
    shearmodulus1 = C44cubic
    shearmodulus2 = (C11cubic - C12cubic)/2.0
    poisson_ratio = 1.0/(1.0 + C11cubic/C12cubic)

    # print results to screen
    print("=========================================")
    print("Components of the Elastic Constant Tensor")
    print("=========================================")

    print("Elastic Constant C11all = {} GPa".format(C11all))
    print("Elastic Constant C22all = {} GPa".format(C22all))
    print("Elastic Constant C33all = {} GPa".format(C33all))

    print("Elastic Constant C12all = {} GPa".format(C12all))
    print("Elastic Constant C13all = {} GPa".format(C13all))
    print("Elastic Constant C23all = {} GPa".format(C23all))

    print("Elastic Constant C44all = {} GPa".format(C44all))
    print("Elastic Constant C55all = {} GPa".format(C55all))
    print("Elastic Constant C66all = {} GPa".format(C66all))

    print("Elastic Constant C14all = {} GPa".format(C14all))
    print("Elastic Constant C15all = {} GPa".format(C15all))
    print("Elastic Constant C16all = {} GPa".format(C16all))

    print("Elastic Constant C24all = {} GPa".format(C24all))
    print("Elastic Constant C25all = {} GPa".format(C25all))
    print("Elastic Constant C26all = {} GPa".format(C26all))

    print("Elastic Constant C34all = {} GPa".format(C34all))
    print("Elastic Constant C35all = {} GPa".format(C35all))
    print("Elastic Constant C36all = {} GPa".format(C36all))

    print("Elastic Constant C45all = {} GPa".format(C45all))
    print("Elastic Constant C46all = {} GPa".format(C46all))
    print("Elastic Constant C56all = {} GPa".format(C56all))

    print("=========================================")
    print("Average properties for a cubic crystal")
    print("=========================================")

    print("Bulk Modulus = {} GPa".format(bulkmodulus))
    print("Shear Modulus 1 = {} GPa".format(shearmodulus1))
    print("Shear Modulus 2 = {} GPa".format(shearmodulus2))
    print("Poisson Ratio = {}".format(poisson_ratio))

    return

if __name__ == "__main__":
    elastic()
