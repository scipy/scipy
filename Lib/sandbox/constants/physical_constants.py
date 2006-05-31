# Compiled by Charles Harris
# Taken from his email message to scipy-dev 
#  dated October 3, 2002
""" Fundamental Physical Constants

    These constants are taken from CODATA Recommended Values of the
    Fundamental Physical Constants: 1998. They may be found at
    physics.nist.gov/constants. The values are stored in the dictionary
    physical_constants as a tuple containing the value, the units, and
    the relative precision, in that order. All constants are in SI units
    unless otherwise stated.

    Several helper functions are provided:

    value(key) returns the value of the physical constant.
    unit(key) returns the units of the physical constant.
    precision(key) returns the relative precision of the physical constant.
    find(sub) prints out a list of keys containing the string sub. 
"""

physical_constants = {
        
    # UNIVERSAL

    'speed of light in vacuum': (2.99792458e8,'m s-1',0.0),
    'magnetic constant': (1.25663706143591730e-6,'N A-2',0.0),
    'electric constant': (8.85418781762038920e-12,'F m-1',0.0),
    'characteristic impedance of vacuum': (3.76730313461770660e2,'Ohm',0.0),
    'Newtonian constant of gravitation': (6.673e-11,'m3 kg-1 s-2',1.5e-3),
    'Planck constant': (6.62606876e-34,'J s',7.8e-8),
    'Planck constant in eV s': (4.13566727e-15,'eV s',3.9e-8),
    'Planck constant divided by 2 pi': (1.054571596e-34,'J s',7.8e-8),
    'Planck constant divided by 2 pi in eV s': (6.58211889e-16,'eV s',3.9e-8),
    'Planck mass': (2.1767e-8,'kg',7.5e-4),
    'Planck length': (1.6160e-35,'m',7.5e-4),
    'Planck time': (5.3906e-44,'s',7.5e-4),

    # ELECTROMAGNETIC

    'elementary charge': (1.602176462e-19,'C',3.9e-8),
    'elementary charge divided by h': (2.417989491e14,'A J-1',3.9e-8),
    'magnetic flux quantum': (2.067833636e-15,'Wb',3.9e-8),
    'conductance quantum': (7.748091696e-5,'S',3.7e-9),
    'inverse of conductance quantum': (12906.403786,'Ohm',3.7e-9),
    'Josephson constant': (483597.898e9,'Hz V-1',3.9e-8),
    'von Klitzing constant': (25812.807572,'Ohm',3.7e-9),
    'Bohr magneton': (927.400899e-26,'J T-1',4.0e-8),
    'Bohr magneton in eV T-1': (5.788381749e-5,'eV T-1',7.3e-9),
    'Bohr magneton divided by h': (13.99624624e9,'Hz T-1',4.0e-8),
    'Bohr magneton divided by hc': (46.6864521,'m-1 T-1',4.0e-8),
    'Bohr magneton divided by k': (0.6717131,'K T-1',1.7e-6),
    'nuclear magneton': (5.05078317e-27,'J T-1',4.0e-8),
    'nuclear magneton in eV T-1': (3.152451238e-8,'eV T-1',7.6e-9),
    'nuclear magneton divided by h': (7.62259396,'MHz T-1',4.0e-8),
    'nuclear magneton divided by hc': (2.54262366e-2,'m-1 T-1',4.0e-8),
    'nuclear magneton divided by k': (3.6582638e-4,'K T-1',1.7e-6),

    # ATOMIC AND NUCLEAR

    # General

    'fine-structure constant': (7.297352533e-3,'',3.7e-9),
    'inverse fine-structure constant': (137.03599976,'',3.7e-9),
    'Rydberg constant': (10973731.568549,'m-1',7.6e-12),
    'Rydberg constant times c': (3.289841960368e15,'Hz',7.6e-12),
    'Rydberg constant times hc': (2.17987190e-18,'J',7.8e-8),
    'Rydberg constant times hc in eV': (13.60569172,'eV',3.9e-8),
    'Bohr radius': (0.5291772083e-10,'m',3.7e-9),
    'Hartree energy': (4.35974381e-18,'J',7.8e-8),
    'Hartree energy in eV': (27.2113834,'eV',3.9e-8),
    'quantum of circulation': (3.636947516e-4,'m2 s-1',7.3e-9),

    # Electroweak

    'Fermi coupling constant': (1.16639e-5,'GeV-2',8.6e-6),
    'weak mixing angle': (0.2224,'',8.7e-3),

    # Electron, e-

    'electron mass': (9.10938188e-31,'kg',7.9e-8),
    'electron mass in u': (5.485799110e-4,'u',2.1e-9),
    'electron mass energy equivalent': (8.18710414e-14,'J',7.9e-8),
    'electron mass in MeV': (0.510998902,'MeV',4.0e-8),
    'electron-muon mass ratio': (4.83633210e-3,'',3.0e-8),
    'electron-tau mass ratio': (2.87555e-4,'',1.6e-4),
    'electron-proton mass ratio': (5.446170232e-4,'',2.1e-9),
    'electron-neutron mass ratio': (5.438673462e-4,'',2.2e-9),
    'electron-deuteron mass ratio': (2.7244371170e-4,'',2.1e-9),
    'electron to alpha particle mass ratio': (1.3709335611e-4,'',2.1e-9),
    'electron charge to mass quotient': (-1.758820174e11,'C kg-1',4.0e-8),
    'electron molar mass': (5.485799110e-7,'kg mol-1',2.1e-9),
    'Compton wavelength': (2.426310215e-12,'m',7.3e-9),
    'Compton wavelength divided by 2 pi': (386.1592642e-15,'m',7.3e-9),
    'classical electron radius': (2.817940285e-15,'m',1.1e-8),
    'Thomson cross section': (0.665245854e-28,'m2',2.2e-8),
    'electron magnetic moment': (-928.476362e-26,'J T-1',4.0e-8),
    'electron magnetic moment to Bohr magneton ratio': (-1.0011596521869,'',4.1e-12),
    'electron magnetic moment to nuclear magneton ratio': (-1838.2819660,'',2.1e-9),
    'electron magnetic moment anomaly': (1.1596521869e-3,'',3.5e-9),
    'electron g-factor': (-2.0023193043737,'',4.1e-12),
    'electron-muon magnetic moment ratio': (206.7669720,'',3.0e-8),
    'electron-proton magnetic moment ratio': (-658.2106875,'',1.0e-8),
    'electron to shielded proton magnetic moment ratio': (-658.2275954,'',1.1e-8),
    'electron-neutron magnetic moment ratio': (960.92050,'',2.4e-7),
    'electron-deuteron magnetic moment ratio': (-2143.923498,'',1.1e-8),
    'electron to shielded helion magnetic moment ratio': (864.058255,'',1.2e-8),
    'electron gyromagnetic ratio': (1.760859794e11,'s-1 T-1',4.0e-8),
    'electron gyromagnetic ratio divided by 2 pi': (2.80249540e10,'s-1 T-1',4.0e-8),

    # Muon, µ-

    'muon mass': (1.88353109e-28,'kg',8.4e-8),
    'muon mass in u': (0.1134289168,'u',3.0e-8),
    'muon mass energy equivalent': (1.69283332e-11,'J',8.4e-8),
    'muon mass in MeV': (105.6583568,'MeV',4.9e-8),
    'muon-electron mass ratio': (206.7682657,'',3.0e-8),
    'muon-tau mass ratio': (5.94572e-2,'',1.6e-4),
    'muon-proton mass ratio': (0.1126095173,'',3.0e-8),
    'muon-neutron mass ratio': (0.1124545079,'',3.0e-8),
    'muon molar mass': (0.1134289168e-3,'kg mol-1',3.0e-8),
    'muon Compton wavelength': (11.73444197e-15,'m',2.9e-8),
    'muon Compton wavelength divided by 2 pi': (1.867594444e-15,'m',2.9e-8),
    'muon magnetic moment': (-4.49044813e-26,'J T-1',4.9e-8),
    'muon magnetic moment to Bohr magneton ratio': (-4.84197085e-3,'',3.0e-8),
    'muon magnetic moment to nuclear magneton ratio': (-8.89059770,'',3.0e-8),
    'muon magnetic moment anomaly': (1.16591602e-3,'',5.5e-7),
    'muon g-factor': (-2.0023318320,'',6.4e-10),
    'muon-proton magnetic moment ratio': (-3.18334539,'',3.2e-8),

    # Tau

    'tau mass': (3.16788e-27,'kg',1.6e-4),
    'tau mass in u': (1.90774,'u',1.6e-4),
    'tau mass energy equivalent': (2.84715e-10,'J',1.6e-4),
    'tau mass in MeV': (1777.05,'MeV',1.6e-4),
    'tau-electron mass ratio': (3477.60,'',1.6e-4),
    'tau-muon mass ratio': (16.8188,'',1.6e-4),
    'tau-proton mass ratio': (1.89396,'',1.6e-4),
    'tau-neutron mass ratio': (1.89135,'',1.6e-4),
    'tau molar mass': (1.90774e-3,'kg mol-1',1.6e-4),
    'tau Compton wavelength': (0.69770e-15,'m',1.6e-4),
    'tau Compton wavelength divided by 2 pi': (0.111042e-15,'m',1.6e-4),

    # Proton, p

    'proton mass': (1.67262158e-27,'kg',7.9e-8),
    'proton mass in u': (1.00727646688,'u',1.3e-10),
    'proton mass energy equivalent': (1.50327731e-10,'J',7.9e-8),
    'proton mass in MeV': (938.271998,'MeV',4.0e-8),
    'proton-electron mass ratio': (1836.1526675,'',2.1e-9),
    'proton-muon mass ratio': (8.88024408,'',3.0e-8),
    'proton-tau mass ratio': (0.527994,'',1.6e-4),
    'proton-neutron mass ratio': (0.99862347855,'',5.8e-10),
    'proton charge to mass quotient': (9.57883408e7,'C kg-1',4.0e-8),
    'proton molar mass': (1.00727646688e-3,'kg mol-1',1.3e-10),
    'proton Compton wavelength': (1.321409847e-15,'m',7.6e-9),
    'proton Compton wavelength divided by 2 pi': (0.2103089089e-15,'m',7.6e-9),
    'proton magnetic moment': (1.410606633e-26,'J T-1',4.1e-8),
    'proton magnetic moment to Bohr magneton ratio': (1.521032203e-3,'',1.0e-8),
    'proton magnetic moment to nuclear magneton ratio': (2.792847337,'',1.0e-8),
    'proton g-factor': (5.585694675,'',1.0e-8),
    'proton-neutron magnetic moment ratio': (-1.45989805,'',2.4e-7),
    'shielded proton magnetic moment': (1.410570399e-26,'J T-1',4.2e-8),
    'shielded proton magnetic moment to Bohr magneton ratio': (1.520993132e-3,'',1.1e-8),
    'shielded proton magnetic moment to nuclear magneton ratio': (2.792775597,'',1.1e-8),
    'proton magnetic shielding correction': (25.687e-6,'',5.7e-4),
    'proton gyromagnetic ratio': (2.67522212e8,'s-1 T-1',4.1e-8),
    'shielded proton gyromagnetic ratio': (2.67515341e8,'s-1 T-1',4.2e-8),

    # Neutron, n

    'neutron mass': (1.67492716e-27,'kg',7.9e-8),
    'neutron mass in u': (1.00866491578,'u',5.4e-10),
    'neutron mass energy equivalent': (1.50534946e-10,'J',7.9e-8),
    'neutron mass energy equivalent in MeV': (939.565330,'MeV',4.0e-8),
    'neutron-electron mass ratio': (1838.6836550,'',2.2e-9),
    'neutron-muon mass ratio': (8.89248478,'',3.0e-8),
    'neutron-tau mass ratio': (0.528722,'',1.6e-4),
    'neutron-proton mass ratio': (1.00137841887,'',5.8e-10),
    'neutron molar mass': (1.00866491578e-3,'kg mol-1',5.4e-10),
    'neutron Compton wavelength': (1.319590898e-15,'m',7.6e-9),
    'neutron Compton wavelength divided by 2 pi': (0.2100194142e-15,'m',7.6e-9),
    'neutron magnetic moment': (-0.96623640e-26,'J T-1',2.4e-7),
    'neutron magnetic moment to Bohr magneton ratio': (-1.04187563e-3,'',2.4e-7),
    'neutron magnetic moment to nuclear magneton ratio': (-1.91304272,'',2.4e-7),
    'neutron g-factor': (-3.82608545,'',2.4e-7),
    'neutron-electron magnetic moment ratio': (1.04066882e-3,'',2.4e-7),
    'neutron-proton magnetic moment ratio': (-0.68497934,'',2.4e-7),
    'neutron to shielded proton magnetic moment ratio': (-0.68499694,'',2.4e-7),
    'neutron gyromagnetic ratio': (1.83247188e8,'s-1 T-1',2.4e-7),
    'neutron gyromagnetic ratio n/': (29.1646958,'MHz T-1',2.4e-7),

    # Deuteron, d

    'deuteron mass': (3.34358309e-27,'kg',7.9e-8),
    'deuteron mass in u': (2.01355321271,'u',1.7e-10),
    'deuteron mass energy equivalent': (3.00506262e-10,'J',7.9e-8),
    'deuteron mass in MeV': (1875.612762,'MeV',4.0e-8),
    'deuteron-electron mass ratio': (3670.4829550,'',2.1e-9),
    'deuteron-proton mass ratio': (1.99900750083,'',2.0e-10),
    'deuteron molar mass': (2.01355321271e-3,'kg mol-1',1.7e-10),
    'deuteron magnetic moment': (0.433073457e-26,'J T-1',4.2e-8),
    'deuteron magnetic moment to Bohr magneton ratio': (0.4669754556e-3,'',1.1e-8),
    'deuteron magnetic moment to nuclear magneton ratio': (0.8574382284,'',1.1e-8),
    'deuteron-electron magnetic moment ratio': (-4.664345537e-4,'',1.1e-8),
    'deuteron-proton magnetic moment ratio': (0.3070122083,'',1.5e-8),
    'deuteron-neutron magnetic moment ratio': (-0.44820652,'',2.4e-7),

    # Helion, h

    'helion mass': (5.00641174e-27,'kg',7.9e-8),
    'helion mass in u': (3.01493223469,'u',2.8e-10),
    'helion mass energy equivalent': (4.49953848e-10,'J',7.9e-8),
    'helion mass in MeV': (2808.39132,'MeV',4.0e-8),
    'helion-electron mass ratio': (5495.885238,'',2.1e-9),
    'helion-proton mass ratio': (2.99315265850,'',3.1e-10),
    'helion molar mass': (3.01493223469e-3,'kg mol-1',2.8e-10),
    'shielded helion magnetic moment': (-1.074552967e-26,'J T-1',4.2e-8),
    'shielded helion magnetic moment to Bohr magneton ratio': (-1.158671474e-3,'',1.2e-8),
    'shielded helion magnetic moment to nuclear magneton ratio': (-2.127497718,'',1.2e-8),
    'shielded helion to proton magnetic moment ratio': (-0.761766563,'',1.5e-8),
    'shielded helion to shielded proton magnetic moment ratio': (-0.7617861313,'',4.3e-9),
    'shielded helion gyromagnetic ratio': (2.037894764e8,'s-1 T-1',4.2e-8),
    'shielded helion gyromagnetic ratio': (32.4341025,'MHz T-1',4.2e-8),

    # Alpha particle

    'alpha particle mass': (6.64465598e-27,'kg',7.9e-8),
    'alpha particle mass in u': (4.0015061747,'u',2.5e-10),
    'alpha particle mass energy equivalent': (5.97191897e-10,'J',7.9e-8),
    'alpha particle mass in MeV': (3727.37904,'MeV',4.0e-8),
    'alpha particle to electron mass ratio': (7294.299508,'',2.1e-9),
    'alpha particle to proton mass ratio': (3.9725996846,'',2.8e-10),
    'alpha particle molar mass': (4.0015061747e-3,'kg mol-1',2.5e-10),

    # PHYSICO-CHEMICAL

    'Avogadro constant': (6.02214199e23,'mol-1',7.9e-8),
    'atomic mass constant': (1.66053873e-27,'kg',7.9e-8),
    'atomic mass constant energy equivalent': (1.49241778e-10,'J',7.9e-8),
    'atomic mass constant in MeV': (931.494013,'MeV',4.0e-8),
    'Faraday constant': (96485.3415,'C mol-1',4.0e-8),
    'molar Planck constant': (3.990312689e-10,'J s mol-1',7.6e-9),
    'molar Planck constant times c': (0.11962656492,'J m mol-1',7.6e-9),
    'molar gas constant': (8.314472,'J mol-1 K-1',1.7e-6),
    'Boltzmann constant': (1.3806503e-23,'J K-1',1.7e-6),
    'Boltzmann constant in eV K-1': (8.617342e-5,'eV K-1',1.7e-6),
    'Boltzmann constant divided by h': (2.0836644e10,'Hz K-1',1.7e-6),
    'Boltzmann constant divided by hc': (69.50356,'m-1 K-1',1.7e-6),
    'molar volume of ideal gas': (22.413996e-3,'m3 mol-1',1.7e-6),
    'Loschmidt constant': (2.6867775e25,'m-3',1.7e-6),
    'Sackur-Tetrode constant': (-1.1648678,'',3.7e-6),
    'Stefan-Boltzmann constant': (5.670400e-8,'W m-2 K-4',7.0e-6),
    'first radiation constant': (3.74177107e-16,'W m2',7.8e-8),
    'first radiation constant for spectral radiance': (1.191042722e-16,'W m2 sr-1',7.8e-8),
    'second radiation constant': (1.4387752e-2,'m K',1.7e-6),
    'Wien displacement law constant': (2.8977686e-3,'m K',1.7e-6),

    # Non-SI Units

    'electron volt': (1.602176462e-19,'J',3.9e-8),
    'atomic mass unit': (1.66053873e-27,'kg',7.9e-8),
    'standard temperature': (273.15,'K',0.0),
    'standard pressure': (1.01325e5,'Pa',0.0),

    # Natural Units

    'n.u. of velocity': (2.99792458e8,'m s-1',0.0),
    'n.u. of action': (1.054571596e-34,'J s',7.8e-8),
    'n.u. of action in eV s': (6.58211889e-16,'eV s',3.9e-8),
    'n.u. of mass': (9.10938188e-31,'kg',7.9e-8),
    'n.u. of energy': (8.18710414e-14,'J',7.9e-8),
    'n.u. of energy in MeV': (0.510998902,'MeV',4.0e-8),
    'n.u. of momentum': (2.73092398e-22,'kg m s-1',7.9e-8),
    'n.u. of momentum in MeV/c': (0.510998902,'MeV/c',4.0e-8),
    'n.u. of length': (386.1592642e-15,'m',7.3e-9),
    'n.u. of time': (1.2880886555e-21,'s',7.3e-9),

    # Atomic Units

    'a.u. of charge': (1.602176462e-19,'C',3.9e-8),
    'a.u. of mass': (9.10938188e-31,'kg',7.9e-8),
    'a.u. of action': (1.054571596e-34,'J s',7.8e-8),
    'a.u. of length': (0.5291772083e-10,'m',3.7e-9),
    'a.u. of energy': (4.35974381e-18,'J',7.8e-8),
    'a.u. of time': (2.418884326500e-17,'s',7.6e-12),
    'a.u. of force': (8.23872181e-8,'N',7.8e-8),
    'a.u. of velocity': (2.1876912529e6,'m s-1',3.7e-9),
    'a.u. of momentum': (1.99285151e-24,'kg m s-1',7.8e-8),
    'a.u. of current': (6.62361753e-3,'A',3.9e-8),
    'a.u. of charge density': (1.081202285e12,'C m-3',4.0e-8),
    'a.u. of electric potential': (27.2113834,'V',3.9e-8),
    'a.u. of electric field': (5.14220624e11,'Vm-1',3.9e-8),
    'a.u. of electric field gradient': (9.71736153e21,'Vm-2',4.0e-8),
    'a.u. of electric dipole moment': (8.47835267e-30,'C m',3.9e-8),
    'a.u. of electric quadrupole moment': (4.48655100e-40,'C m2',4.0e-8),
    'a.u. of electric polarizability': (1.648777251e-41,'C2 m2 J-1',1.1e-8),
    'a.u. of 1st hyperpolarizability': (3.20636157e-53,'C3 m3 J-2',4.2e-8),
    'a.u. of 2nd hyperpolarizability': (6.23538112e-65,'C4 m4 J-3',8.1e-8),
    'a.u. of magnetic flux density': (2.350517349e5,'T',4.0e-8),
    'a.u. of magnetic dipole moment': (1.854801799e-23,'J T-1',4.0e-8),
    'a.u. of magnetizability': (7.89103641e-29,'J T-2',1.8e-8),
    'a.u. of permittivity': (1.112650056e-10,'F m-1',0.0)
    }

import string

def value(key) :
    """value indexed by key"""
    return physical_constants[key][0]

def unit(key) :
    """unit indexed by key"""
    return physical_constants[key][1]

def precision(key) :
    """relative precision indexed by key"""
    return physical_constants[key][2]

def find(sub) :
    """list all keys containing the string sub"""
    l_sub = string.lower(sub)
    result = []
    for key in physical_constants :
        l_key = string.lower(key)
        if string.find(l_key,l_sub) >= 0 :
            result.append(key)
    result.sort()
    for key in result :
        print key
