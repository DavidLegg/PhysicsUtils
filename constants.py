import math
from utils.quantity import *

class Constant(Quantity):
    '''Defines a physics constant with units, and optionally a canonical symbol (variable) and a description.'''
    _ALL_CONSTANTS = []

    def __init__(self, value, symbol='', desc='', unit=None):
        '''Constructs a constant. If unit is None, deduces unit from value, which should be a Quantity, or becomes a scalar.'''
        super().__init__(value, unit)
        assert isinstance(symbol, str), 'Symbol must be a string'
        assert isinstance(desc, str), 'Description must be a string'
        self._sym  = symbol
        self._desc = desc
        Constant._ALL_CONSTANTS.append(self)

    def description(self):
        '''A natural text description of the constant.'''
        # "if not empty" formatting function
        ine = lambda s, fmt: (fmt.format(s) if s else '')
        return ine(self._sym, '{} = ') + super().__str__() + ine(self._desc, ' ({})')

    @staticmethod
    def lookup(query):
        '''Does a fuzzy search against all known constants.'''
        regex = re.compile('.*'.join(query.split()))
        return [x for x in Constant._ALL_CONSTANTS if regex.search(x.description())]

def find_constant(query):
    for c in Constant.lookup(query):
        print(c.description())

pi    = Constant(math.pi, 'π', 'π (unitless constant)')
c     = Constant(299792458, 'c', 'Speed of light', METER / SECOND)
e     = Constant(1.60217662e-19, 'e', 'Electron charge', COULOMB)
m_e   = Constant(9.10938356e-31, 'm_e', 'Electron mass', KILOGRAM)
m_p   = Constant(1.6726219e-27, 'm_p', 'Proton mass', KILOGRAM)
mu_0  = Constant(1.2566370621219e-6, 'μ_0', 'Permeability of free space', HENRY / METER)
eps_0 = Constant(1 / (mu_0 * c**2), 'ε_0', 'Permittivity of free space')
k     = Constant(1/(4*pi*eps_0), 'k', 'Coulomb constant')
h     = Constant(6.62607004e-34, 'h', 'Planck constant', KILOGRAM * METER**2 / SECOND)
hbar  = Constant(h / (2 * pi), 'ħ', 'Reduced Planck constant')
alpha = Constant(k * e**2 / (hbar * c), 'α', 'Fine-structure constant')
a_0   = Constant(hbar / (m_e * c * alpha), 'a_0', 'Bohr radius')
