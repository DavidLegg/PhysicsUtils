from collections import defaultdict
from utils.number import ScientificNumber

class Dimension:
    '''The dimension that a Unit can belong to, e.g. Length or Time.

Uses a collection of base dimensions, and allows arbitrary composition of dimensions through normal mathematical operators.'''
    def __init__(self):
        '''Initializes a new dimension. For internal use only. Reference pre-defined dimensions instead of declaring new ones.'''
        self._components = defaultdict(int)

    def _copy(self):
        output = Dimension()
        output._components = self._components.copy()
        return output

    @staticmethod
    def base(symbol):
        '''Declare a new base dimension. This dimension is incommensurate with all other base dimensions.'''
        output = Dimension()
        output._components[symbol] = 1
        return output

    def components(self):
        '''Yields each base Dimension, and the power that it's raised to.'''
        d = Dimension()
        for k,v in self._components.items():
            d._components[k] = 1
            yield (d, v)
            del d._components[k]

    def standard_unit(self):
        for u in Unit._FOR_DIM.get(self, []):
            if u._si_factor == 1:
                return u
        # No pre-defined standard unit. Compose base units:
        u = SCALAR
        for k,v in self.components():
            u *= k.standard_unit() ** v
        return u

    def __abs__(self):
        return self

    def __neg__(self):
        return self

    def __pos__(self):
        return self

    def __add__(self, other):
        assert isinstance(other, Dimension), 'Cannot add Dimension and {}'.format(type(other))
        assert self == other, 'Cannot add unlike Dimensions {} and {}'.format(self, other)
        return self

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self + other

    def __rsub__(self, other):
        return self + other

    def __mul__(self, other):
        assert isinstance(other, Dimension), 'Cannot multiply Dimension and {}'.format(type(other))
        output = self._copy()
        for k,v in other._components.items():
            output._components[k] += v
        return output

    def __rmul__(self, other):
        # Dimension multiplication is symmetric
        return self * other

    def __floordiv__(self, other):
        return self / other

    def __rfloordiv__(self, other):
        return other / self

    def __truediv__(self, other):
        assert isinstance(other, Dimension), 'Cannot divide Dimension and {}'.format(type(other))
        output = self._copy()
        for k,v in other._components.items():
            output._components[k] -= v
        return output

    def __rtruediv__(self, other):
        assert isinstance(other, Dimension), 'Cannot divide Dimension and {}'.format(type(other))
        output = other._copy()
        for k,v in self._components.items():
            output._components[k] -= v
        return output

    def __pow__(self, exp):
        assert isinstance(exp, int) or isinstance(exp, float), 'Cannot raise Dimension to a {}'.format(type(other))
        output = self._copy()
        for k in self._components.keys():
            output._components[k] *= exp
        return output

    def __eq__(self, other):
        if not isinstance(other, Dimension):
            return False
        self._trim_components()
        other._trim_components()
        return self._components == other._components

    def __neq__(self, other):
        return not (self == other)

    def __str__(self):
        self._trim_components()
        if self._components:
            return ' '.join((k if v == 1 else '{}^{}'.format(k, v)) for k,v in sorted(self._components.items()))
        else:
            return 'SCALAR_DIM'

    def __repr__(self):
        return str(self)

    def __hash__(self):
        self._trim_components()
        return hash(tuple(sorted(self._components)))

    def _trim_components(self):
        for k,v in list(self._components.items()):
            if v == 0:
                del self._components[k]
            elif isinstance(v, float) and v.is_integer():
                self._components[k] = int(v)

SCALAR_DIM    = Dimension()
TIME          = Dimension.base('T')
LENGTH        = Dimension.base('L')
MASS          = Dimension.base('M')
CURRENT       = Dimension.base('I')
TEMPERATURE   = Dimension.base('Temp')
AMOUNT        = Dimension.base('N')
LUM_INTENSITY = Dimension.base('J')

VELOCITY = LENGTH / TIME
ACCELERATION = VELOCITY / TIME
FORCE = MASS * ACCELERATION
ENERGY = FORCE * LENGTH
POWER = ENERGY / TIME
MOMENTUM = MASS * VELOCITY
FREQUENCY = TIME**(-1)
AREA = LENGTH**2
VOLUME = LENGTH**3
PRESSURE = FORCE / AREA
CHARGE = CURRENT * TIME
VOLTAGE = POWER / CURRENT
CAPACITANCE = CHARGE / VOLTAGE
RESISTANCE = VOLTAGE / CURRENT
MAGNETIC_FLUX = ENERGY / CURRENT

class Unit:
    '''The Unit for a Quantity.

Comprises base Units, corresponding to base Dimensions, and derived Units defined in terms of base Units.
Units can be composed using standard mathematical operations, though errors will be raised if Units of different Dimensions are combined in inappropriate ways.'''

    # Dictionary from dimension to standard units for that dimension
    _PREFIXES = []
    _FOR_DIM  = defaultdict(list)

    def __init__(self):
        '''Initializes a new unit. For internal use only. Reference pre-defined dimensions instead of declaring new ones.'''
        self._short_name = None
        self._long_name  = None
        self.dimension  = SCALAR_DIM
        self._si_factor  = 1

    @staticmethod
    def base(short_name, long_name, base_dimension):
        '''Define a base unit for base_dimension.'''
        unit = Unit()
        unit._short_name = short_name
        unit._long_name  = long_name
        unit.dimension   = base_dimension
        unit._si_factor  = 1
        Unit._FOR_DIM[base_dimension].append(unit)
        return unit

    @staticmethod
    def define(short_name, long_name, quantity, make_standard=True):
        '''Define a new named Unit, equaling one of Quantity.

    The Dimension for this Unit is derived automatically from quantity.
    If make_standard is True, will set this named Unit as the standard Unit for that Dimension, used for automatic naming.'''
        # Force quantity to be a Quantity, even if it comes in as a number or Unit instead
        quantity = Quantity(quantity, None)
        unit = Unit()
        unit._short_name = short_name
        unit._long_name  = long_name
        unit.dimension   = quantity.unit.dimension
        unit._si_factor  = quantity.value * quantity.unit._si_factor
        if make_standard:
            if unit.dimension == SCALAR_DIM:
                Unit._PREFIXES.append((unit._si_factor, unit))
            else:
                Unit._FOR_DIM[unit.dimension].append(unit)
        return unit

    def descale(self):
        '''Returns scaling coefficient for this unit, and a version of this unit without scaling coefficient.'''
        unit = Unit()
        unit.dimension = self.dimension
        return (self._si_factor, unit)

    def autoname(self):
        '''Picks a name for this unit if one can be automatically determined.'''
        def find_prefix(value):
            found = False
            for v,p in Unit._PREFIXES:
                if v.tolerant_equal(value):
                    return p
            return None

        if self.dimension == SCALAR_DIM:
            p = find_prefix(self._si_factor)
            if p is not None:
                self._short_name = p._short_name
                self._long_name  = p._long_name
        else:
            for u in Unit._FOR_DIM[self.dimension]:
                q = Quantity(1, self).to(u)
                if q.value == 1:
                    self._short_name = u._short_name
                    self._long_name  = u._long_name
                    break
                p = find_prefix(q.value)
                if p is not None:
                    self._short_name = p._short_name + u._short_name
                    self._long_name  = p._long_name + '-' + u._long_name
                    break

    def __abs__(self):
        return self

    def __neg__(self):
        return self

    def __pos__(self):
        return self

    def __add__(self, other):
        assert isinstance(other, Unit), 'Cannot add/subtract Unit and {}'.format(type(other))
        assert self == other, 'Cannot add/subtract unlike Units {} and {}'.format(self, other)
        return self

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self + other

    def __rsub__(self, other):
        return self + other

    def __mul__(self, other):
        if isinstance(other, Unit):
            unit = Unit()
            unit.dimension = self.dimension * other.dimension
            unit._si_factor = self._si_factor * other._si_factor
            return unit
        else:
            return other * Quantity(1, self)

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        if isinstance(other, Unit):
            return self * (other ** -1)
        else:
            return Quantity(1, self) / other

    def __rtruediv__(self, other):
        if isinstance(other, Unit):
            return other / self
        else:
            return other / Quantity(1, self)

    def __pow__(self, exp):
        assert isinstance(exp, int) or isinstance(exp, float), 'Cannot raise Unit to a {}'.format(type(other))
        unit = Unit()
        unit.dimension = self.dimension ** exp
        unit._si_factor = self._si_factor ** exp
        return unit

    def __eq__(self, other):
        return isinstance(other, Unit) and self.dimension == other.dimension and self._si_factor == other._si_factor

    def __neq__(self, other):
        return not (self == other)

    def __str__(self):
        if self._short_name:
            return self._short_name
        elif self._long_name:
            return self._long_name

        # if name isn't yet given, attempt to determine a name
        self.autoname()
        if self._short_name:
            return self._short_name
        elif self._long_name:
            return self._long_name

        s = self._derivation()
        return (s if s else '(scalar)')

    def __repr__(self):
        return str(self)

    def _derivation(self):
        x = []
        for k,v in self.dimension.components():
            if v == 0:
                continue
            x.append((str(k.standard_unit()), v))

        x.sort()
        return (('' if self._si_factor == 1 else 'x {} '.format(self._si_factor))
                + ' '.join(str(k if v == 1 else '{}^{}'.format(k, v)) for k,v in x))

    def description(self):
        '''Returns a long-form description of the unit, including the dimensions and derivation if appropriate.'''
        # "if not empty" formatting function
        ine = lambda s, fmt: (fmt.format(s) if s else '')

        derivation = self._derivation()
        if not (self._short_name or self._long_name):
            # If anonymous, attempt to name
            self.autoname()

        if not (self._short_name or self._long_name):
            # Anonymous unit, so name it by the derivation:
            return 'Unit ' + derivation
        elif self._long_name:
            # Fully named unit
            return 'Unit ' + self._long_name + ine(self._short_name, ' [{}]') + ' = ' + derivation
        else:
            # Only a short name
            return 'Unit ' + self._short_name + ' = ' + derivation

# Need to define a scalar unit for Quantity to function properly
SCALAR = Unit()

class Quantity:
    '''Describes a numerical amount in a specified Unit.

Quantities can participate in normal mathematical operations, and will calculate units and scaling factors automatically.
Errors will be raised if incommensurate Units are combined inappropriately.'''
    def __init__(self, value, unit):
        if unit is None:
            value = Quantity._force(value)
            self.value = ScientificNumber(value.value)
            self.unit  = value.unit
        else:
            self.value = ScientificNumber(value)
            self.unit  = unit

    @staticmethod
    def _force(thing):
        if isinstance(thing, Quantity):
            return thing
        if isinstance(thing, int) or isinstance(thing, float) or isinstance(thing, complex):
            return Quantity(thing, SCALAR)
        if isinstance(thing, Unit):
            return Quantity(1, thing)
        raise ValueError('{} is not a Quantity, int or float'.format(thing))

    def standardize(self):
        '''Converts this Quantity to standard SI base units.'''
        scale, unit = self.unit.descale()
        return Quantity(self.value * scale, unit)

    def to(self, unit):
        '''Returns this quantity, re-expressed in given units.'''
        assert self.unit.dimension == unit.dimension, 'Cannot convert between unlike dimensions {} and {}'.format(self.unit.dimension, unit.dimension)
        standard = self.standardize()
        return Quantity(standard.value / unit._si_factor, unit)

    def dimension(self):
        return self.unit.dimension

    def __abs__(self):
        return Quantity(abs(self.value), abs(self.unit))

    def __neg__(self):
        return Quantity(-self.value, -self.unit)

    def __pos__(self):
        return Quantity(+self.value, +self.unit)

    def __add__(self, other):
        other = Quantity._force(other).to(self.unit)
        return Quantity(self.value + other.value, self.unit + other.unit)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return other + (-self)

    def __mul__(self, other):
        other = Quantity._force(other)
        return Quantity(self.value * other.value, self.unit * other.unit)

    def __rmul__(self, other):
        return self * other

    def __floordiv__(self, other):
        other = Quantity._force(other)
        return Quantity(self.value // other.value, self.unit // other.unit)

    def __rfloordiv__(self, other):
        other = Quantity._force(other)
        return Quantity(other.value // self.value, other.unit // self.unit)

    def __truediv__(self, other):
        other = Quantity._force(other)
        return Quantity(self.value / other.value, self.unit / other.unit)

    def __rtruediv__(self, other):
        other = Quantity._force(other)
        return Quantity(other.value / self.value, other.unit / self.unit)

    def __pow__(self, exp):
        assert isinstance(exp, int) or isinstance(exp, float), 'Cannot raise Quantity to a {}'.format(type(other))
        return Quantity(self.value ** exp, self.unit ** exp)

    def __eq__(self, other):
        other = Quantity._force(other).to(self.unit)
        return self.value == other.value

    def __neq__(self, other):
        return not (self == other)

    def __lt__(self, other):
        other = Quantity._force(other).to(self.unit)
        return self.value < other.value

    def __le__(self, other):
        other = Quantity._force(other).to(self.unit)
        return self.value <= other.value

    def __gt__(self, other):
        other = Quantity._force(other).to(self.unit)
        return self.value > other.value

    def __ge__(self, other):
        other = Quantity._force(other).to(self.unit)
        return self.value >= other.value

    def __int__(self):
        return int(self.to(SCALAR).value)

    def __float__(self):
        return float(self.to(SCALAR).value)

    def __complex__(self):
        return complex(self.to(SCALAR).value)

    def conjugate(self):
        return Quantity(self.value.conjugate(), self.unit)

    def __str__(self):
        return '{} {}'.format(self.value, self.unit)

    def __repr__(self):
        return str(self)

# SI Base units

SECOND   = Unit.base('s',   'second',   TIME)
METER    = Unit.base('m',   'meter',    LENGTH)
KILOGRAM = Unit.base('kg',  'kilogram', MASS)
AMPERE   = Unit.base('A',   'ampere',   CURRENT)
KELVIN   = Unit.base('K',   'kelvin',   TEMPERATURE)
MOLE     = Unit.base('mol', 'mole',     AMOUNT)
CANDELA  = Unit.base('cd',  'candela',  LUM_INTENSITY)

# Common SI Derived Units, for convenience

HERTZ   = Unit.define('Hz', 'hertz',   1 / SECOND)
NEWTON  = Unit.define('N',  'newton',  KILOGRAM * METER / SECOND**2)
PASCAL  = Unit.define('Pa', 'pascal',  NEWTON / METER**2)
JOULE   = Unit.define('J',  'joule',   NEWTON * METER)
WATT    = Unit.define('W',  'watt',    JOULE / SECOND)
COULOMB = Unit.define('C',  'coulomb', AMPERE * SECOND)
VOLT    = Unit.define('V',  'volt',    JOULE / COULOMB)
FARAD   = Unit.define('F',  'farad',   COULOMB / VOLT)
OHM     = Unit.define('Ω',  'ohm',     VOLT / AMPERE)
SIEMENS = Unit.define('S',  'siemens', 1 / OHM)
WEBER   = Unit.define('Wb', 'weber',   JOULE / AMPERE)
TESLA   = Unit.define('T',  'tesla',   VOLT * SECOND / METER**2)
HENRY   = Unit.define('H',  'henry',   OHM * SECOND)
LUMEN   = Unit.define('lm', 'lumen',   CANDELA)
LUX     = Unit.define('lx', 'lux',     CANDELA / METER**2)

# Prefixes

YOTTA = Unit.define('Y',  'yotta', 1e24)
ZETTA = Unit.define('Z',  'zetta', 1e21)
EXA   = Unit.define('E',  'exa',   1e18)
PETA  = Unit.define('P',  'peta',  1e15)
TERA  = Unit.define('T',  'tera',  1e12)
GIGA  = Unit.define('G',  'giga',  1e9)
MEGA  = Unit.define('M',  'mega',  1e6)
KILO  = Unit.define('k',  'kilo',  1e3)
HECTO = Unit.define('h',  'hecto', 1e2)
DEKA  = Unit.define('da', 'deka',  1e1)
DECI  = Unit.define('d',  'deci',  1e-1)
CENTI = Unit.define('c',  'centi', 1e-2)
MILLI = Unit.define('m',  'milli', 1e-3)
MICRO = Unit.define('mu', 'micro', 1e-6)
NANO  = Unit.define('n',  'nano',  1e-9)
PICO  = Unit.define('p',  'pico',  1e-12)
FEMTO = Unit.define('f',  'femto', 1e-15)
ATTO  = Unit.define('a',  'atto',  1e-18)
ZEPTO = Unit.define('z',  'zepto', 1e-21)
YOCTO = Unit.define('y',  'yocto', 1e-24)

# Extra non-SI units and their conversion factors

ELECTRONVOLT = Unit.define('eV', 'electronvolt', 1.602176634e-19*JOULE)
AMU          = Unit.define('u', 'atomic mass unit', 1.6605390666050e-27*KILOGRAM)
DALTON       = Unit.define('Da', 'dalton', AMU)
ANGSTROM     = Unit.define('Å', 'angstrom', 1e-10*METER)
