SCI_DECIMALS = 2

def set_precision(n):
    global SCI_DECIMALS
    assert isinstance(n, int) and n >= 0, "Precision must be a non-negative integer"
    SCI_DECIMALS = n

def sci(number):
    return '{:.{d}e}'.format(number, d=SCI_DECIMALS)

class _Scientific:
    EQUAL_THRESH = 1e-12

    # Reverse factory method: Converts _Scientific numbers to regular ones
    @staticmethod
    def _undo(x):
        if isinstance(x, complex):
            return complex(x)
        else:
            return float(x)
    
    # True overrides

    def __init__(self, over_type):
        self._over_type = over_type

    def tolerant_equal(self, value):
        return float.__lt__(abs(self - value), abs(self * ScientificFloat.EQUAL_THRESH))

    def approx(self, value, eps=1e-6):
        return float.__lt__(abs(self - value), abs(self * eps))


    def __eq__(self, value):
        return self.tolerant_equal(value)

    def __neq__(self, value):
        return not self.tolerant_equal(value)

    def __lt__(self, value):
        return self._over_type.__lt__(self._over_type(self), _Scientific._undo(value)) and not self.tolerant_equal(value)

    def __gt__(self, value):
        return self._over_type.__gt__(self._over_type(self), _Scientific._undo(value)) and not self.tolerant_equal(value)

    def __le__(self, value):
        return self._over_type.__le__(self._over_type(self), _Scientific._undo(value)) or self.tolerant_equal(value)

    def __ge__(self, value):
        return self._over_type.__ge__(self._over_type(self), _Scientific._undo(value)) or self.tolerant_equal(value)

    # Type up-conversions, that propagate the Scientific types through calculations

    def __abs__(self):
        return ScientificNumber( super().__abs__() )

    def __add__(self, value):
        return ScientificNumber(self._over_type(self) + value )

    def __divmod__(self, value):
        return ScientificNumber( divmod(self._over_type(self), value) )

    def __floordiv__(self, value):
        return ScientificNumber( self._over_type(self) // value )

    def __mod__(self, value):
        return ScientificNumber( self._over_type(self) % value )

    def __mul__(self, value):
        return ScientificNumber( self._over_type(self) * value )

    def __neg__(self):
        return ScientificNumber( -self._over_type(self) )

    def __pos__(self):
        return ScientificNumber( +self._over_type(self) )

    def __pow__(self, value, mod=None):
        return ScientificNumber( pow(self._over_type(self), value, mod) )

    def __radd__(self, value):
        return ScientificNumber( value + self._over_type(self) )

    def __rdivmod__(self, value):
        return ScientificNumber( divmod(value, self._over_type(self)) )

    def __rfloordiv__(self, value):
        return ScientificNumber( value // self._over_type(self) )

    def __rmod__(self, value):
        return ScientificNumber( value % self._over_type(self) )

    def __rmul__(self, value):
        return ScientificNumber( value * self._over_type(self) )

    def __rpow__(self, value, mod=None):
        return ScientificNumber( pow(value, self._over_type(self), mod) )

    def __rsub__(self, value):
        return ScientificNumber( value - self._over_type(self) )

    def __rtruediv__(self, value):
        return ScientificNumber( value / self._over_type(self) )

    def __sub__(self, value):
        return ScientificNumber( self._over_type(self) - value )

    def __truediv__(self, value):
        return ScientificNumber( self._over_type(self) / value )

    def conjugate(*args, **kwargs):
        return ScientificNumber( self._over_type(self).conjugate(*args, **kwargs) )

class ScientificFloat(_Scientific, float):
    '''Floating point number that prints itself in scientific notation, and uses tolerant equality.'''

    # Actual overrides:
    def __new__(self, value):
        return float.__new__(self, value)

    def __init__(self, value):
        _Scientific.__init__(self, float)
        float.__init__(self)

    def __str__(self):
        return sci(self)

    def __repr__(self):
        return sci(self)

    def fromhex(*args, **kwargs):
        return ScientificNumber( super().fromhex(*args, **kwargs) )

class ScientificComplex(_Scientific, complex):
    '''Complex floating point number that prints itself in scientific notation, and uses tolerant equality.'''

    # Actual overrides:
    def __new__(self, value):
        return complex.__new__(self, value)

    def __init__(self, value):
        _Scientific.__init__(self, complex)
        complex.__init__(self)

    def __str__(self):
        return '({} + {}j)'.format(sci(self.real), sci(self.imag))

    def __repr__(self):
        return str(self)

def ScientificNumber(x):
    if isinstance(x, complex):
        return ScientificComplex(x)
    else:
        return ScientificFloat(x)
