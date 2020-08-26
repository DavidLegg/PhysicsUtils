from utils.quantity import SECOND, LENGTH, ENERGY
from utils.constants import hbar
from utils.graphics import MovingPlot
import numpy as np
import time

def evolve_wavefn_free(psi_0, a, b, m, dt=(1e-4*SECOND), delay=(0*SECOND), N=100, show_components=True):
    '''Generates a moving plot of the free wavefunction over time, in the interval [a, b], for a particle of mass m.
psi_0 should be a function returning a Quantity, in appropriate units.
a and b should be Quantities with units of length.
dt should be a Quantity with units of time.
N is an integer, the number of sample points to evolve with.
delay is a Quantity, the amount of time to sleep between updates.'''
    delay = delay.to(SECOND).value
    
    step_quantity = (b - a) / (N - 1)
    X_quantity = [a + i*step_quantity for i in range(N)]
    Y_quantity = [psi_0(x).standardize() for x in X_quantity]
    if Y_quantity[0].dimension() != LENGTH**(-1/2):
        raise ValueError('psi_0 should return Quantity with dimension {}, not {}'.format(LENGTH**(-1/2), Y_quantity[0].dimension()))

    step = float( step_quantity.standardize().value )
    multiplier = complex( (dt * 1j * hbar / (2 * m)).standardize().value ) / (step ** 2)
    X_display = np.array([float(x.value) for x in X_quantity])
    X = np.array([float(x.standardize().value) for x in X_quantity])
    Y = np.array([complex(y.value) for y in Y_quantity])
    Y /= ( np.abs(Y[:-1])**2 * np.diff(X) ).sum()**(1/2)
    if show_components:
        plot = MovingPlot(X_display, np.abs(Y), np.real(Y), np.imag(Y))
        plot.lines[1].set_color('blue')
        plot.lines[1].set_linewidth(1)
        plot.lines[2].set_color('red')
        plot.lines[2].set_linewidth(1)
    else:
        plot = MovingPlot(X_display, np.abs(Y))
    plot.lines[0].set_color('black')
    plot.lines[0].set_linewidth(2)
    plot.lines[0].set_zorder(10)
    plot.ax.set_xlabel('Position ({})'.format(step_quantity.unit))
    plot.ax.set_ylabel('Wavefunction Magnitude ({})'.format(Y_quantity[0].unit))

    t = 0*dt
    plot.ax.set_title('Evolving Wavefunction... (t = {})'.format(t))
    # Temporary for efficiently calculating second-derivatives
    change = np.zeros(Y.shape, complex)
    while True:
        time.sleep(delay)
        # Implements the free Schrodinger equation as the update step
        change[1:-1] = (Y[:-2] + Y[2:] - 2*Y[1:-1])
        change[0]    = change[1] / 2
        change[-1]   = change[-2] / 2
        Y += multiplier * change
        # Re-normalize; not theoretically necessary, but helps avoid piling-up of error in practice
        Y /= (( np.abs(Y[:-1])**2 ).sum() * step)**(1/2)
        t += dt
        plot.ax.set_title('Evolving Wavefunction... (t = {})'.format(t))
        if show_components:
            plot.batch_update(np.abs(Y), np.real(Y), np.imag(Y))
        else:
            plot.update(np.abs(Y))

def evolve_wavefn(psi_0, V, a, b, m, dt=(1e-4*SECOND), delay=(0*SECOND), N=100, show_components=True):
    '''Generates a moving plot of the wavefunction over time, under potential V, in the interval [a, b], for a particle of mass m.
psi_0 and V should be functions returning a Quantity, in appropriate units.
a and b should be Quantities with units of length.
dt should be a Quantity with units of time.
N is an integer, the number of sample points to evolve with.
delay is a Quantity, the amount of time to sleep between updates.'''
    delay = delay.to(SECOND).value
    
    step_quantity = (b - a) / (N - 1)
    X_quantity = [a + i*step_quantity for i in range(N)]
    Y_quantity = [psi_0(x).standardize() for x in X_quantity]
    if Y_quantity[0].dimension() != LENGTH**(-1/2):
        raise ValueError('psi_0 should return Quantity with dimension {}, not {}'.format(LENGTH**(-1/2), Y_quantity[0].dimension()))
    if V(a).dimension() != ENERGY:
        raise ValueError('V should return Quantity with dimension {}, not {}'.format(ENERGY, V(a).dimension()))

    step = float( step_quantity.standardize().value )
    free_multiplier = complex( (dt * 1j * hbar / (2 * m)).standardize().value ) / (step ** 2)
    V_quantity = [V(x) for x in X_quantity]
    V_multiplier = np.array([complex( (dt * v * (-1j) / hbar).standardize().value ) for v in V_quantity])

    X_display = np.array([float(x.value) for x in X_quantity])
    X = np.array([float(x.standardize().value) for x in X_quantity])
    Y = np.array([complex(y.value) for y in Y_quantity])
    Y /= ( np.abs(Y[:-1])**2 * np.diff(X) ).sum()**(1/2)

    # Re-scale potential for display:
    V_display  = np.array([float(v.value) for v in V_quantity])
    V_display -= V_display.min()
    V_display = np.abs(Y).max() * (2 * V_display / V_display.max() - 1)

    if show_components:
        plot = MovingPlot(X_display, np.abs(Y), np.real(Y), np.imag(Y), V_display)
        plot.lines[1].set_color('blue')
        plot.lines[1].set_linewidth(1)
        plot.lines[2].set_color('red')
        plot.lines[2].set_linewidth(1)
        plot.lines[3].set_color('grey')
        plot.lines[3].set_linewidth(1)
        plot.lines[3].set_zorder(-10)
    else:
        plot = MovingPlot(X_display, np.abs(Y), V_display)
        plot.lines[1].set_color('grey')
        plot.lines[1].set_linewidth(1)
    plot.lines[0].set_color('black')
    plot.lines[0].set_linewidth(2)
    plot.lines[0].set_zorder(10)
    plot.ax.set_xlabel('Position ({})'.format(step_quantity.unit))
    plot.ax.set_ylabel('Wavefunction Magnitude ({})'.format(Y_quantity[0].unit))

    t = 0*dt
    plot.ax.set_title('Evolving Wavefunction... (t = {})'.format(t))
    # Temporary for efficiently calculating second-derivatives
    change = np.zeros(Y.shape, complex)
    error_acc = np.zeros(Y.shape, complex)
    while True:
        time.sleep(delay)
        # Implements the Schrodinger equation as the update step
        change[1:-1] = (Y[:-2] + Y[2:] - 2*Y[1:-1])
        change[0]    = change[1] / 2
        change[-1]   = change[-2] / 2

        # Conceptually, want this update:
        # Y += free_multiplier * change + V_multiplier * Y
        # Converted to Kahan sum for higher precision:
        change = free_multiplier * change + V_multiplier * Y - error_acc
        Y_new = Y + change
        error_acc = (Y_new - Y) - change
        Y = Y_new

        # Re-normalize; not theoretically necessary, but helps avoid piling-up of error in practice
        Y /= (( np.abs(Y[:-1])**2 ).sum() * step)**(1/2)
        t += dt
        plot.ax.set_title('Evolving Wavefunction... (t = {})'.format(t))
        if show_components:
            plot.batch_update(np.abs(Y), np.real(Y), np.imag(Y))
        else:
            plot.update(np.abs(Y))
