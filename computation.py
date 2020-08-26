from utils.quantity import Quantity

def sqrt(value):
    '''Square root that works for Unit or Quantity as well.'''
    return value**0.5

def kahan_sum(iterable, start = 0):
    '''Low-error Kahan Summation. Taken from https://en.wikipedia.org/wiki/Kahan_summation_algorithm'''
    total = start
    c = 0*start
    for x in iterable:
        y = x - c
        t = total + y
        c = (t - total) - y
        total = t
    return total

def integrate(f, a, b, N = 30, eps = 1e-4):
    '''Integral from a to b of f(x) wrt x, calculated adaptively to maintain an estimated relative error < eps, and using an initial rough estimate to avoid unnecessarily precise work.'''

    fa = f(a)
    step = (b - a) / (N - 1)

    points = [(a + i*step, f(a + i*step), 1) for i in range(N - 1, -1, -1)]

    # Use a simple trapezoidal approximation with no subdivisions as a rough estimate of the integral
    abs_eps = eps * abs((points[0][1] + points[-1][1]) / 2 + sum((y for (_, y, _) in points[1:-1]), 0*fa)) / N

    def riemann_estimates():
        while len(points) > 1:
            x1,y1,w1 = points.pop()
            x2,y2,w2 = points[-1]
            # Sample a test point between x1, x2, and do intermediate calculations
            xt  = (x1 + x2) / 2
            yt  = f(xt)
            wt  = 2 * max(w1, w2)
            y12 = y1 + y2
            ytt = 2 * yt
            ytotal = y12 + ytt
            while abs(y12 - ytt) > eps * abs(ytotal) + abs_eps: # unacceptable error
                # Make the test point the new x2
                x2 = xt
                y2 = yt
                w2 = wt
                points.append((x2,y2,w2))
                # Sample a test point between x1, x2, and redo intermediate calculations
                xt  = (x1 + x2) / 2
                yt  = f(xt)
                wt *= 2
                y12 = y1 + y2
                ytt = 2 * yt
                ytotal = y12 + ytt
            # Error is acceptable, yield the Riemann estimate
            yield (ytotal / 2) / wt
    return kahan_sum(riemann_estimates(), 0*fa) * step


