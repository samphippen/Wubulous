from mfi import mfi

import gmpy
from time import time
import random

itercount = 0
tests     = 0

class ecc:
    def __init__(self, a4, a6, mod):
        self.a4 = mfi(a4, mod)
        self.a6 = mfi(a6, mod)
        self.mod = mod        

    def is_quadratic_residue(self, x):
        l = x.jacobi()
        return l == 1

    def is_valid_point(self, x, y):
        x = mfi(x, self.mod)
        y = mfi(y, self.mod)
        if x == 0 and y == 1:
            return True
        l = y**2
        r = x**3 + x*self.a4 + self.a6
        return l == r

    def get_point(self, n):
        x = mfi(n, self.mod)
        while x != 0:
            t = (x ** 3) + (self.a4 * x) + self.a6
            if t == 0:
                return [(x, mfi(0, self.mod))]
            elif self.is_quadratic_residue(t):
                vals = t.modsqrt()
                points = []
                for k in vals:
                    points.append((x,k))
                return points
            x += 1
        return []

    def get_random_point(self):
        r = gmpy.mpz(random.random() * self.mod)
        return random.choice(self.get_point(r))

    def get_points_list(self):
        points = [(0,1)]
        x = mfi(1, self.mod)
        while x != 0:
            p = self.get_point(x)
            points += p
            x = p[0][0] + 1
        return points
   
    def negate(self, x,y):
        x = mfi(x, self.mod)
        y = mfi(y, self.mod)
        return x, y + x

    def modinv(self, n):
        n = mfi(n, self.mod)
        return n.modinv()

    def double(self, x1, y1):
        x1 = mfi(x1, self.mod)
        y1 = mfi(y1, self.mod)

        if (y1 == 1 and x1 == 0):
            return 0,1
        if (y1 + y1) == 0:
            return 0,1

        mi = self.modinv(y1 * 2)
        a4 = mfi(self.a4, self.mod)
        ld = (x1*x1*3 + a4) * mi
        x3 = ld * ld - x1 - x1
        y3 = ld * (x1 - x3) - y1
        return x3, y3

    def add(self, x1, y1, x2, y2):
        ld = 0
        x1 = mfi(x1, self.mod)
        y1 = mfi(y1, self.mod)
        x2 = mfi(x2, self.mod)
        y2 = mfi(y2, self.mod)

        if (y1 == 1 and x1 == 0):
            return x2, y2
        elif (x2 == 0 and y2 == 1):
            return x1,y1

        if x1 == x2:
            if y1 == y2:
                return self.double(x1,y1)
            else:
                return (0,1)

        else:
            ld = (y2 - y1) * (x2-x1).modinv()
            x3 = ld * ld - x2 - x1
            y3 = ld * (x1 - x3) - y1
            return x3, y3

    def multiply(self, x1, y1, factor):
        accum = (x1, y1)
        for i in xrange(0, factor - 1):
            accum = self.add(accum[0], accum[1], x1, y1)
        return accum

class ecc_pollardrho:
    def __init__(self, curve, order):
        self.curve = curve
        self.order = order  

    def next_x(self, xi, p, q):
        x = xi[0]
        y = xi[1]
        if x >= 0 and x <= self.curve.mod / 3: return self.curve.add(q[0], q[1], x, y)
        elif x > self.curve.mod / 3 and x <= (self.curve.mod * 2)/3: return self.curve.double(x, y)
        else: return self.curve.add(p[0], p[1], x, y)

    def next_a(self, a, xi):
        x = xi[0]
        a = mfi(a, self.order)
        if x >= 0 and x <= self.curve.mod / 3: return a
        elif x > self.curve.mod / 3 and x <= (self.curve.mod * 2) /3: return a * 2
        else: return a + 1 
  
    def next_b(self, b, xi):
        x = xi[0]
        b = mfi(b, self.order)
        if x >= 0 and x <= self.curve.mod / 3: return b+1
        elif x > self.curve.mod / 3 and x <= (self.curve.mod * 2) / 3: return b * 2
        else: return b 
    
    def solve(self, p, q):
        global itercount
        next_x = self.next_x
        next_a = self.next_a
        next_b = self.next_b
        global tests
        tests += 1
        xi = (0,1)
        a = mfi(0, self.order)
        b = mfi(0, self.order)
        x2i = (0,1)
        a2 = mfi(0, self.order)
        b2 = mfi(0, self.order)
        while True:
            old_xi = xi
            xi = next_x(xi, p, q)
            a = next_a(a, old_xi)
            b = next_b(b, old_xi)

            old_x2i = x2i
            x2i = next_x(next_x(x2i, p, q), p, q)
            a2  = next_a(next_a(a2, old_x2i), next_x(old_x2i, p, q))
            b2  = next_b(next_b(b2, old_x2i), next_x(old_x2i, p, q))
            itercount += 1
            if (xi == x2i):
                print a,b,a2,b2
                r = mfi(b-b2, self.order).modinv()
                p = mfi(a2-a, self.order)
                print 'iterations', itercount
                return p * r

if __name__ == "__main__":
    x = []
    y = []
    curve = ecc(2147483656,2060571714,2147483659)
    p = (1466500883, 1666020463)
    assert p in curve.get_point(p[0])
    bigorder = 2147440849
    pr = ecc_pollardrho(curve, bigorder)
    a = p
    for i in xrange(0,10):
      d = random.randint(0, 100)    
      print 'multiplying'
      b = curve.multiply(a[0], a[1], d)
      print 'done multiplying'
      assert curve.is_valid_point(b[0], b[1])
      print 'starting solve'
      print pr.solve(a, b), d
      print 'finishing solve'
    print itercount
    print tests
