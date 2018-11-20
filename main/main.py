from crack.crack import Crack
import math


class Main:

    def __init__(self):
        self.cracks = []
        self.N = int(input('Enter number of cracks: '))
        self.G = int(input('Enter the value of the elastic: '))
        self.N1 = 40
        self.N2 = 800
        self.init_cracks()
        self.solve()

    def init_cracks(self):
        for i in range(self.N + 1):
            self.cracks.append(Crack(i + 1))

    def solve(self):
        equations = []
        right_side = []

        # main loop
        for crack in self.cracks:
            # use extra conditions
            equations.append(self.use_extra_conditions(crack))
            right_side.append(0)

            # main formula
            collocation_points = self.get_collocation_points()
            for s_j in collocation_points:
                equations.append(self.calculate(crack, s_j))
                right_side.append(-2 * crack.tao / self.G)

        self.solve_equations(equations, right_side)

    def use_extra_conditions(self, crack):
        equation = []

        def integrand(x, n, a):
            return math.cos(n * math.acos(x / a)) / math.sqrt(1 - (x / a)**2)

        for n in range(self.N1 + 1):
            equation.append(quad(integrand, 0, 1, args=(n, crack.a)))

        return equation

    @staticmethod
    def solve_equations(a, b):
        a = np.array(a)
        b = np.array(b)

        return np.linalg.solve(a, b)

    def get_collocation_points(self):
        points = []
        for l in range(self.N1 + 1):
            points.append(math.cos(math.pi * l / (self.N1 + 1)))

        return points

    def calculate(self, crack_j, s_j):
        equation = []
        for n in range(self.N1):
            # first part of the formula
            equation.append(self.calc_U(s_j / crack_j.a, n))

            for crack_k in self.cracks:
                if self.is_current_crack(crack_j, crack_k):
                    continue
                for i in range(1, self.N2 + 1):
                    # second part of the formula
                    ksi = math.cos(math.pi * (2 * i - 1) / 2 * self.N2)
                    equation.append(crack_k.a / self.N2 * self.calc_K(crack_j, crack_k, ksi, s_j) * self.calc_T(ksi, n))

        return equation

    @staticmethod
    def is_current_crack(crack_j, crack_k):
        return crack_j.index == crack_k.index

    @staticmethod
    def calc_U(u, n):
        if n == 0:
            return 0

        return math.sin(n * math.acos(u))/math.cos(math.acos(u))

    def calc_K(self, crack_j, crack_k, ksi, s_j):
        C1 = crack_k.a * ksi * math.cos(crack_k.alpha - crack_j.alpha) - s_j + \
             (crack_k.x - crack_j.x) * math.cos(crack_j.alpha) + (crack_k.y - crack_j.y) * math.sin(crack_j.alpha)
        C2 = crack_k.a * ksi * math.sin(crack_k.alpha - crack_j.alpha) + \
             (crack_k.y - crack_j.y) * math.cos(crack_j.alpha) + (crack_j.x - crack_k.x) * math.sin(crack_j.alpha)

        return (C1 * self.G)/(C1 * C1 - C2 * C2)

    def calc_T(self, ksi, n):
        return math.cos(n * math.acos(ksi))
