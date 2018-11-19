import math


class Crack:
    def __init__(self):
        self.x = int(input('Enter coordinate x of the center crack'))
        self.y = int(input('Enter coordinate y of the center crack'))
        self.a = int(input('Enter a - length from center to end crack'))
        self.alpha = math.radians(int(input('Enter angle of inclination')))
        self.tao = int(input('Enter the value of load crack'))

