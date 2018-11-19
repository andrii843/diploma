from crack.crack import Crack


class Main:

    def __init__(self):
        self.cracks = []
        self.N = int(input('Enter number of cracks: '))
        self.G = int(input('Enter the value of the elastic: '))
        self.init_cracks()

    def init_cracks(self):
        for i in range(0, self.N):
            self.cracks.append(Crack())
