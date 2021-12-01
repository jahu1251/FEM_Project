import math


class Gauss:

    N1 = [-(1/math.sqrt(3)), (1/math.sqrt(3))]
    N2 = [-(math.sqrt(3/5)), 0, (math.sqrt(3/5))]
    N3 = [-0.861136, -0.339981, 0.339981, 0.861136]
    N4 = [-0.906180, -0.538469, 0, 0.538469, 0.906180]

    A1 = [1, 1]
    A2 = [5/9, 8/9, 5/9]
    A3 = [0.347855, 0.652145, 0.652145, 0.347855]
    A4 = [0.236927, 0.478629, 0.568889, 0.478629, 0.236927]

    Nodes_Data = [N1, N2, N3, N4]
    Factors_Data = [A1, A2, A3, A4]

    @staticmethod
    def fx(x):
        return 5*pow(x, 2) + 3*x + 6

    @staticmethod
    def fxy(x, y):
        return 5 * (x * x) * (y * y) + (3 * x * y) + 6

    def two_dimensions_gauss_integral(self, N):
        integral = 0
        j = 0
        for i in self.Factors_Data[N]:
            integral += i * self.fx(self.Nodes_Data[N][j])
            j += 1
        return integral

    def three_dimensions_gauss_integral(self, N):
        integral2 = 0
        l = 0
        for i in self.Factors_Data[N]:
            k = 0
            for j in self.Factors_Data[N]:
                integral2 += i * j * self.fxy(self.Nodes_Data[N][k], self.Nodes_Data[N][l])
                k += 1
            l += 1
        return integral2

