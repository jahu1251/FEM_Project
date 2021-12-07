import math
import numpy as np


class Element4_2D:

    def N1(x):
        return (-1) * 0.25 * (1.0 - x)

    def N2(x):
        return 0.25 * (1.0 - x)

    def N3(x):
        return 0.25 * (1.0 + x)

    def N4(x):
        return (-1) * 0.25 * (1.0 + x)

    def N1_hbc(x, y):
        return 0.25 * (1. - x) * (1. - y)

    def N2_hbc(x, y):
        return 0.25 * (1. + x) * (1. - y)

    def N3_hbc(x, y):
        return 0.25 * (1. + x) * (1. + y)

    def N4_hbc(x, y):
        return 0.25 * (1. - x) * (1. + y)

    #devariatives things
    Ksi_2 = [-(1 / math.sqrt(3)), (1 / math.sqrt(3)), (1 / math.sqrt(3)), -(1 / math.sqrt(3))]
    Eta_2 = [-(1 / math.sqrt(3)), -(1 / math.sqrt(3)), (1 / math.sqrt(3)), (1 / math.sqrt(3))]
    Ksi_3 = [-(math.sqrt(3/5)), 0, (math.sqrt(3/5)), -(math.sqrt(3/5)), 0, math.sqrt(3/5), -(math.sqrt(3/5)), 0, math.sqrt(3/5)]
    Eta_3 = [-(math.sqrt(3/5)), -(math.sqrt(3/5)), -(math.sqrt(3/5)), 0, 0, 0, math.sqrt(3/5), math.sqrt(3/5), math.sqrt(3/5)]
    functions = [N1, N2, N3, N4]
    functions2 = [N1, N4, N3, N2]
    derivativeKsi = []
    derivativeEta = []

    #hbc things
    Ksi_2_wall_points_1 = [-(1 / math.sqrt(3)), (1 / math.sqrt(3))]
    Ksi_2_wall_points_2 = [1, 1]
    Ksi_2_wall_points_3 = [(1 / math.sqrt(3)), -(1 / math.sqrt(3))]
    Ksi_2_wall_points_4 = [-1, -1]

    Eta_2_wall_points_1 = [-1, -1]
    Eta_2_wall_points_2 = [-(1 / math.sqrt(3)), (1 / math.sqrt(3))]
    Eta_2_wall_points_3 = [1, 1]
    Eta_2_wall_points_4 = [(1 / math.sqrt(3)), -(1 / math.sqrt(3))]

    Ksi_3_wall_points_1 = [-(math.sqrt(3 / 5)), 0, (math.sqrt(3 / 5))]
    Ksi_3_wall_points_2 = [1, 1, 1]
    Ksi_3_wall_points_3 = [(math.sqrt(3 / 5)), 0, -(math.sqrt(3 / 5))]
    Ksi_3_wall_points_4 = [-1, -1, -1]

    Eta_3_wall_points_1 = [-1, -1, -1]
    Eta_3_wall_points_2 = [-(math.sqrt(3 / 5)), 0, (math.sqrt(3 / 5))]
    Eta_3_wall_points_3 = [1, 1, 1]
    Eta_3_wall_points_4 = [(math.sqrt(3 / 5)), 0, -(math.sqrt(3 / 5))]

    functions_hbc = [N1_hbc, N2_hbc, N3_hbc, N4_hbc]
    point_weights_3 = [5/9, 8/9, 5/9, 5/9, 8/9, 5/9, 5/9, 8/9, 5/9, 5/9, 8/9, 5/9]
    walls_hbc = np.zeros((4, 4, 4))

    Ksi_wall_points_2_list = [Ksi_2_wall_points_1, Ksi_2_wall_points_2, Ksi_2_wall_points_3, Ksi_2_wall_points_4]
    Eta_wall_points_2_list = [Eta_2_wall_points_1, Eta_2_wall_points_2, Eta_2_wall_points_3, Eta_2_wall_points_4]

    Ksi_wall_points_3_list = [Ksi_3_wall_points_1, Ksi_3_wall_points_2, Ksi_3_wall_points_3, Ksi_3_wall_points_4]
    Eta_wall_points_3_list = [Eta_3_wall_points_1, Eta_3_wall_points_2, Eta_3_wall_points_3, Eta_3_wall_points_4]

    p_vect = np.zeros((4,4))

    N_values = np.zeros((4, 4))
    N_values_3 = np.zeros((4, 9))

    def __init__(self, cond, integ_points, c, ro):
        self.cond = cond
        self.ro = ro
        self.c = c

        if integ_points == 2:
            self.temp_matrix_for_hbc = np.zeros((4, 2, 4))
        if integ_points == 3:
            self.temp_matrix_for_hbc = np.zeros((4, 3, 4))

    def two_dim_2(self, tab, n):
        final = []
        for i in tab:
            final.append(n(i))
        return final

    def calculate_2(self):
        for i in range(4):
            for j in range(4):
                self.N_values[i, j] = self.functions_hbc[i](self.Ksi_2[j], self.Eta_2[j])

        output = []
        for i in self.functions:
            output.append(self.two_dim_2(self.Eta_2, i))

        self.derivativeKsi = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]

        for z in range(4):
            for x in range(4):
                self.derivativeKsi[z][x] = output[x][z]

        output2 = []
        for i in self.functions2:
            output2.append(self.two_dim_2(self.Ksi_2, i))

        # print("Ksi")
        # for i in self.derivativeKsi:
        #     print(i)

        self.derivativeEta = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]

        for z in range(4):
            for x in range(4):
                self.derivativeEta[z][x] = output2[x][z]

        # print("Eta")
        # for i in self.derivativeEta:
        #     print(i)

    def calculate_3(self):

        for i in range(4):
            for j in range(9):
                self.N_values_3[i, j] = self.functions_hbc[i](self.Ksi_3[j], self.Eta_3[j])

        output = []
        for i in self.functions:
            output.append(self.two_dim_2(self.Eta_3, i))

        self.derivativeKsi = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0],
                  [0, 0, 0, 0], [0, 0, 0, 0]]

        for z in range(9):
            for x in range(4):
                self.derivativeKsi[z][x] = output[x][z]

        # print("Ksi")
        # for i in self.derivativeKsi:
        #     print(i)

        output2 = []
        for i in self.functions2:
            output2.append(self.two_dim_2(self.Ksi_3, i))

        self.derivativeEta = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0],
                              [0, 0, 0, 0],
                              [0, 0, 0, 0], [0, 0, 0, 0]]

        for z in range(9):
            for x in range(4):
                self.derivativeEta[z][x] = output2[x][z]

        # print("Eta")
        # for i in self.derivativeEta:
        #     print(i)

    def calculate_temp_matrix_for_hbc_2(self):
        for i in range(4):
            for j in range(4):
                for k in range(2):
                    self.temp_matrix_for_hbc[i, k, j] = self.functions_hbc[j](self.Ksi_wall_points_2_list[i][k], self.Eta_wall_points_2_list[i][k])
            # print("element ", i)
            # print(self.temp_matrix_for_hbc[i])


    def calculate_temp_matrix_for_hbc_3(self):
        for i in range(4):
            for j in range(4):
                for k in range(3):
                    self.temp_matrix_for_hbc[i, k, j] = self.functions_hbc[j](self.Ksi_wall_points_3_list[i][k], self.Eta_wall_points_3_list[i][k])
            # print("element ", i)
            # print(self.temp_matrix_for_hbc[i])
