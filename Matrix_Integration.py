import numpy as np


class Matrix_Integration:

    pw3 = np.array([5 / 9, 8 / 9, 5 / 9, 5 / 9, 8 / 9, 5 / 9, 5 / 9, 8 / 9, 5 / 9, 5 / 9, 8 / 9, 5 / 9])
    pwm = np.array([pw3[0] * pw3[0], pw3[0] * pw3[1], pw3[1] * pw3[1]])
    pwp = np.array([pwm[0], pwm[1], pwm[0], pwm[1], pwm[2], pwm[1], pwm[0], pwm[1], pwm[0]])

    def __init__(self, ksi_der, eta_der, jacobian_matrix_list, cond, c, ro, N_values):

        self.ksi_derivatives = np.array(ksi_der)
        self.eta_derivatives = np.array(eta_der)
        self.jacobian_matrix_list = np.array(jacobian_matrix_list)
        self.cond = cond
        self.ro = ro
        self.c = c
        self.N_values = np.array(N_values)

    def calculate_temp_table_dx(self):

        dx_tab = np.zeros((len(self.ksi_derivatives), 4))

        for i in range(len(self.ksi_derivatives)):
            for j in range(len(self.ksi_derivatives[0])):
                dx_tab[i][j] = (self.jacobian_matrix_list[i][j][0][0] * self.ksi_derivatives[i][j]) + (self.jacobian_matrix_list[i][j][0][1] * self.eta_derivatives[i][j])

        # print("Tabela przejsciowa dx :")
        # print(dx_tab)

        return dx_tab

    def calculate_temp_table_dy(self):

        dy_tab = np.zeros((len(self.eta_derivatives), 4))

        for i in range(len(self.ksi_derivatives)):
            for j in range(len(self.ksi_derivatives[0])):
                dy_tab[i][j] = (self.jacobian_matrix_list[i][j][1][0] * self.ksi_derivatives[i][j]) + (self.jacobian_matrix_list[i][j][1][1] * self.eta_derivatives[i][j])

        # print("Tabela przejsciowa dy :")
        # print(dy_tab)

        return dy_tab

    def calculate_matrix_h_for_element(self, dx_tab, dy_tab, jacobian):

        y = 1/np.linalg.det(jacobian)

        h_pc_list = np.zeros((4, 4, 4))

        for i in range(4):
            temp_matrix_x = np.matrix(dx_tab[i])
            temp_matrix_y = np.matrix(dy_tab[i])
            temp_matrix_x_transposed = temp_matrix_x.transpose()
            temp_matrix_y_transposed = temp_matrix_y.transpose()
            h_pc_list[i] = self.cond*(np.outer(temp_matrix_x, temp_matrix_x_transposed) + np.outer(temp_matrix_y, temp_matrix_y_transposed))*y[i]

        # z = 1
        # for x in H_pc_list:
        #     print("[H]Pc", z)
        #     z = z + 1
        #     print(np.round(x, 2))

        print("Macierz H: ")
        h = h_pc_list[0] + h_pc_list[1] + h_pc_list[2] + h_pc_list[3]
        print(np.round(h, 2))

        return h

    def calculate_matrix_h_for_element_3(self, dx_tab, dy_tab, jacobian):

        y = 1/np.linalg.det(jacobian)

        h_pc_list = np.zeros((9, 4, 4))

        for i in range(9):
            temp_matrix_x = np.matrix(dx_tab[i])
            temp_matrix_y = np.matrix(dy_tab[i])
            temp_matrix_x_transposed = temp_matrix_x.transpose()
            temp_matrix_y_transposed = temp_matrix_y.transpose()
            h_pc_list[i] = self.pwp[i]*(self.cond*(np.outer(temp_matrix_x, temp_matrix_x_transposed) + np.outer(temp_matrix_y, temp_matrix_y_transposed))*y[i])

        # z = 1
        # for x in H_pc_list:
        #     print("[H]Pc", z)
        #     z = z + 1
        #     print(np.round(x, 2))

        print("Macierz H: ")
        h = h_pc_list[0] + h_pc_list[1] + h_pc_list[2] + h_pc_list[3]+ h_pc_list[4] + h_pc_list[5] + h_pc_list[6]+ h_pc_list[7] + h_pc_list[8]
        print(np.round(h, 2))

        return h

    def calculate_matrix_c_for_element(self, jacobian):

        det_j = 1/np.linalg.det(jacobian)

        c_pc_list = np.zeros((4, 4, 4))

        for i in range(4):
            temp_N_Values = np.matrix(self.N_values[i])
            c_pc_list[i] = self.c * self.ro * (np.outer(temp_N_Values, temp_N_Values.transpose())) * det_j

        # z = 1
        # for x in C_pc_list:
        #     print("[C]Pc", z)
        #     z = z + 1
        #     print(np.round(x, 2))

        print("Macierz C: ")
        c = c_pc_list[0] + c_pc_list[1] + c_pc_list[2] + c_pc_list[3]
        print(np.round(c, 2))

        return c

    def calculate_matrix_c_for_element_3(self, jacobian):

        det_j = 1/np.linalg.det(jacobian)

        c_pc_list = np.zeros((9, 4, 4))

        for i in range(9):
            temp_N_Values = np.matrix(self.N_values[i])
            c_pc_list[i] = self.pwp[i]*(self.c * self.ro * (np.outer(temp_N_Values, temp_N_Values.transpose())) * det_j[i])

        # z = 1
        # for x in C_pc_list:
        #     print("[C]Pc", z)
        #     z = z + 1
        #     print(np.round(x, 2))

        print("Macierz C: ")
        c = c_pc_list[0] + c_pc_list[1] + c_pc_list[2] + c_pc_list[3] + c_pc_list[4] + c_pc_list[5] + c_pc_list[6]+ c_pc_list[7]+ c_pc_list[8]
        print(np.round(c, 2))

        return c

    def calculate_H_matrixes_for_grid(self):

        global_h_matrix_list = np.zeros((len(self.jacobian_matrix_list), 4, 4))
        global_c_matrix_list = np.zeros((len(self.jacobian_matrix_list), 4, 4))

        for i in range(len(self.jacobian_matrix_list)):

            print("=========================== element", i, "=================================")
            temp_dx_table = self.calculate_temp_table_dx()
            temp_dy_table = self.calculate_temp_table_dy()
            global_h_matrix_list[i] = self.calculate_matrix_h_for_element(temp_dx_table, temp_dy_table, self.jacobian_matrix_list[i])
            global_c_matrix_list[i] = self.calculate_matrix_c_for_element(self.jacobian_matrix_list[i])
            print(global_h_matrix_list[i])
            print(global_c_matrix_list[i])

        self.global_h_matrix_list = global_h_matrix_list
        self.global_c_matrix_list = global_c_matrix_list

    def calculate_H_matrixes_for_grid_3(self):

        global_h_matrix_list = np.zeros((len(self.jacobian_matrix_list), 4, 4))
        global_c_matrix_list = np.zeros((len(self.jacobian_matrix_list), 4, 4))

        for i in range(len(self.jacobian_matrix_list)):

            print("=========================== element", i, "=================================")
            temp_dx_table = self.calculate_temp_table_dx()
            temp_dy_table = self.calculate_temp_table_dy()
            global_h_matrix_list[i] = self.calculate_matrix_h_for_element_3(temp_dx_table, temp_dy_table, self.jacobian_matrix_list[i])
            global_c_matrix_list[i] = self.calculate_matrix_c_for_element_3(self.jacobian_matrix_list[i])
            print(global_h_matrix_list[i])
            print(global_c_matrix_list[i])

        self.global_h_matrix_list = global_h_matrix_list
        self.global_c_matrix_list = global_c_matrix_list
