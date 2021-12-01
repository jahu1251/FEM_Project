import numpy as np

class Matrix_Integration:

    def __init__(self, ksi_der, eta_der, jacobian_matrix_list, cond):

        self.ksi_derivatives = ksi_der
        self.eta_derivatives = eta_der
        self.jacobian_matrix_list = jacobian_matrix_list
        self.cond = cond

    def calculate_temp_table_dx(self):

        dx_tab = np.zeros((len(self.ksi_derivatives), 4))

        for i in range(len(self.ksi_derivatives)):
            for j in range(len(self.ksi_derivatives[0])):
                dx_tab[i][j] = (self.jacobian_matrix_list[i][j][0][0] * self.ksi_derivatives[i][j]) + (self.jacobian_matrix_list[i][j][0][1] * self.eta_derivatives[i][j])

        print("Tabela przejsciowa dx :")
        print(dx_tab)

        return dx_tab


    def calculate_temp_table_dy(self):

        dy_tab = np.zeros((len(self.eta_derivatives), 4))

        for i in range(len(self.ksi_derivatives)):
            for j in range(len(self.ksi_derivatives[0])):
                dy_tab[i][j] = (self.jacobian_matrix_list[i][j][1][0] * self.ksi_derivatives[i][j]) + (self.jacobian_matrix_list[i][j][1][1] * self.eta_derivatives[i][j])

        print("Tabela przejsciowa dy :")
        print(dy_tab)

        return dy_tab

    def calculate_matrix_for_element(self, dx_tab, dy_tab, jacobian):

        y = 1/np.linalg.det(jacobian)

        H_pc1 = np.zeros((4, 4))
        H_pc2 = np.zeros((4, 4))
        H_pc3 = np.zeros((4, 4))
        H_pc4 = np.zeros((4, 4))
        H_pc_list = [H_pc1, H_pc2, H_pc3, H_pc4]

        for i in range(4):
            temp_matrix_x = np.matrix(dx_tab[i])
            temp_matrix_y = np.matrix(dy_tab[i])
            temp_matrix_x_transposed = temp_matrix_x.transpose()
            temp_matrix_y_transposed = temp_matrix_y.transpose()
            H_pc_list[i] = self.cond*(np.outer(temp_matrix_x, temp_matrix_x_transposed) + np.outer(temp_matrix_y, temp_matrix_y_transposed))*y

        z = 1
        for x in H_pc_list:
            print("[H]Pc", z)
            z = z + 1
            print(np.round(x, 2))

        print("Macierz H: ")
        H = H_pc_list[0] + H_pc_list[1] + H_pc_list[2] + H_pc_list[3]
        print(np.round(H, 2))

        return H

    def calculate_H_matrixes_for_grid(self):

        global_h_matrix_list = np.zeros((len(self.jacobian_matrix_list), 4, 4))

        for i in range(len(self.jacobian_matrix_list)):

            print("=========================== element", i, "=================================")
            temp_dx_table = self.calculate_temp_table_dx()
            temp_dy_table = self.calculate_temp_table_dy()
            global_h_matrix_list[i] = self.calculate_matrix_for_element(temp_dx_table, temp_dy_table, self.jacobian_matrix_list[i])
            print(global_h_matrix_list[i])

        self.global_h_matrix_list = global_h_matrix_list
