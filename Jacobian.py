import numpy as np


class Jacobian:

    def __init__(self, ksi_der, eta_der, _elements):
        self.ksi_derivatives = np.array(ksi_der)
        self.eta_derivatives = np.array(eta_der)

        x = len(_elements)

        self.auxiliary_matrix_list = np.zeros((x, 2, 2))
        self.jacobian_matrix_list = np.zeros((x, 2, 2))

        self.elements = _elements

    @staticmethod
    def calculate_matrix_cell_x(derivatives, values):
        result_x = 0.0
        for i in range(len(derivatives)):
            result_x = result_x + (derivatives[i] * values[i].x)
        return result_x

    @staticmethod
    def calculate_matrix_cell_y(derivatives, values):
        result_y = 0.0
        for i in range(len(derivatives)):
            result_y = result_y + (derivatives[i] * values[i].y)
        return result_y

    def calc_jacobian(self, nodes):

        auxiliary_element_matrix_list = np.zeros((len(self.ksi_derivatives), 2, 2))
        element_matrix_list = np.zeros((len(self.ksi_derivatives), 2, 2))

        for i in range(len(auxiliary_element_matrix_list)):
            auxiliary_element_matrix_list[i][1][1] = self.calculate_matrix_cell_x(self.ksi_derivatives[i], nodes)
            auxiliary_element_matrix_list[i][0][0] = self.calculate_matrix_cell_y(self.eta_derivatives[i], nodes)

        j = 0
        for i in auxiliary_element_matrix_list:
            determinant = 1/np.linalg.det(i)
            element_matrix_list[j] = i * determinant
            j = j + 1

        return element_matrix_list

    def calculate_grid_jacobians(self):
        jacobian_result = []
        j = 0
        for i in self.elements:
            jacobian_result.append(self.calc_jacobian(i.nodes))

        # for i in jacobian_result:
        #     print(i)

        return jacobian_result

