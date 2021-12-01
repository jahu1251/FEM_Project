import numpy as np


class Jacobian:

    auxiliary_matrix_list = []
    jacobian_matrix_list = []

    ksi_derivatives = []
    eta_derivatives = []

    rigid_values = [[0, 0.025, 0.025, 0], [0, 0, 0.025, 0.025]]

    elements = []

    def __init__(self, ksi_der, eta_der, _elements):
        self.ksi_derivatives = ksi_der
        self.eta_derivatives = eta_der

        x = len(_elements)

        self.auxiliary_matrix_list = np.zeros((x, 2, 2))
        self.jacobian_matrix_list = np.zeros((x, 2, 2))

        self.elements = _elements

    @staticmethod
    def calculate_matrix_cell_rigid(derivatives, rigid_values):
        result = 0
        for i in range(len(derivatives)):
            result = result + (derivatives[i] * rigid_values[i])

        return result

    def calculate_matrix_cell_x(self, derivatives, values):
        result_x = 0.0
        for i in range(len(derivatives)):
            result_x = result_x + (derivatives[i] * values[i].x)
        return result_x

    def calculate_matrix_cell_y(self, derivatives, values):
        result_y = 0.0
        for i in range(len(derivatives)):
            result_y = result_y + (derivatives[i] * values[i].y)
        return result_y

    def calc_jacobian_rigid(self):

        for i in range(len(self.auxiliary_matrix_list)):
            print()
            self.auxiliary_matrix_list[i][0][0] = self.calculate_matrix_cell_rigid(self.eta_derivatives[i], self.rigid_values[1])
            self.auxiliary_matrix_list[i][1][1] = self.calculate_matrix_cell_rigid(self.ksi_derivatives[i], self.rigid_values[0])

        j = 0
        for i in self.auxiliary_matrix_list:
            print(i)
            determinant = 1/np.linalg.det(i)
            self.jacobian_matrix_list[j] = i * determinant
            j = j + 1

        print("Pochodne po x i pochodne po y dla sztynego przykladu :")
        for i in self.jacobian_matrix_list:
            print(i)

        return self.jacobian_matrix_list

    def calc_jacobian(self, nodes):

        auxiliary_element_matrix_list = np.zeros((4, 2, 2))
        element_matrix_list = np.zeros((4, 2, 2))

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

        for i in jacobian_result:
            print(i)

        return jacobian_result

