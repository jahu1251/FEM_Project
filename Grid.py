from Element import Element
from Node import Node
import numpy as np


class Grid:
    H, B = None, None  # wysokość i szerokość siatki
    nH, nB, nN, nE = None, None, None, None  # liczba elementów w pionie, poziomie, liczba węzłów oraz liczba elementów
    stepX, stepY = None, None  # krok w pionie i poziomie
    elements = []  # lista elementów
    temp_x = 0  # pomocnicze x
    temp_y = 0  # pomocnicze y
    temp_id = 1  # pomocnicze id
    weights_3 = np.array([5 / 9, 8 / 9, 5 / 9])  # wagi dla 3 punktowego schematu całkowania

    def __init__(self, nH, nB, H, B, temp, time_step, init_temp):
        self.nH = nH
        self.nB = nB
        self.H = H
        self.B = B
        self.nN = self.nH * self.nB  # łączna liczba węzłów = liczba w pionie x liczba w poziomie
        self.nodes_list = [Node] * self.nN  # lista pomocnicza do przechoeywania node'ów
        self.global_H_matrix = np.zeros((self.nN, self.nN))  # globalna macierz H
        self.global_C_matrix = np.zeros((self.nN, self.nN))  # globalna macierz C
        self.global_P_matrix = np.zeros(self.nN)  # globalna macierz P
        self.temp = temp  # ciepło właściwe
        self.time_step = time_step  # krok czasowy
        self.init_temp = init_temp  # temperatura początkowa węzłów
        self.calculate_parameters()
        self.generate_elements()

    def calculate_parameters(self):

        self.nE = (self.nH - 1) * (
                self.nB - 1)  # liczba elementów = (liczba węzłów w pionie -1) x (liczba węzłów w poziomie -1)
        self.stepX = self.B / (self.nB - 1)  # krok w poziomie = szerokość / liczba elementów
        self.stepY = self.H / (self.nH - 1)  # krok w pionie = wysokość / liczba elementów w pionie

        print("Liczba węzłów : ", self.nN, ", liczba elementów : ", self.nE, ", długość kroku X : ", self.stepX,
              ", długość kroku Y : ", self.stepY)

    def generate_elements(self):

        k = 0
        for i in range(self.nB - 1):  # iteracja po wierszu
            for j in range(self.nH - 1):  # iteracja po kolumnie
                bc0, bc1, bc2, bc3 = False, False, False, False  # pomocnicze zmienne do określania warunku brzegowego

                # sprawdzanie warunku brzegowego
                if self.temp_x == 0:
                    bc0 = True
                    bc3 = True

                if self.temp_x + self.stepX == self.B:
                    bc1 = True
                    bc2 = True

                if self.temp_y == 0:
                    bc0 = True
                    bc1 = True

                if self.temp_y + self.stepY == self.H:
                    bc2 = True
                    bc3 = True

                # definicja nowych node'ów
                node_1 = Node(self.temp_x, self.temp_y, bc0, self.init_temp)
                node_2 = Node(self.temp_x + self.stepX, self.temp_y, bc1, self.init_temp)
                node_3 = Node(self.temp_x + self.stepX, self.temp_y + self.stepY, bc2, self.init_temp)
                node_4 = Node(self.temp_x, self.temp_y + self.stepY, bc3, self.init_temp)

                # tworzenie nowego elementu z odpowiednimi zmiennymi
                el_obj = Element(self.temp_id,
                                 self.temp_id + self.nH,
                                 self.temp_id + self.nH + 1,
                                 self.temp_id + 1,
                                 node_1,
                                 node_2,
                                 node_3,
                                 node_4)

                # dodanie nowego elementu do siatki, inkrementacja zmiennych
                self.elements.append(el_obj)
                self.temp_y = self.temp_y + self.stepY
                self.temp_id = self.temp_id + 1

                k = k + 1

            # zerowanie y (nowa kolumna), inkrementacja x oraz zwiększenie id
            self.temp_y = 0
            self.temp_x = self.temp_x + self.stepX
            self.temp_id = self.temp_id + 1

        # przypisanie node'ów do listy pomocnieczej
        for i in self.elements:
            self.nodes_list[i.ID[0] - 1] = i.nodes[0]
            self.nodes_list[i.ID[1] - 1] = i.nodes[1]
            self.nodes_list[i.ID[2] - 1] = i.nodes[2]
            self.nodes_list[i.ID[3] - 1] = i.nodes[3]

    # wypisanie siatki na ekran
    def print_grid(self):
        temp = 0
        for i in self.elements:
            temp += 1
            print(temp, ": ",
                  i.ID[0], i.ID[1], i.ID[2], i.ID[3], ": (",
                  round(i.nodes[0].x, 3), ",", round(i.nodes[0].y, 3), i.nodes[0].is_bc, "), (",
                  round(i.nodes[1].x, 3), ",", round(i.nodes[1].y, 3), i.nodes[1].is_bc, "), (",
                  round(i.nodes[2].x, 3), ",", round(i.nodes[2].y, 3), i.nodes[2].is_bc, "), (",
                  round(i.nodes[3].x, 3), ",", round(i.nodes[3].y, 3), i.nodes[3].is_bc, ")"
                  )

    # przpisanie macierzy h dla każdego elementu (liczone są one w innnej klasie)
    def add_h_matrixes_to_elements(self, h_matrix):

        for i in range(len(self.elements)):
            self.elements[i].h_matrix = h_matrix[i]

    # przypisanie macierz C do elementów (liczone w innej klasie)
    def add_c_matrixes_to_elements(self, c_matrix):

        for i in range(len(self.elements)):
            self.elements[i].c_matrix = c_matrix[i]

    # obliczanie macierzy hbc dla elementów
    def calculate_hbc_matrixes_for_elements(self, hbc_matrix, cond):

        j = 0
        for i in self.elements:

            temp_hbc_matrix = np.zeros((4, 4, 4))
            temp_vect_p = np.zeros((4, 4))

            # sprawdzenie warunków brzegowych węzłów w celu obliczenia odpowiednich macierzy hbc oraz wektorow p
            if i.nodes[0].is_bc is True and i.nodes[1].is_bc is True:
                temp_hbc_matrix[0] = cond * (
                        np.outer(hbc_matrix[0, 0], hbc_matrix[0, 0].transpose()) +
                        np.outer(hbc_matrix[0, 1], hbc_matrix[0, 1].transpose())) * (
                                             i.nodes[1].x - i.nodes[0].x) / 2

                temp_vect_p[0] = cond * (
                        hbc_matrix[0, 0] * self.temp + hbc_matrix[0, 1] * self.temp) * (i.nodes[1].x - i.nodes[0].x) / 2

            if i.nodes[1].is_bc is True and i.nodes[2].is_bc is True:
                temp_hbc_matrix[1] = cond * (
                        np.outer(hbc_matrix[1, 0], hbc_matrix[1, 0].transpose()) +
                        np.outer(hbc_matrix[1, 1], hbc_matrix[1, 1].transpose())) * (
                                             i.nodes[2].y - i.nodes[1].y) / 2

                temp_vect_p[1] = cond * (hbc_matrix[1, 0] * self.temp + hbc_matrix[1, 1] * self.temp) * (i.nodes[2].y - i.nodes[1].y) / 2

            if i.nodes[2].is_bc is True and i.nodes[3].is_bc is True:
                temp_hbc_matrix[2] = cond * (
                        np.outer(hbc_matrix[2, 0], hbc_matrix[2, 0].transpose()) +
                        np.outer(hbc_matrix[2, 1], hbc_matrix[2, 1].transpose())) * (
                                             i.nodes[2].x - i.nodes[3].x) / 2

                temp_vect_p[2] = cond * (hbc_matrix[2, 0] * self.temp + hbc_matrix[2, 1] * self.temp) * (i.nodes[2].x - i.nodes[3].x) / 2

            if i.nodes[3].is_bc is True and i.nodes[0].is_bc is True:
                temp_hbc_matrix[3] = cond * (
                        np.outer(hbc_matrix[3, 0], hbc_matrix[3, 0].transpose()) +
                        np.outer(hbc_matrix[3, 1], hbc_matrix[3, 1].transpose())) * (
                                             i.nodes[3].y - i.nodes[0].y) / 2

                temp_vect_p[3] = cond * (hbc_matrix[3, 0] * self.temp + hbc_matrix[3, 1] * self.temp) * (i.nodes[3].y - i.nodes[0].y) / 2

            # sumowanie macierzy hbc oraz p
            i.hbc_matrix = temp_hbc_matrix[0] + temp_hbc_matrix[1] + temp_hbc_matrix[2] + temp_hbc_matrix[3]
            i.p_vector = temp_vect_p[0] + temp_vect_p[1] + temp_vect_p[2] + temp_vect_p[3]

            print(i.p_vector)
            j = j + 1

        # k = 0
        # for i in self.elements:
        #     print("element ", k)
        #     print(i.hbc_matrix)
        #     k = k + 1

    # obliczanie macierzy hbc dla elementów w 3 punktowym schemacie całkowania
    def calculate_hbc_matrixes_for_elements_3(self, hbc_matrix, cond):

        j = 0
        for i in self.elements:
            temp_hbc_matrix = np.zeros((4, 4, 4))
            temp_vect_p = np.zeros((4, 4))

            # sprawdzenie warunków brzegowych węzłów w celu obliczenia odpowiednich macierzy hbc oraz wektorow p
            if i.nodes[0].is_bc is True and i.nodes[1].is_bc is True:
                temp_hbc_matrix[0] = cond * (
                        self.weights_3[0] * np.outer(hbc_matrix[0, 0], hbc_matrix[0, 0].transpose()) +
                        self.weights_3[1] * np.outer(hbc_matrix[0, 1], hbc_matrix[0, 1].transpose()) +
                        self.weights_3[2] * np.outer(hbc_matrix[0, 2], hbc_matrix[0, 2].transpose())) * (
                                             i.nodes[1].x - i.nodes[0].x) / 2
                temp_vect_p[0] = cond * (
                        self.weights_3[0] * hbc_matrix[0, 0] * self.temp +
                        self.weights_3[1] * hbc_matrix[0, 1] * self.temp +
                        self.weights_3[2] * hbc_matrix[0, 2] * self.temp) * (
                                         i.nodes[1].x - i.nodes[0].x) / 2

            if i.nodes[1].is_bc is True and i.nodes[2].is_bc is True:
                temp_hbc_matrix[1] = cond * (
                        self.weights_3[0] * np.outer(hbc_matrix[1, 0], hbc_matrix[1, 0].transpose()) +
                        self.weights_3[1] * np.outer(hbc_matrix[1, 1], hbc_matrix[1, 1].transpose()) +
                        self.weights_3[2] * np.outer(hbc_matrix[1, 2], hbc_matrix[1, 2].transpose())) * (
                                             i.nodes[2].y - i.nodes[1].y) / 2

                temp_vect_p[1] = cond * (
                        self.weights_3[0] * hbc_matrix[1, 0] * self.temp +
                        self.weights_3[1] * hbc_matrix[1, 1] * self.temp +
                        self.weights_3[2] * hbc_matrix[1, 2] * self.temp) * (
                                         i.nodes[2].y - i.nodes[1].y) / 2

            if i.nodes[2].is_bc is True and i.nodes[3].is_bc is True:
                temp_hbc_matrix[2] = cond * (
                        self.weights_3[0] * np.outer(hbc_matrix[2, 0], hbc_matrix[2, 0].transpose()) +
                        self.weights_3[1] * np.outer(hbc_matrix[2, 1], hbc_matrix[2, 1].transpose()) +
                        self.weights_3[2] * np.outer(hbc_matrix[2, 2], hbc_matrix[2, 2].transpose())) * (
                                             i.nodes[2].x - i.nodes[3].x) / 2

                temp_vect_p[2] = cond * (
                        self.weights_3[0] * hbc_matrix[2, 0] * self.temp +
                        self.weights_3[1] * hbc_matrix[2, 1] * self.temp +
                        self.weights_3[2] * hbc_matrix[2, 2] * self.temp) * (
                                         i.nodes[2].x - i.nodes[3].x) / 2

            if i.nodes[3].is_bc is True and i.nodes[0].is_bc is True:
                temp_hbc_matrix[3] = cond * (
                        self.weights_3[0] * np.outer(hbc_matrix[3, 0], hbc_matrix[3, 0].transpose()) +
                        self.weights_3[1] * np.outer(hbc_matrix[3, 1], hbc_matrix[3, 1].transpose()) +
                        self.weights_3[2] * np.outer(hbc_matrix[3, 2], hbc_matrix[3, 2].transpose())) * (
                                             i.nodes[3].y - i.nodes[0].y) / 2

                temp_vect_p[3] = cond * (
                        self.weights_3[0] * hbc_matrix[3, 0] * self.temp +
                        self.weights_3[1] * hbc_matrix[3, 1] * self.temp +
                        self.weights_3[2] * hbc_matrix[3, 2] * self.temp) * (
                                         i.nodes[3].y - i.nodes[0].y) / 2

            i.hbc_matrix = temp_hbc_matrix[0] + temp_hbc_matrix[1] + temp_hbc_matrix[2] + temp_hbc_matrix[3]
            i.p_vector = temp_vect_p[0] + temp_vect_p[1] + temp_vect_p[2] + temp_vect_p[3]

            j = j + 1

        k = 0
        for i in self.elements:
            print("element ", k, "macierz HBC : ")
            print(i.hbc_matrix)
            k = k + 1

            print("element ", k, "wektor P : ")
            print(i.p_vector)
            k = k + 1

    # agregacja macierzy H oraz wektora P
    def aggregation(self):

        for i in range(self.nE):
            for j in range(4):
                self.global_P_matrix[self.elements[i].ID[j] - 1] = self.global_P_matrix[self.elements[i].ID[j] - 1] + self.elements[i].p_vector[j]
                for k in range(4):
                    self.global_H_matrix[self.elements[i].ID[j] - 1, self.elements[i].ID[k] - 1] = \
                        self.global_H_matrix[self.elements[i].ID[j] - 1, self.elements[i].ID[k] - 1] + \
                        self.elements[i].h_matrix[j, k] + \
                        self.elements[i].hbc_matrix[j, k]

                    self.global_C_matrix[self.elements[i].ID[j] - 1, self.elements[i].ID[k] - 1] = \
                        self.global_C_matrix[self.elements[i].ID[j] - 1, self.elements[i].ID[k] - 1] + \
                        self.elements[i].c_matrix[j, k]

        for i in range(self.nN):
            for j in range(self.nN):
                self.global_H_matrix[i, j] = self.global_H_matrix[i, j] + self.global_C_matrix[i, j] / self.time_step

    # obliczanie globalnej macierzy H oraz wektora P ze wzoru
    def calculate_p_and_h_matrixes(self):

        vector = np.zeros(self.nN)

        for i in range(len(self.global_P_matrix)):
            vector[i] = self.global_P_matrix[i]

        for i in range(self.nN):
            for j in range(self.nN):
                vector[i] = vector[i] + (self.global_C_matrix[i, j] / self.time_step) * self.nodes_list[j].init_temp

        temp_vector = np.linalg.solve(self.global_H_matrix, vector)

        for o in range(len(temp_vector)):
            self.nodes_list[o].init_temp = temp_vector[o]

        print("Maksymalna temperatura : ", np.amax(temp_vector), "Minimalna tempareatura : ", np.amin(temp_vector))

    @staticmethod
    def equation_solve(a, b):

        num = len(a)
        n = np.empty(num)
        x1 = np.empty(num)
        x2 = np.empty(num)
        M = np.empty((num, num))

        i = 0
        while i < num:
            n[i] = 1 / a[i][i]
            i += 1

        # Calculate M = -D^-1 (L + U)
        i = 0
        while i < num:
            j = 0
            while j < num:
                if i == j:
                    M[i][j] = 0
                else:
                    M[i][j] = -(a[i][j] * n[i])
                j += 1
            i += 1

        for k in range(0, 100):
            i = 0
            while i < num:
                x2[i] = n[i] * b[i]
                j = 0
                while j < num:
                    x2[i] += M[i][j] * x1[j]
                    j += 1
                i += 1
            i = 0
            while i < num:
                x1[i] = x2[i]
                i += 1
        return list(x1)
