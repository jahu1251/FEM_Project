from Element import Element
from Node import Node
import numpy as np


class Grid:
    H, B = None, None   # wysokość i szerokość siatki
    nH, nB, nN, nE = None, None, None, None     # liczba elementów w pionie, poziomie, liczba węzłów oraz liczba elementów
    stepX, stepY = None, None       # krok w pionie i poziomie
    elements = []       #lista elementów
    temp_x = 0      # pomocnicze x
    temp_y = 0      # pomocnicze y
    temp_id = 1     # pomocnicze id
    weights_3 = [5 / 9, 8 / 9, 5 / 9]


    def __init__(self, nH, nB, H, B, temp, time_step, init_temp):
        self.nH = nH
        self.nB = nB
        self.H = H
        self.B = B
        self.nN = self.nH * self.nB  # łączna liczba węzłów = liczba w pionie x liczba w poziomie
        self.nodes_list = np.zeros((self.nN))
        self.global_H_matrix = np.zeros((self.nN, self.nN))
        self.global_C_matrix = np.zeros((self.nN, self.nN))
        self.global_P_matrix = np.zeros((self.nN))
        self.temp = temp
        self.time_step = time_step
        self.init_temp = init_temp
        self.calculate_parameters()
        self.generate_elements()


    def calculate_parameters(self):
        self.nE = (self.nH - 1) * (self.nB - 1)     # liczba elementów = (liczba węzłów w pionie -1) x
                                                            #                   (liczba węzłów w poziomie -1)
        self.stepX = self.B / (self.nB - 1)     # krok w poziomie = szerokość / liczba elementów
        self.stepY = self.H / (self.nH - 1)     # krok w pionie = wysokość / liczba elementów w pionie

        print("Liczba węzłów : ", self.nN, ", liczba elementów : ", self.nE, ", długość kroku X : ", self.stepX,
              ", długość kroku Y : ", self.stepY)

    def generate_elements(self):

        for i in range(self.nB-1):      # iteracja po wierszu
            for j in range(self.nH-1):      # iteracja po kolumnie
                bc0, bc1, bc2, bc3 = False, False, False, False     # pomocnicze zmienne do określania warunku brzegowego

                #sprawdzanie warunku brzegowego
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

                node_1 = Node(self.temp_x, self.temp_y, bc0, self.init_temp)
                node_2 = Node(self.temp_x + self.stepX, self.temp_y, bc1, self.init_temp)
                node_3 = Node(self.temp_x + self.stepX, self.temp_y + self.stepY, bc2, self.init_temp)
                node_4 = Node(self.temp_x, self.temp_y + self.stepY, bc3, self.init_temp)

                self.nodes_list[i] = node_1
                self.nodes_list[i + self.nH + 1] = node_2
                self.nodes_list[i + self.nH + 2] = node_3
                self.nodes_list[i + 1] = node_4


                #tworzenie nowego elementu z odpowiednimi zmiennymi
                el_obj = Element(self.temp_id,
                                 self.temp_id + self.nH,
                                 self.temp_id + self.nH + 1,
                                 self.temp_id + 1,
                                 node_1,
                                 node_2,
                                 node_3,
                                 node_4)

                #dodanie nowego elementu do siatki, inkrementacja zmiennych
                self.elements.append(el_obj)
                self.temp_y = self.temp_y + self.stepY
                self.temp_id = self.temp_id + 1

            #zerowanie y (nowa kolumna), inkrementacja x oraz zwiększenie id
            self.temp_y = 0
            self.temp_x = self.temp_x + self.stepX
            self.temp_id = self.temp_id + 1

    #wypisanie siatki na ekran
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

    def add_h_matrixes_to_elements(self, h_matrix):
        for i in range(len(self.elements)):
            self.elements[i].h_matrix = h_matrix[i]

    def add_c_matrixes_to_elements(self, c_matrix):
        for i in range(len(self.elements)):
            self.elements[i].c_matrix = c_matrix[i]

    def calculate_hbc_matrixes_for_elements(self, hbc_matrix, cond):

        j = 0
        for i in self.elements:

                temp_hbc_matrix = np.zeros((4, 4, 4))
                temp_vect_P = np.zeros((4, 4))

                if i.nodes[0].is_bc == True and i.nodes[1].is_bc == True:
                    temp_hbc_matrix[0] = cond * (np.outer(hbc_matrix[0, 0], hbc_matrix[0, 0].transpose()) + np.outer(hbc_matrix[0, 1], hbc_matrix[0, 1].transpose())) * (i.nodes[1].x - i.nodes[0].x)/2
                    temp_vect_P[0] = cond * (hbc_matrix[0, 0]*self.temp + hbc_matrix[0, 1]*self.temp) * (i.nodes[1].x - i.nodes[0].x)/2

                if i.nodes[1].is_bc == True and i.nodes[2].is_bc == True:
                    temp_hbc_matrix[1] = cond * (np.outer(hbc_matrix[1, 0], hbc_matrix[1, 0].transpose()) + np.outer(hbc_matrix[1, 1], hbc_matrix[1, 1].transpose())) * (i.nodes[2].y - i.nodes[1].y) / 2
                    temp_vect_P[1] = cond * (hbc_matrix[1, 0]*self.temp + hbc_matrix[1, 1]*self.temp) * (i.nodes[2].y - i.nodes[1].y) / 2

                if i.nodes[2].is_bc == True and i.nodes[3].is_bc == True:
                    temp_hbc_matrix[2] = cond * (np.outer(hbc_matrix[2, 0], hbc_matrix[2, 0].transpose()) + np.outer(hbc_matrix[2, 1],hbc_matrix[2, 1].transpose())) * (i.nodes[2].x - i.nodes[3].x) / 2
                    temp_vect_P[2] = cond * (hbc_matrix[2, 0]*self.temp + hbc_matrix[2, 1]*self.temp) * (i.nodes[2].x - i.nodes[3].x) / 2

                if i.nodes[3].is_bc == True and i.nodes[0].is_bc == True:
                    temp_hbc_matrix[3] = cond * (np.outer(hbc_matrix[3, 0], hbc_matrix[3, 0].transpose()) + np.outer(hbc_matrix[3, 1], hbc_matrix[3, 1].transpose())) * (i.nodes[3].y - i.nodes[0].y)/2
                    temp_vect_P[3] = cond * (hbc_matrix[3, 0]*self.temp + hbc_matrix[3, 1]*self.temp) * (i.nodes[3].y - i.nodes[0].y)/2

                i.hbc_matrix = temp_hbc_matrix[0] + temp_hbc_matrix[1] + temp_hbc_matrix[2] + temp_hbc_matrix[3]
                i.p_vect = temp_vect_P[0] + temp_vect_P[1] + temp_vect_P[2] + temp_vect_P[3]

                print(i.p_vect)
                j = j + 1

        k = 0
        for i in self.elements:
            print("element ", k)
            print(i.hbc_matrix)
            k = k + 1


    def calculate_hbc_matrixes_for_elements_3(self, hbc_matrix, cond):

        j = 0
        for i in self.elements:

                temp_hbc_matrix = np.zeros((4, 4, 4))
                temp_vect_P = np.zeros((4, 4))

                if i.nodes[0].is_bc == True and i.nodes[1].is_bc == True:
                    temp_hbc_matrix[0] = cond * (self.weights_3[0]*np.outer(hbc_matrix[0, 0], hbc_matrix[0, 0].transpose()) + self.weights_3[1]*np.outer(hbc_matrix[0, 1], hbc_matrix[0, 1].transpose()) + self.weights_3[2]*np.outer(hbc_matrix[0, 2], hbc_matrix[0, 2].transpose())) * (i.nodes[1].x - i.nodes[0].x)/2
                    temp_vect_P[0] = cond * (self.weights_3[0]*hbc_matrix[0, 0]*self.temp + self.weights_3[1]*hbc_matrix[0, 1]*self.temp+ self.weights_3[2]*hbc_matrix[0, 2]*self.temp) * (i.nodes[1].x - i.nodes[0].x)/2

                if i.nodes[1].is_bc == True and i.nodes[2].is_bc == True:
                    temp_hbc_matrix[1] = cond * (self.weights_3[0]*np.outer(hbc_matrix[1, 0], hbc_matrix[1, 0].transpose()) + self.weights_3[1]*np.outer(hbc_matrix[1, 1], hbc_matrix[1, 1].transpose()) + self.weights_3[2]*np.outer(hbc_matrix[1, 2], hbc_matrix[1, 2].transpose())) * (i.nodes[2].y - i.nodes[1].y) / 2
                    temp_vect_P[1] = cond * (self.weights_3[0]*hbc_matrix[1, 0]*self.temp + self.weights_3[1]*hbc_matrix[1, 1]*self.temp+ self.weights_3[2]*hbc_matrix[1, 2]*self.temp) * (i.nodes[2].y - i.nodes[1].y) / 2

                if i.nodes[2].is_bc == True and i.nodes[3].is_bc == True:
                    temp_hbc_matrix[2] = cond * (self.weights_3[0]*np.outer(hbc_matrix[2, 0], hbc_matrix[2, 0].transpose()) + self.weights_3[1]*np.outer(hbc_matrix[2, 1],hbc_matrix[2, 1].transpose()) + self.weights_3[2]*np.outer(hbc_matrix[2, 2], hbc_matrix[2, 2].transpose())) * (i.nodes[2].x - i.nodes[3].x) / 2
                    temp_vect_P[2] = cond * (self.weights_3[0]*hbc_matrix[2, 0]*self.temp + self.weights_3[1]*hbc_matrix[2, 1]*self.temp+ self.weights_3[2]*hbc_matrix[2, 2]*self.temp) * (i.nodes[2].x - i.nodes[3].x) / 2

                if i.nodes[3].is_bc == True and i.nodes[0].is_bc == True:
                    temp_hbc_matrix[3] = cond * (self.weights_3[0]*np.outer(hbc_matrix[3, 0], hbc_matrix[3, 0].transpose()) + self.weights_3[1]*np.outer(hbc_matrix[3, 1], hbc_matrix[3, 1].transpose()) + self.weights_3[2]*np.outer(hbc_matrix[3, 2], hbc_matrix[3, 2].transpose())) * (i.nodes[3].y - i.nodes[0].y)/2
                    temp_vect_P[3] = cond * (self.weights_3[0]*hbc_matrix[3, 0]*self.temp + self.weights_3[1]*hbc_matrix[3, 1]*self.temp+ self.weights_3[2]*hbc_matrix[3, 2]*self.temp) * (i.nodes[3].y - i.nodes[0].y)/2

                i.hbc_matrix = temp_hbc_matrix[0] + temp_hbc_matrix[1] + temp_hbc_matrix[2] + temp_hbc_matrix[3]
                i.p_vect = temp_vect_P[0] + temp_vect_P[1] + temp_vect_P[2] + temp_vect_P[3]
                print("gowno")
                print(i.p_vect)
                j = j + 1

        k = 0
        for i in self.elements:
            print("element ", k)
            print(i.hbc_matrix)
            k = k + 1


    def aggregation(self):

        for i in range(self.nE):
            for j in range(4):
                self.global_P_matrix[self.elements[i].ID[j] - 1] = self.global_P_matrix[self.elements[i].ID[j] - 1] + self.elements[i].p_vect[j]
                for k in range(4):
                    self.global_H_matrix[self.elements[i].ID[j] - 1, self.elements[i].ID[k] - 1] = self.global_H_matrix[self.elements[i].ID[j] - 1, self.elements[i].ID[k] - 1] + self.elements[i].h_matrix[j, k] + self.elements[i].hbc_matrix[j, k]
                    self.global_C_matrix[self.elements[i].ID[j] - 1, self.elements[i].ID[k] - 1] = self.global_C_matrix[self.elements[i].ID[j] - 1, self.elements[i].ID[k] - 1] + self.elements[i].c_matrix[j, k]

        for i in range(self.nN):
            for j in range(self.nN):
                self.global_H_matrix[i, j] = self.global_H_matrix[i, j] + self.global_C_matrix[i, j]/self.time_step
                self.global_P_matrix[i] = self.global_P_matrix[i] + (self.global_C_matrix[i, j]/self.time_step)*self.nodes_list[i]

        for l in self.global_H_matrix:
            print(l)

        for m in self.global_P_matrix:
            print(m)

        for n in self.global_C_matrix:
            print(n)





