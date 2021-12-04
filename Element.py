import numpy as np

class Element:
    def __init__(self, id1, id2, id3, id4, node1, node2, node3, node4):
        self.ID = []
        self.ID.append(id1)
        self.ID.append(id2)
        self.ID.append(id3)
        self.ID.append(id4)
        self.nodes = []
        self.nodes.append(node1)
        self.nodes.append(node2)
        self.nodes.append(node3)
        self.nodes.append(node4)
        self.h_matrix = np.zeros((4, 4))
        self.c_matrix = np.zeros((4, 4))
        self.hbc_matrix = np.zeros((4, 4))
        self.p_vect = np.zeros((4))