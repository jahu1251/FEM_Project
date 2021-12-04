

class Node:
    x, y = None, None

    def __init__(self, x, y, bc, init_temp):
        self.x = x
        self.y = y
        self.is_bc = bc
        self.init_temp = init_temp
