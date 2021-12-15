from Element4_2D import Element4_2D
from Grid import Grid
from Jacobian import Jacobian
from Matrix_Integration import Matrix_Integration


def main():

    f = open("datafile.txt", "r")
    data = []
    nodes_x = []
    nodes_y = []
    elements_list = []
    bc_list = []
    is_node = False
    is_element = False
    is_bc = False

    for i in f:
        if "*Node" in i:
            is_node = True
            is_element = False
            is_bc = False
            continue
        elif "*Element" in i:
            is_element = True
            is_node = False
            is_element = False
            continue
        elif "*BC" in i:
            is_bc = True
            is_node = False
            is_element = False
            continue
        elif "*Node" not in i and "*Element" not in i and "*BC" not in i:
            is_bc = False
            is_node = False
            is_element = False

        if is_node is True:
            nodes_x.append(i.split(", ")[1])
            nodes_y.append(i.split(", ")[2])
        elif is_element is True:
            for l in i.split(", "):
                elements_list.append(l)
        elif is_bc is True:
            for k in i.split(", "):
                bc_list.append(k)
        else:
            data.append(i.split(" ")[-1])

    simulation_time = data[0].split(" ")[1]
    init_temp = data[5].split(" ")[1]
    time_step = data[1].split(" ")[1]
    alfa = data[3].split(" ")[1]
    ro = data[6].split(" ")[1]
    c = data[7].split(" ")[1]
    temp = data[4].split(" ")[1]
    cond = data[2].split(" ")[1]
    nN = data[8].split(" ")[1]
    nE = data[9].split(" ")[1]

    #simulate(4, 4, 0.1, 0.1, 3, 25, 1200, 700, 7800, 300, 1, 100, 100)

    simulate_from_file(nE, nN, 3, cond, temp, c, ro, alfa, time_step, init_temp, simulation_time, nodes_list, id_list)


def simulate(nH, nB, H, B, integ_points_num, cond, temp, c, ro, alfa, time_step, init_temp, simulation_time):

    print("================ Siatka ================")
    grid_obj = Grid(nH, nB, H, B, temp, time_step, init_temp)
    grid_obj.print_grid()

    # print("================ Całkowanie Gauss'a ================")
    # gauss_obj = Gauss()
    # print("Całka 1d dla 2-punktowego schematu całkowania :", gauss_obj.two_dimensions_gauss_integral(0))
    # print("Całka 1d dla 3-punktowego schematu całkowania :", gauss_obj.two_dimensions_gauss_integral(1))
    # print("Całka 2d dla 2-punktowego schematu całkowania :", gauss_obj.three_dimensions_gauss_integral(0))
    # print("Całka 2d dla 3-punktowego schematu całkowania :", gauss_obj.three_dimensions_gauss_integral(1))

    universal_element_obj = Element4_2D(cond, integ_points_num, c, ro)
    if integ_points_num == 2:
        universal_element_obj.calculate_2()
        universal_element_obj.calculate_temp_matrix_for_hbc_2()
    elif integ_points_num == 3:
        universal_element_obj.calculate_3()
        universal_element_obj.calculate_temp_matrix_for_hbc_3()

    jacobian_obj = Jacobian(universal_element_obj.derivativeKsi, universal_element_obj.derivativeEta, grid_obj.elements)

    if integ_points_num == 2:
        matrix_integration_obj = Matrix_Integration(universal_element_obj.derivativeKsi,
                                                    universal_element_obj.derivativeEta,
                                                    jacobian_obj.calculate_grid_jacobians(), cond, c, ro,
                                                    universal_element_obj.N_values)
        matrix_integration_obj.calculate_H_matrixes_for_grid()
    elif integ_points_num == 3:
        matrix_integration_obj = Matrix_Integration(universal_element_obj.derivativeKsi,
                                                    universal_element_obj.derivativeEta,
                                                    jacobian_obj.calculate_grid_jacobians(), cond, c, ro,
                                                    universal_element_obj.N_values_3)
        matrix_integration_obj.calculate_H_matrixes_for_grid_3()

    grid_obj.add_h_matrixes_to_elements(matrix_integration_obj.global_h_matrix_list)
    grid_obj.add_c_matrixes_to_elements(matrix_integration_obj.global_c_matrix_list)

    if integ_points_num == 2:
        grid_obj.calculate_hbc_matrixes_for_elements(universal_element_obj.temp_matrix_for_hbc, alfa)
    elif integ_points_num == 3:
        grid_obj.calculate_hbc_matrixes_for_elements_3(universal_element_obj.temp_matrix_for_hbc, alfa)

    grid_obj.aggregation()

    current_time = 0

    while current_time < simulation_time:
        print("Czas : ", current_time)
        grid_obj.calculate_p_and_h_matrixes()
        current_time += time_step

def simulate_from_file(nE, nN, integ_points_num, cond, temp, c, ro, alfa, time_step, init_temp, simulation_time, nodes_list, id_list):

    print("================ Siatka ================")
    grid_obj = Grid(3, 3, nE, nN, nodes_list, id_list, temp, time_step, init_temp)
    grid_obj.print_grid()

    # print("================ Całkowanie Gauss'a ================")
    # gauss_obj = Gauss()
    # print("Całka 1d dla 2-punktowego schematu całkowania :", gauss_obj.two_dimensions_gauss_integral(0))
    # print("Całka 1d dla 3-punktowego schematu całkowania :", gauss_obj.two_dimensions_gauss_integral(1))
    # print("Całka 2d dla 2-punktowego schematu całkowania :", gauss_obj.three_dimensions_gauss_integral(0))
    # print("Całka 2d dla 3-punktowego schematu całkowania :", gauss_obj.three_dimensions_gauss_integral(1))

    universal_element_obj = Element4_2D(cond, integ_points_num, c, ro)
    if integ_points_num == 2:
        universal_element_obj.calculate_2()
        universal_element_obj.calculate_temp_matrix_for_hbc_2()
    elif integ_points_num == 3:
        universal_element_obj.calculate_3()
        universal_element_obj.calculate_temp_matrix_for_hbc_3()

    jacobian_obj = Jacobian(universal_element_obj.derivativeKsi, universal_element_obj.derivativeEta, grid_obj.elements)

    if integ_points_num == 2:
        matrix_integration_obj = Matrix_Integration(universal_element_obj.derivativeKsi,
                                                    universal_element_obj.derivativeEta,
                                                    jacobian_obj.calculate_grid_jacobians(), cond, c, ro,
                                                    universal_element_obj.N_values)
        matrix_integration_obj.calculate_H_matrixes_for_grid()
    elif integ_points_num == 3:
        matrix_integration_obj = Matrix_Integration(universal_element_obj.derivativeKsi,
                                                    universal_element_obj.derivativeEta,
                                                    jacobian_obj.calculate_grid_jacobians(), cond, c, ro,
                                                    universal_element_obj.N_values_3)
        matrix_integration_obj.calculate_H_matrixes_for_grid_3()

    grid_obj.add_h_matrixes_to_elements(matrix_integration_obj.global_h_matrix_list)
    grid_obj.add_c_matrixes_to_elements(matrix_integration_obj.global_c_matrix_list)

    if integ_points_num == 2:
        grid_obj.calculate_hbc_matrixes_for_elements(universal_element_obj.temp_matrix_for_hbc, alfa)
    elif integ_points_num == 3:
        grid_obj.calculate_hbc_matrixes_for_elements_3(universal_element_obj.temp_matrix_for_hbc, alfa)

    grid_obj.aggregation()

    current_time = 0

    while current_time < simulation_time:
        print("Czas : ", current_time)
        grid_obj.calculate_p_and_h_matrixes()
        current_time += time_step


if __name__ == "__main__":
    main()
