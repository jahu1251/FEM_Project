from Element4_2D import Element4_2D
from Grid import Grid
from Jacobian import Jacobian
from Matrix_Integration import Matrix_Integration


def main():

    simulate(4, 4, 0.1, 0.1, 3, 25, 1200, 700, 7800, 300, 50, 100, 500)


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

    matrix_integration_obj = Matrix_Integration(universal_element_obj.derivativeKsi,
                                                universal_element_obj.derivativeEta,
                                                jacobian_obj.calculate_grid_jacobians(), cond, c, ro,
                                                universal_element_obj.N_values)
    matrix_integration_obj.calculate_H_matrixes_for_grid()

    grid_obj.add_h_matrixes_to_elements(matrix_integration_obj.global_h_matrix_list)
    grid_obj.add_c_matrixes_to_elements(matrix_integration_obj.global_c_matrix_list)

    if integ_points_num == 2:
        grid_obj.calculate_hbc_matrixes_for_elements(universal_element_obj.temp_matrix_for_hbc, alfa)
    elif integ_points_num == 3:
        grid_obj.calculate_hbc_matrixes_for_elements_3(universal_element_obj.temp_matrix_for_hbc, alfa)

    grid_obj.aggregation()

    current_time = 0

    while current_time < simulation_time:
        print("czas : ", current_time)
        grid_obj.calculate_p_and_h_matrixes()
        current_time += time_step

    #grid_obj.calculate_p_and_h_matrixes(time_step)


if __name__ == "__main__":
    main()
