from Element4_2D import Element4_2D
from Grid import Grid
from Jacobian import Jacobian
from Matrix_Integration import Matrix_Integration


def main():
    simulate(5, 5, 0.1, 0.1, 2, 25, 1200)


def simulate(nH, nB, H, B, integ_points_num, cond, temp):
    print("================ Siatka ================")
    grid_obj = Grid(nH, nB, H, B, temp)
    grid_obj.print_grid()

    # print("================ Całkowanie Gauss'a ================")
    # gauss_obj = Gauss()
    # print("Całka 1d dla 2-punktowego schematu całkowania :", gauss_obj.two_dimensions_gauss_integral(0))
    # print("Całka 1d dla 3-punktowego schematu całkowania :", gauss_obj.two_dimensions_gauss_integral(1))
    # print("Całka 2d dla 2-punktowego schematu całkowania :", gauss_obj.three_dimensions_gauss_integral(0))
    # print("Całka 2d dla 3-punktowego schematu całkowania :", gauss_obj.three_dimensions_gauss_integral(1))

    print("================ Pochodne ================")

    universal_element_obj = Element4_2D(cond, integ_points_num)
    if integ_points_num == 2:
        universal_element_obj.calculate_2()
        universal_element_obj.calculate_temp_matrix_for_hbc_2()
    elif integ_points_num == 3:
        universal_element_obj.calculate_3()
        universal_element_obj.calculate_temp_matrix_for_hbc_3()

    print("================ Jacobian ================")
    jacobian_obj = Jacobian(universal_element_obj.derivativeKsi, universal_element_obj.derivativeEta, grid_obj.elements)

    print("================ Macierz ================")
    matrix_integration_obj = Matrix_Integration(universal_element_obj.derivativeKsi, universal_element_obj.derivativeEta, jacobian_obj.calculate_grid_jacobians(), cond)
    matrix_integration_obj.calculate_H_matrixes_for_grid()

    grid_obj.add_h_matrixes_to_elements(matrix_integration_obj.global_h_matrix_list)
    grid_obj.calculate_hbc_matrixes_for_elements(universal_element_obj.temp_matrix_for_hbc, cond)

    grid_obj.aggregation()

if __name__ == "__main__":
    main()
