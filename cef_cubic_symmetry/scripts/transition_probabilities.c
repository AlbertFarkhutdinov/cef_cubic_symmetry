#include <stdio.h>
#include <math.h>

float get(float **j_z, float **j_plus, float **j_minus, float **transition_probability, float **eigenfunctions, float j, float squared_j, int size){
	int row, row_j, column;
    for (row = 0; row < size; row++){
        j_z[row][row] += (pow(eigenfunctions[size - 1][row], 2) * (size - 1 - j));
        for (row_j = 0; row_j < (size - 1); row_j++){
            j_z[row][row] += pow(eigenfunctions[row_j][row], 2) * (row_j - j);
            j_plus[row][row] += (eigenfunctions[row_j + 1][row] *
                                 eigenfunctions[row_j][row] *
                                 sqrt(squared_j - (row_j - j) * (row_j - j + 1)));
        }

        j_minus[row][row] = j_plus[row][row];
        for (column = row + 1; column < size; column++){
            float mqn_1 = size - 1 - j;
            j_z[row][column] += (eigenfunctions[size - 1][row] *
                                 eigenfunctions[size - 1][column] * mqn_1);
            for (row_j = 0; row_j < (size - 1); row_j++){
                mqn_1 = row_j - j;
                j_z[row][column] += (eigenfunctions[row_j][row] *
									eigenfunctions[row_j][column] * mqn_1);
                int column_j = row_j + 1;
                float mqn_2 = column_j - j;
                float common_root = sqrt(squared_j - mqn_1 * mqn_2);
                j_plus[row][column] += (eigenfunctions[column_j][row] *
                                        eigenfunctions[row_j][column] *
                                        common_root);
                j_minus[row][column] += (eigenfunctions[row_j][row] *
                                         eigenfunctions[column_j][column] *
                                         common_root);
            }
            transition_probability[row][column] = ((pow(2 * j_z[row][column], 2) +
                                                    pow(j_plus[row][column], 2) +
                                                    pow(j_minus[row][column], 2)) / 3);
            j_z[column][row] = j_z[row][column];
            j_plus[column][row] = j_minus[row][column];
            j_minus[column][row] = j_plus[row][column];
            transition_probability[column][row] = transition_probability[row][column];
        }
    }
    float result[4] = {**j_z, **j_plus, **j_minus, **transition_probability};
    return *result;
}

main()
{
}
