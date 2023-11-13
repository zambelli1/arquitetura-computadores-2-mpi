#include <iostream>
#include <cmath>
#include <mpi.h>

void compute(int total_count, int my_count, float my_points[][3]) {
  // total_count is the total number of points
  // my_count is the number of points for this process
  // my_points is a float table of shape [my_count][3]

  // MPI variables
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // 1- Sum over all the points in local_sum
  float local_sum[3] = {0.0f, 0.0f, 0.0f};
  for (int i = 0; i < my_count; ++i) {
    for (int j = 0; j < 3; ++j) {
      local_sum[j] += my_points[i][j];
    }
  }

  // 2- Reduce the sum of all the points on the variable barycentre
  float barycentre[3];
  MPI_Allreduce(local_sum, barycentre, 3, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

  // 3- Divide every component of the barycentre by the total number of points
  for (int i = 0; i < 3; ++i) {
    barycentre[i] /= total_count;
  }

  // For every point
  for (int i = 0; i < my_count; ++i) {
    float dist = 0.0f;

    // 4- Compute the distance for every point
    for (int j = 0; j < 3; ++j) {
      dist += std::pow(my_points[i][j] - barycentre[j], 2);
    }
    dist = std::sqrt(dist);

    // And printing the result
    std::cout << rank << " " << dist << std::endl;
  }
}

int main_compute(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  // MPI variables
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Example usage of compute function
  const int total_count = 1000;
  const int my_count = total_count / MPI::COMM_WORLD.Get_size();
  float my_points[my_count][3];

  // Fill my_points with some example values
  for (int i = 0; i < my_count; ++i) {
    for (int j = 0; j < 3; ++j) {
      my_points[i][j] = rank * my_count + i + 1; // some example values
    }
  }

  compute(total_count, my_count, my_points);

  MPI_Finalize();

  return 0;
}
