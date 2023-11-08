#include <iostream>
#include <mpi.h>
#include <cstdlib>

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Lê o valor local do processo
  // local_value conterá um número inteiro específico para o processo 0 e outro para o processo 1
  int local_value;
  local_value = atoi(argv[1]);

  int other_value;
  if (rank == 0) {
    // Aqui, insira o código para o primeiro processo:
    // 1- Envia o valor para o processo 1
    MPI_Send(&local_value, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    // 2- Recebe o valor do processo 1 (em other_value)
    MPI_Recv(&other_value, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // 3- Imprime a soma dos dois valores na saída padrão
    std::cout << "Soma: " << local_value + other_value << std::endl;
  }
  else {
    // Aqui, insira o código para o segundo processo:
    // 1- Recebe o valor do processo 0 (em other_value)
    MPI_Recv(&other_value, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // 2- Envia o valor para o processo 0
    MPI_Send(&local_value, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    // 3- Imprime o produto dos dois valores na saída padrão
    std::cout << "Produto: " << local_value * other_value << std::endl;
  }

  MPI_Finalize();

  return 0;
}
