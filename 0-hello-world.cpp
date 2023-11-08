#include <iostream>
#include <sstream>
#include <mpi.h>

int main(int argc, char** argv) {
    int rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Use stringstream para converter o número em uma string
    std::stringstream ss;
    ss << "Hello world, from process #" << rank;
    std::string mensagem = ss.str();

    // Imprima a mensagem
    std::cout << mensagem << std::endl;

    MPI_Finalize();

    return 0;
}
