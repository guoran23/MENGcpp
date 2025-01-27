#ifndef MPIManager_H
#define MPIManager_H

#include <mpi.h>
class MPIManager {
public:
    // 禁止拷贝构造和赋值操作
    MPIManager(const MPIManager&) = delete;
    MPIManager& operator=(const MPIManager&) = delete;

    // 获取唯一实例
    static MPIManager& getInstance(int argc = 0, char** argv = nullptr) {
        static MPIManager instance(argc, argv);
        return instance;
    }

    int getRank() const { return rank; }
    int getSize() const { return size; }

private:
    MPIManager(int argc, char** argv) {
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
    }

    ~MPIManager() { MPI_Finalize(); }

    int rank, size;
};

#endif // MPIManager_H