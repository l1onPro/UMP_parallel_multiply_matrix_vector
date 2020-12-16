// UMP_test_1.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include"mpi.h"
#include <vector>
#include <iostream>
#include <chrono>
#include <functional>


using namespace std;

int rank_MPI, size_MPI;


void PrintMatrix(int* matrix, int heigth, int width)
{
	for (int i = 0, k = 1; i < heigth * width; i++, k++)
	{
		std::cout << " " << matrix[i];	
		
		if (k == width)
		{
			k = 0;
			std::cout << endl;			
		}		
	}	
	std::cout << endl;
}
void InitMatrixData(int*& matrix, int size)
{
	for (int i = 0; i < size * size; i++)
	{
		matrix[i] = rand() % 10 + 1;
	}
}

void PrintVector(int* vector, int size)
{	
	for (int i = 0; i < size; i++)
	{
		std::cout << " " << vector[i] << endl;
	}
	std::cout << endl;
}
void InitVectorData(int*& vector, int size)
{
	for (int i = 0; i < size; i++)
	{
		vector[i] = rand() % 10 + 1;
	}
}

void InitParallel(int*& matrix, int*& vector, int*& result, int& size, int*& procRows, int*& procResult, int& RowNum)
{
	if (rank_MPI == 0)
	{
		// Ввод размера матрицы и вектора
		do
		{
			cout << endl;
			cout << "Enter size: ";
			cin >> size;

			if (size < size_MPI) cout << "Size of the objects must be greater than number of processes!" << endl;
			if (size % size_MPI != 0) cout << "Size of objects must be divisible by number of processes!" << endl;

		} while ((size < size_MPI) || (size % size_MPI != 0));
	}
	MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);

	RowNum = size / size_MPI; //кол-во полос матрицы на всех процессах
	vector = new int[size];
	result = new int[size];
	procRows = new int[RowNum * size];
	procResult = new int[RowNum];

	if (rank_MPI == 0)
	{
		matrix = new int[size * size];
		InitMatrixData(matrix, size);
		InitVectorData(vector, size);		
	}
}
void InitLinear(int*& matrix, int*& vector, int*& result, int& size)
{
	do
	{
		cout << endl;
		cout << "Enter size: ";
		cin >> size;

		if (size < 0) cout << "Size of the objects must be greater than 0s!" << endl;
	} while (size < 0);

	vector = new int[size];
	result = new int[size];
	matrix = new int[size * size];

	InitMatrixData(matrix, size);
	InitVectorData(vector, size);
}
void ProcessParallelTermination(int*& matrix, int*& vector, int*& result, int*& procRow, int*& procResult)
{
	if (rank_MPI == 0)
	{
		delete[] matrix;
		delete[] vector;
		delete[] result;
		delete[] procRow;
		delete[] procResult;
	}
}
void ProcessLinearTermination(int*& matrix, int*& vector, int*& result)
{
	delete[] matrix;
	delete[] vector;
	delete[] result;
}

//Распределение данных между процессами
void DataDistribution(int* matrix, int* procRows, int* vector, int size, int RowNum)
{
	MPI_Bcast(vector, size, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(matrix, RowNum * size, MPI_INT, procRows, RowNum * size, MPI_INT, 0, MPI_COMM_WORLD);
}
//Проверка правильности разделения данных между процессами
void TestDistribution(int* matrix, int* vector, int* procRows, int size, int RowNum) {
	
	if (rank_MPI == 0) {
		cout << "Initial Matrix. " << endl; PrintMatrix(matrix, size, size);
		cout << "Initial Vector: " << endl;	PrintVector(vector, size);
		cout << endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i = 0; i < size_MPI; i++) {
		if (rank_MPI == i) {
			cout << "*****************************" << endl;
			cout << "ProcRank = " << rank_MPI << endl;
			cout << "Matrix Stripe:" << endl;			
			PrintMatrix(procRows, RowNum, size);
			cout << "Vector: " << endl;			
			PrintVector(vector, size);	
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
}

//параллельный алгоритм
void ParallelMultiply(int*& procRows, int*& vector, int*& procResult, int size, int RowNum)
{
	for (int i = 0; i < RowNum; i++)
	{
		procResult[i] = 0;
		for (int j = 0; j < size; j++)
		{
			procResult[i] += procRows[i * size + j] * vector[j];
		}
	}
}
//проверка параллельного алгоритма
void TestPartialResults(int* pProcResult, int RowNum) {
	for (int i = 0; i < size_MPI; i++) {
		if (rank_MPI == i) {
			cout << "ProcRank = " << rank_MPI << endl;
			cout << "Part of result vector: " << endl;
			PrintVector(pProcResult, RowNum);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
}
//Последовательный алгоритм
void LinearMultiply(int* matrix, int* vector, int* result, int size)
{	
	for (int i = 0; i < size; i++)
	{
		result[i] = 0;
		for (int j = 0; j < size; j++)
		{
			result[i] += matrix[i * size + j] * vector[j];
		}
	}
}

int main(int argc, char** argv)
{
	int size;
	int* mat;
	int* vec;
	int* result;

	int* procRows; //Полоса матрицы на текущем процессе
	int* procResult; //Полоса матрицы на текущем процессе
	int RowNum; //Количество строк в матричной полосе

	bool parallel = false;
	if (parallel)
	{
		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &size_MPI);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank_MPI);

		if (rank_MPI == 0)
		{
			cout << "Parallel matrix-vector multiply program" << endl;
		}
		
		InitParallel(mat, vec, result, size, procRows, procResult, RowNum);

		auto startTime = std::chrono::high_resolution_clock::now();

		//Распределение данных между процессами
		//DataDistribution(mat, procRows, vec, size, RowNum);	
		MPI_Bcast(vec, size, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Scatter(mat, RowNum * size, MPI_INT, procRows, RowNum * size, MPI_INT, 0, MPI_COMM_WORLD);
		//TestDistribution(mat, vec, procRows, size, RowNum);

		ParallelMultiply(procRows, vec, procResult, size, RowNum);
		//TestPartialResults(procResult, RowNum);
		
		//Сбор результатов в один вектор
		MPI_Allgather(procResult, RowNum, MPI_INT, result, RowNum, MPI_INT, MPI_COMM_WORLD);

		if (rank_MPI == 0)
		{
			auto stopTime = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double, std::milli> executionTime = stopTime - startTime;
			
			//PrintMatrix(mat, size, size);
			//PrintVector(vec, size);
			//PrintVector(result, size);
			
			std::cout << executionTime.count() << " ms" << std::endl;		
		}

		ProcessParallelTermination(mat, vec, result, procRows, procResult);		

		MPI_Finalize();
	}
	else
	{
		cout << "Linear matrix-vector multiply program" << endl;
		
		InitLinear(mat, vec, result, size);		
		
		auto startTime = std::chrono::high_resolution_clock::now();

		LinearMultiply(mat, vec, result, size);	

		auto stopTime = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::milli> executionTime = stopTime - startTime;
		
		//PrintMatrix(mat, size, size);
		//PrintVector(vec, size);
		//PrintVector(result, size);

		std::cout << executionTime.count() << " ms" << std::endl;

		ProcessLinearTermination(mat, vec, result);
	}
	
	
	
	return 0;
}
