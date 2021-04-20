#include <stdio.h> //printf
#include <stdlib.h>//srand, rand
#include <time.h>  //time

const unsigned BLOCKS_PER_GRID = 5;
const unsigned THREADS_PER_BLOCK = 20;

template <typename T>
class Matrix{
	unsigned M,N; //rozmiary macierzy
	T** tab;      //wskaznik na dane macierzy
public:
	Matrix(unsigned M, unsigned N){
		this->M = M;
		this->N = N;
		tab = new T*[M];
		for(unsigned i=0;i<M;i++){
			tab[i] = new T[N];
			for(unsigned j=0;j<N;j++){
				tab[i][j] = (i==j);
			}
		}
	}
	
	~Matrix(){
		for(unsigned i=0;i<M;i++){
			delete[] tab[i];	
		}
		delete[] tab;
	}
	unsigned getRows() const{ return M;}
	unsigned getColumns() const{ return N;}
	T get(unsigned m, unsigned n) const{
		if(m<M && n<N) return tab[m][n];
		throw "Nie znaleziono elementu pod podana para indeksow.";
	}
	void set(unsigned m, unsigned n, T num){
		if(m<M && n<N) tab[m][n] = num;
		else throw "Nie mozna ustawic wartosci elementu pod podana para indeksow.";
	}
};

template<typename T>
void printMatrixCPU(const char* format, const Matrix<T>& matrix){
	for(unsigned i=0; i<matrix.getRows(); i++){
		for(unsigned j=0; j<matrix.getColumns(); j++){
			printf(format, matrix.get(i,j));
		}
		printf("\n");
	}
}
template<typename T>
void printMatrixGPU(const char* format, void* matrix){
	for(unsigned i=0; i<*(unsigned*)matrix; i++){  //wartosc *(unsigned*)matrix przechowuje liczbe wierszy
		for(unsigned j=0; j<*((unsigned*)matrix+1); j++){ //wartosc *(unsigned*)matrix+1 przechowuje liczbe kolumn 
			//duzo castowania
			printf(format, *((T*)((unsigned*)matrix+2)+i*(*(unsigned*)matrix)+j));
		}
		printf("\n");
	}
}

int main(){
	srand(time(NULL));
	//reserve matrix of ints
	const unsigned M = 2; //test values
	const unsigned N = 5;
	Matrix<int> matrixInt(M,N);
	for(unsigned i=0; i<M; i++){
		for(unsigned j=0; j<N; j++){
			matrixInt.set(i, j, 2*i+j);				
		}
	}
	printMatrixCPU<int>("%d ", matrixInt);

	//copy of the matrix in device memory
	void* d_ptr_matrixInt;
	cudaMallocManaged(&d_ptr_matrixInt, sizeof(matrixInt));
	//cudaMemcpy((void*)d_ptr_matrixInt, (void*)&matrixInt, sizeof(matrixInt), cudaMemcpyDefault);						
	//Serializacja	
	*((unsigned*)d_ptr_matrixInt/*+0*sizeof(unsigned)*/) = matrixInt.getRows();
	*((unsigned*)d_ptr_matrixInt+1/**sizeof(unsigned)*/) = matrixInt.getColumns();
	for(unsigned i=0;i<matrixInt.getRows();i++){
		for(unsigned j=0;j<matrixInt.getColumns();j++){
			*((int*)((unsigned*)d_ptr_matrixInt+2)+(i*matrixInt.getRows()+j)) = matrixInt.get(i,j);
			//*((int*)d_ptr_matrixInt+2*sizeof(unsigned)+(i*matrixInt.getRows()+j)*sizeof(int)) = matrixInt.get(i,j);
		}
	}

	printMatrixGPU<int>("%d ", d_ptr_matrixInt);
	printf("%d\n", sizeof(unsigned));
	printf("%d\n", sizeof(int));
	cudaFree(d_ptr_matrixInt);
	return 0;
}
