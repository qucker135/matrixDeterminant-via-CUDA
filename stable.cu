#include <stdio.h> //printf
#include <stdlib.h>//srand, rand
#include <time.h>  //time
#include <fstream>
#include <string>

#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>

//globalne ustawienia
const unsigned BLOCKS_PER_GRID = 16;
const unsigned THREADS_PER_BLOCK = 256;

class Matrix{
	unsigned M,N; //rozmiary macierzy
	double** tab;      //wskaznik na dane macierzy
	Matrix minor(unsigned m,unsigned n){
		Matrix result(M-1,N-1);
		for(unsigned i=0;i<M;i++){
			for(unsigned j=0;j<N;j++){
				if(i<m && j<n) result.set(i,j,this->get(i,j));
				else if(i<m && j>n) result.set(i,j-1,this->get(i,j));
				else if(i>m && j<n) result.set(i-1,j,this->get(i,j));
				else if(i>m && j>n) result.set(i-1,j-1,this->get(i,j));
			}
		}
		return result;
	}
public:
	Matrix(unsigned M, unsigned N){
		this->M = M;
		this->N = N;
		tab = new double*[M];
		for(unsigned i=0;i<M;i++){
			tab[i] = new double[N];
			for(unsigned j=0;j<N;j++){
				tab[i][j] = (i==j); //domyslnie macierz jednostkowa
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

	double get(unsigned m, unsigned n) const{
		if(m<M && n<N) return tab[m][n];
		throw "Nie znaleziono elementu pod podana para indeksow.";
	}

	void set(unsigned m, unsigned n, double num){
		if(m<M && n<N) tab[m][n] = num;
		else throw "Nie mozna ustawic wartosci elementu pod podana para indeksow.";
	}

	Matrix operator+(const Matrix& mat) const{
		if(this->M != mat.getRows() || this->N != mat.getColumns()) throw "Nie mozna dodawac macierzy o roznych wymiarach!";
		Matrix result = Matrix(M,N);
		for(unsigned i=0;i<M;i++){
			for(unsigned j=0;j<N;j++){
				result.set(i,j,this->get(i,j)+mat.get(i,j));
			}
		}
		return result;
	}

	Matrix operator*(const Matrix& mat) const{
		if(this->N != mat.getRows()) throw "Pierwsza macierz nie ma tylu kolumn, co druga wierszy - odmowa wykonania mnozenia!";
		Matrix result = Matrix(M, mat.getColumns());
		for(unsigned i=0;i<M;i++){
			for(unsigned j=0;j<mat.getColumns();j++){
				//result.set(i,j,this->get(i,j)+mat.get(i,j));
				result.set(i,j,0.0);
				for(unsigned k=0;k<this->N;k++){
					result.set(i,j,result.get(i,j)+this->get(i,k)*mat.get(k,j));
				}
			}
		}
		return result;
	}

	double det(){
		if(M != N) throw "Wyznacznik mozna obliczyc tylko dla macierzy kwadratowej!";
		if(M == 1) return get(0,0);
		else{
			double result = 0.0;
			for(unsigned i=0;i<M;i++){
				//result += ((i%2==0)?(1.0):(-1.0))*get(i,0)*minor(i,0).det();
				Matrix minr = minor(i,0);
				//printf("%.1f \n", get(i,0));
				result += ((i%2==0)?(1.0):(-1.0))*get(i,0)*minr.det();
			}
			return result;
		}
	}

	Matrix inverse(){
		double mainDet = det();
		if(mainDet == 0.0 || mainDet == -0.0) throw "Wyznacznik macierzy rowny zero!";
		Matrix result = Matrix(M,N);
		if(M==1) result.set(0,0,1.0/mainDet);
		else{
			for(unsigned i=0;i<M;i++){
				for(unsigned j=0;j<N;j++){
					double detElement = (((i%2)+(j%2))%2==0?1.0:-1.0)*minor(i,j).det();

					result.set(j,i,detElement/mainDet);
				}
			}
		}
		return result;

	}
  	void setSize(unsigned M, unsigned N) {
		for(int i =0 ;i<this->M;i++){
			delete[] tab[i];
		}
		delete[] tab;

		this->M = M;
		this->N = N;

		this->tab = new double* [this->M];
		for (unsigned i = 0; i < this->M; i++) {
			tab[i] = new double[this->N];
			for(unsigned j=0;j< this->N; j++){
				tab[i][j] = (i==j); //domyslnie macierz jednostkowa	
			}
		}
	}
	/*
	Matrix& operator=(const Matrix& m){
		this->setSize(m.getRows(),m.getColumns());
		for(unsigned i=0;i<m.getRows();i++){
			for(unsigned j=0;j<m.getColumns();j++){
				this->set(i,j,m.get(i,j));
			}
		}
		return &m;
	}
	*/

};

void printMatrixCPU(const char* format, const Matrix& matrix){
	for(unsigned i=0; i<matrix.getRows(); i++){
		for(unsigned j=0; j<matrix.getColumns(); j++){
			printf(format, matrix.get(i,j));
		}
		printf("\n");
	}
}
void printMatrixGPU(const char* format, void* matrix){
	//debug
	printf("%d\n",*(unsigned*)matrix);
	printf("%d\n",*((unsigned*)matrix+1));
	
	for(unsigned i=0; i<*(unsigned*)matrix; i++){  //wartosc *(unsigned*)matrix przechowuje liczbe wierszy
		for(unsigned j=0; j<*((unsigned*)matrix+1); j++){ //wartosc *(unsigned*)matrix+1 przechowuje liczbe kolumn
			//duzo castowania
			printf(format, *((double*)((unsigned*)matrix+2)+i*(*((unsigned*)matrix+1))+j));
		}
		printf("\n");
	}
}

void copyMatrixToGPU( void* ptr, const Matrix& matrix){
	//Kopia macierzy w pamięci GPU
	*((unsigned*)ptr) = matrix.getRows();
	*((unsigned*)ptr+1) = matrix.getColumns();
	for(unsigned i=0;i<matrix.getRows();i++){
		for(unsigned j=0;j<matrix.getColumns();j++){
			//*((double*)((unsigned*)ptr+2)+(i*matrix.getRows()+j)) = matrix.get(i,j);
			*((double*)((unsigned*)ptr+2)+(i*matrix.getColumns()+j)) = matrix.get(i,j);
		}
	}
}

__global__
void add(void* matrix1, void* matrix2, void* result){//PRZY ZALOZENIU, ZE PAMIEC jest juz zaalokowana, w a matrix1 i matrix 2 są juz gotowe wartosci, zakladamy tez, ze zgadzaja sie wymiary macierzy i ewentualne inne warunki, NOTE: wymiary macierzy trzeba wpisac poza kernelem
	unsigned idx = threadIdx.x + blockDim.x * blockIdx.x;
	unsigned stride = blockDim.x * gridDim.x;

	for(unsigned i = idx; i<(*(unsigned*)matrix1)*(*((unsigned*)matrix1+1)); i+=stride){
		*((double*)((unsigned*)result+2)+i) = *((double*)((unsigned*)matrix1+2)+i) + *((double*)((unsigned*)matrix2+2)+i); 
	}
}

__global__
void mul(void* matrix1, void* matrix2, void* result){
	unsigned idx = threadIdx.x + blockDim.x * blockIdx.x;
	unsigned stride = blockDim.x * gridDim.x;

	
	for(unsigned i = idx; i<(*(unsigned*)matrix1)*(*((unsigned*)matrix2+1)); i+=stride){ //przechodzimy po elementach macierzy wynikowej
		double product = 0.0;
		for(unsigned j=0; j<*(unsigned*)matrix2; j++){
			product += (*((double*)((unsigned*)matrix1+2)+(i/(*((unsigned*)matrix2+1)))*(*((unsigned*)matrix1+1))+j))
				*  (*((double*)((unsigned*)matrix2+2)+j*(*((unsigned*)matrix2+1))+i%(*((unsigned*)matrix2+1))));
		}
		*((double*)((unsigned*)result+2)+i) = product;
	}
}


//BLOCKS_PER_GRID * THREADS_PER_BLOCK >= N!

__global__
void detHelper(void* matrix, double* tab_helper/*, unsigned N*/){
	unsigned idx = threadIdx.x + blockDim.x * blockIdx.x;
	unsigned stride = blockDim.x * gridDim.x;
	unsigned N = *(unsigned*)matrix;

	//zapisujemy idx w silniowym systemie liczbowym:
	int* idxTab = new int[N];
	for(int i=1;i<=N;i++){
		idxTab[N-i] = idx%i;
		idx/=i;
	}
	//upewniamy sie, ze udzial w obliczeniach biora tylko te watki, dla ktorych oryginalny index nalezal do przedzialu <0, N!), czyli ze mamy co najwyzej N! watkow wykonujacych jakas prace
	if(idx==0){
		//zapisujemy stride (skok dla pojedynczego watku) w silniowym systemie liczbowym
		int* strideTab = new int[N];
		for(int i=1;i<=N;i++){
			strideTab[N-i] = stride%i;
			stride/=i;
		}

		//przywracamy oryginalne wartosci
		idx = threadIdx.x + blockDim.x * blockIdx.x;
		stride = blockDim.x * gridDim.x;

		//warunek konca pracy watku (wyjscie poza zakres na skutek stalego dodawania strideTab do idxTab)
		while(idxTab[0]<N){
			//wskaznik parzystosci permutacji, 0 - parzysta, 1 - nieparzysta
			int parz = 0;
			for(int i=0;i<N;i++) parz = (parz + idxTab[i])%2;
			//konwersja na interesujaca permutacje
			for(int i=0;i<N;i++){
				for(int j=1;j<=N;j++){
					bool niepojawilo = true;
					for(int k=0;k<i;k++){
						if(idxTab[k] == j){niepojawilo = false; break;}
					}
					if(niepojawilo){
						if(idxTab[i]==0){
							idxTab[i]=j;break;
						}
						else idxTab[i]--;
					}
				}
			}
			//iloczyn czastkowy (jeden z N! z wzoru na wyznacznik), z odpowiednim znakiem
			double product = ((parz%2==0) ? 1.0 : (-1.0));

			for(int i=0;i<N;i++) product*= *(((double*)((unsigned*)matrix+2))+i*N+(idxTab[i]-1));//element z i-tego wiersza i (idxTab[i]-1) kolumny

			//product*=M[i*N+(idxTab[i]-1)+1]; //here we have a product, one of N!

			tab_helper[idx] += product;

			//printf("%lf\n",tab_helper[idx]);

			//konwersja odwrotna (przywrocenie numeru indexu)
			for(int i=0;i<N;i++){
				int ile = 0;
				for(int j=i+1;j<N;j++){
					if(idxTab[j]<idxTab[i]) ile++;
				}
				idxTab[i] = ile;
			}
			//dodanie idxTab+=strideTab, (patrz: warunek konca petli)
			int ak=0;
			for(int i=1;i<=N;i++){
				idxTab[N-i]=idxTab[N-i]+strideTab[N-i]+ak;
				ak=idxTab[N-i]/i;
				if(i!=N) idxTab[N-i]%=i;
			}
			idxTab[0]+=ak;


		}
		delete[] strideTab;
	}
	delete[] idxTab;
}

//do usuniecia
/*
void det(void* matrix, double* tab_helper, double* wsk){ //tab_helper - tablica czastkowych wynikow z pracy wszystkich watkow, wsk - wynik koncowy
	//pobranie wymiarow macierzy
	unsigned M = *(unsigned*)matrix;
	unsigned N = *((unsigned*)matrix+1);
	//jesli macierz nie jest kwadratowa, rzuc wyjatkiem:
	if(M!=N) throw "Macierz niekwadratowa";
	detHelper<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(matrix, tab_helper, N);
	cudaDeviceSynchronize();
	(*wsk)=0.0;// do usuniecia
	//debug
	for(unsigned i=0;i<BLOCKS_PER_GRID*THREADS_PER_BLOCK;i++) printf("%f ",tab_helper[i]);
	for(unsigned i=0;i<BLOCKS_PER_GRID*THREADS_PER_BLOCK;i++) *wsk+=tab_helper[i];
}
*/

//do usuniecia
/*
__global__
void add(void* matrix1, void* matrix2, void* prod, int M, int N){
		const int tidx = blockDim.x * blockIdx.x + threadIdx.x;
		const int tidy = blockDim.y * blockIdx.y + threadIdx.y;

		if(tidx<M && tidy<N){
		*(((double*)(unsigned*)prod+2)+tidx + tidy*N ) = *(((double*)(unsigned*)matrix1+2)+tidx + tidy*N ) + *(((double*)(unsigned*)matrix2+2)+tidx + tidy*N );
		}
}


//do usuniecia

__global__
void mull(void* A, void* B, void* C, int N) {

  int col = blockIdx.x * blockDim.x + threadIdx.x;
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	double sum = 0;
	if (row < N && col < N) {
		for (int i = 0; i < N; i++){
		 sum += *(((double*)(unsigned*)A+2)+row*N+i) * *(((double*)(unsigned*)B+2)+i*N+col);
		}
		*(((double*)(unsigned*)C+2)+row*N+col) = sum;
	}
}
*/


void loadMatrixFromFile(Matrix & matrix) {
	char filename[50];
	printf("Podaj nazwe pliku: ");
	scanf("%s", filename);

	std::ifstream myfile;
	myfile.open(filename);
	
	if (!myfile.good()){
		printf("Nie udalo sie otworzyc pliku o podanej nazwie!\n");
	}
	if (myfile.is_open())
	{
		std::string line;
		getline(myfile, line);
		unsigned M = stoul(line);
		getline(myfile, line);
		unsigned N = stoul(line);
		matrix.setSize(M, N);
			for (unsigned i = 0; i < M; i++) {
				for (unsigned j = 0; j < N; j++) {
					getline(myfile, line);
					matrix.set(i, j, stod(line));
				}
			}

		myfile.close();
		//printMatrixCPU("%.1f ", matrix);
	}
}

void saveMatrixToFile(Matrix& matrixInt) {
	char filename[50];
	printf("Podaj nazwe pliku: ");
	scanf("%s", filename);

	std::ofstream myfile;
	myfile.open(filename);
	myfile << matrixInt.getRows()<<"\n";
	myfile << matrixInt.getColumns()<<"\n";
	for (int i = 0; i < matrixInt.getRows(); i++) {
		for (int j = 0; j < matrixInt.getColumns(); j++) {
			myfile << matrixInt.get(i, j)<< "\n";
		}
	}
	myfile.close();
}

//do konsultacji, i prawdopodobnie usuniecia (bo mamy copyMatrixToGPU)
/*
void executeAdding(Matrix matrix1, Matrix matrix2){
	//Pamięć potrzebna do przechowania macierzy
		size_t sizeOfMatrix1 = sizeof(matrix1) + matrix1.getRows() * matrix1.getColumns() * sizeof(double);
		size_t sizeOfMatrix2 = sizeof(matrix2) + matrix2.getRows() * matrix2.getColumns() * sizeof(double);

		//host pointers
		void* h_ptr_matrixInt = malloc(sizeOfMatrix1);
		void* h_ptr_matrixInt2 = malloc(sizeOfMatrix2);
		void* h_ptr_product = malloc(sizeOfMatrix1);

		//zapisanie macierzy w postaci przystępnej dla GPU
		copyMatrixToGPU(h_ptr_matrixInt, matrix1);
		copyMatrixToGPU(h_ptr_matrixInt2, matrix2);

		//device pointers
		void* d_ptr_matrixInt;
		void* d_ptr_matrixInt2;
		void* d_ptr_product;

		cudaMalloc(&d_ptr_matrixInt, sizeOfMatrix1);
		cudaMalloc(&d_ptr_matrixInt2, sizeOfMatrix1);
		cudaMalloc(&d_ptr_product, sizeOfMatrix1);

		cudaMemcpy(d_ptr_matrixInt, h_ptr_matrixInt, sizeOfMatrix1, cudaMemcpyHostToDevice);
		cudaMemcpy(d_ptr_matrixInt2, h_ptr_matrixInt2, sizeOfMatrix1, cudaMemcpyHostToDevice);

		const size_t BLOCK_DIM = 16;
		dim3 dimBlock(BLOCK_DIM, BLOCK_DIM);
		dim3 dimGrid((int)ceil((double)matrix1.getColumns()/dimBlock.x),(int)ceil((double)matrix1.getRows()/dimBlock.y));
		add<<<dimGrid, dimBlock>>>(d_ptr_matrixInt, d_ptr_matrixInt2, d_ptr_product, matrix1.getRows(), matrix1.getColumns());
		cudaMemcpy(h_ptr_product, d_ptr_product, sizeOfMatrix1, cudaMemcpyDeviceToHost);

		*(unsigned*)h_ptr_product = matrix1.getRows();
		*((unsigned*)h_ptr_product+1) = matrix1.getColumns();
		printf("%u \n",*((unsigned*)h_ptr_product+1));
		printMatrixGPU("%.1lf ", h_ptr_matrixInt);
		printf("\n\n");
		printMatrixGPU("%.1lf ", h_ptr_matrixInt2);
			printf("\n\n");
		printMatrixGPU("%.1lf ", h_ptr_product);
		cudaFree(d_ptr_matrixInt);
		cudaFree(d_ptr_matrixInt2);
		cudaFree(d_ptr_product);

		free(h_ptr_matrixInt);
		free(h_ptr_matrixInt2);
		free(h_ptr_product);
}
*/

//do konsultacji i prawdopodobnie usuniecia (bo mamy copyMatrixToGPU)

/*
void executeMultiplying(Matrix matrix1, Matrix matrix2){
		Matrix prodHost = matrix1*matrix2;
		printMatrixCPU("%.1lf ", prodHost);

		Matrix productMatrix(matrix1.getRows(), matrix1.getColumns());
	//Pamięć potrzebna do przechowania macierzy
		size_t sizeOfMatrix1 = sizeof(matrix1) + matrix1.getRows() * matrix1.getColumns() * sizeof(double);
		size_t sizeOfMatrix2 = sizeof(matrix2) + matrix2.getRows() * matrix2.getColumns() * sizeof(double);
		size_t sizeOfProdMatrix = sizeof(matrix2) + productMatrix.getRows() * productMatrix.getColumns() * sizeof(double);

		//host pointers
		void* h_ptr_matrix1 = malloc(sizeOfMatrix1);
		void* h_ptr_matrix2 = malloc(sizeOfMatrix2);
		void* h_ptr_product = malloc(sizeOfProdMatrix);

		//zapisanie macierzy w postaci przystępnej dla GPU
		copyMatrixToGPU(h_ptr_matrix1, matrix1);
		copyMatrixToGPU(h_ptr_matrix2, matrix2);

		//device pointers
		void* d_ptr_matrix1;
		void* d_ptr_matrix2;
		void* d_ptr_product;

		cudaMalloc(&d_ptr_matrix1, sizeOfMatrix1);
		cudaMalloc(&d_ptr_matrix2, sizeOfMatrix2);
		cudaMalloc(&d_ptr_product, sizeOfProdMatrix);

		cudaMemcpy(d_ptr_matrix1, h_ptr_matrix1, sizeOfMatrix1, cudaMemcpyHostToDevice);
		cudaMemcpy(d_ptr_matrix2, h_ptr_matrix2, sizeOfMatrix2, cudaMemcpyHostToDevice);

		const size_t BLOCK_DIM = 16;
		dim3 dimBlock(BLOCK_DIM, BLOCK_DIM);
		dim3 dimGrid((int)ceil((double)matrix1.getRows()/dimBlock.x),(int)ceil((double)matrix1.getColumns()/dimBlock.y));

		int Arows = matrix1.getRows();
		// int Acols = matrix1.getColumns();
		// int Brows = matrix2.getRows();
		// int Bcols = matrix2.getColumns();
		// int Prows = productMatrix.getRows();
		// int Pcols = productMatrix.getColumns();
		mull<<<dimGrid, dimBlock>>>(d_ptr_matrix1, d_ptr_matrix2, d_ptr_product, Arows);
		cudaMemcpy(h_ptr_product, d_ptr_product, sizeOfProdMatrix, cudaMemcpyDeviceToHost);

		*(unsigned*)h_ptr_product = matrix1.getRows();
		*((unsigned*)h_ptr_product+1) = matrix1.getColumns();

		printf("%u \n",*((unsigned*)h_ptr_product+1));
		printMatrixGPU("%.1lf ", h_ptr_matrix1);
		printf("\n\n");
		printMatrixGPU("%.1lf ", h_ptr_matrix2);
			printf("\n\n");
		printMatrixGPU("%.1lf ", h_ptr_product);
		cudaFree(d_ptr_matrix1);
		cudaFree(d_ptr_matrix2);
		cudaFree(d_ptr_product);

		free(h_ptr_matrix1);
		free(h_ptr_matrix2);
		free(h_ptr_product);
}
*/

/*
void startMenu() {
	//reserve matrix of ints
	const unsigned M = 3; //test values
	const unsigned N = 3;
	Matrix matrix1(M,N);
 	Matrix matrix2(M,N);

	for(int i = 0 ;i<M;i++){
		for(int j = 0 ; j<N ; j++){
			matrix1.set(i,j,i);
		}
	}
	for(int i = 0 ;i<M;i++){
		for(int j = 0 ; j<N ; j++){
			matrix2.set(i,j,j);
		}
	}
	int x = 11;
	while (x != 0)
	{
		printf( "*******************MENU********************\n");
		printf("1.Wczytaj macierz z pliku\n");
		printf("2.Zapisz macierz do pliku\n");
		printf("3.Wyznacznik\n");
		printf("4.Dodaj macierze\n");
		printf("5.Przemnóż macierze\n");
		printf("6.Macierz odwrotna\n");
    		printf("0.Wyjdź\n");

		printf("Wybieram : \n");
		scanf("%d", &x);
		system("clear");

    		int subselect= 0;
		switch (x){
		case 1: //ladowanie
			printf("1 - Macierz 1\n2 - Macierz 2\n");
			scanf("%d", &subselect);
			subselect==1?loadMatrixFromFile(matrix1):loadMatrixFromFile(matrix2);
		break;
		case 2: //zapisywanie
			printf("1 - Macierz 1\n2 - Macierz 2\n");
			scanf("%d", &subselect);
			subselect==1?saveMatrixToFile(matrix1):saveMatrixToFile(matrix2);
		break;
		case 3: //wyznacznik
		//	detOfMatrix(matrixInt);
			break;
		case 4: //dodawanie
			if(matrix1.getRows() != matrix2.getRows() || matrix1.getColumns()!=matrix2.getColumns()){
				printf("Macierze różnych rozmiarów!\n");
			}
			//else executeAdding(matrix1, matrix2);
		break;
		case 5: //mnozenie
			if(matrix1.getColumns() != matrix2.getRows()){
				printf("\n");
				break;
			}
			//else executeMultiplying(matrix1, matrix2);
		break;
    		case 6: //macierz odwrotna

      		break;
    		case 0: //wyjscie
      			exit(0);
    		default:
      		break;

		}
	}
}
*/

using namespace std::chrono;

int main(){

	// 	Matrix test(3,3);
// 	for(int i = 0 ;i<test.getRows();i++){
// 		for(int j=0; j< test.getColumns(); j++){
// 			test.set(i,j,i*2+j+1);
// 		}
// 	}
// void* pointerToTest;
// void* d_pointerToTest;
// pointerToTest = malloc(9*sizeof(double)+2*sizeof(unsigned));
//
// 	copyMatrixToGPU(pointerToTest, test);
// 	printMatrixGPU("%.1lf ", pointerToTest);
//
// 	cudaMalloc(&d_pointerToTest, 9*sizeof(double)+2*sizeof(unsigned));
//
// 	cudaMemcpy(d_pointerToTest, pointerToTest, 9*sizeof(double)+2*sizeof(unsigned), cudaMemcpyHostToDevice);
//
// 	const size_t BLOCK_DIM = 16;
// 	dim3 dimBlock(BLOCK_DIM, BLOCK_DIM);
// 	dim3 dimGrid((int)ceil((double)test.getRows()/dimBlock.x),(int)ceil((double)test.getColumns()/dimBlock.y));
//
// 	double* tab_helper = new double[9];
// 	detHelper<<<dimGrid, dimBlock>>>(d_pointerToTest, tab_helper, 3);
// 	cudaDeviceSynchronize();





	/*
	double* tab;// = new int[BLOCKS_PER_GRID*THREADS_PER_BLOCK];
	unsigned N =10;
		 high_resolution_clock::time_point t11 = high_resolution_clock::now();
		cudaMallocManaged(&tab, sizeof(double)*BLOCKS_PER_GRID*THREADS_PER_BLOCK);
		for(int i=0;i<BLOCKS_PER_GRID*THREADS_PER_BLOCK;i++) tab[i] = 0;

		void* A;

		cudaMallocManaged(&A, (N*N)*sizeof(double)+2*sizeof(unsigned));

		*(unsigned*)A = N;
		*((unsigned*)A+1)  = N;

		for(int i=0;i<*(unsigned*)A;i++){
			for(int j=0;j<*(unsigned*)A ;j++){
				*((double*)((unsigned*)A+2)+i*N+j) = ((i==j) ? 3.0 : 2.0);//rand()%21-10;
			}
		}
		double* w = new double;
		 *w = 0;
		detHelper<<<BLOCKS_PER_GRID,THREADS_PER_BLOCK>>>(A,tab,N);
		cudaDeviceSynchronize();
		high_resolution_clock::time_point t22 = high_resolution_clock::now();

	  duration<double> time_span1 = duration_cast<duration<double>>(t22 - t11);

	  std::cout << "It took me " << time_span1.count() << " seconds.";
	  std::cout << std::endl;



		Matrix test(N,N);
		for(int i = 0 ;i< N; i++){
			for(int j =0 ; j<N;j++){
				test.set(i,j,(i==j) ? 3.0 : 2.0);
			}
		}

 		 high_resolution_clock::time_point t1 = high_resolution_clock::now();
		 test.det();
		 high_resolution_clock::time_point t2 = high_resolution_clock::now();

		 duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

		 std::cout << "(host)It took me " << time_span.count() << " seconds.";
		 std::cout << std::endl;

		for(int i=0;i<BLOCKS_PER_GRID*THREADS_PER_BLOCK;i++) (*w)+=tab[i];
		printf("%lf\n",*w);

		cudaFree(A);
		cudaFree(tab);
		delete w;
	*/
	/*
	const unsigned Ma = 2;
	const unsigned Na = 3;
	Matrix A(Ma,Na);
	for(unsigned i=0;i<M;i++){
		for(unsigned j=0;j<N;j++){
			A.set(i,j,M-i);	
		}
	}
	A.set(0,0,53);
	A.set(0,1,-3);
	A.set(0,2,2.5);
	//A.set(0,3,1.3);
	//A.set(0,4,1.2);
	A.set(1,0,-53);
	A.set(1,1,4);
	A.set(1,2,1.5);
	//A.set(1,3,0.3);
	//A.set(1,4,-1.2);

	const unsigned Mb = Na;
	const unsigned Nb = 5;

	Matrix B(Mb,Nb);
	for(unsigned i=0;i<M;i++){
		for(unsigned j=0;j<N;j++){
			B.set(i,j,N-j);	
		}
	}

	B.set(0,0,17);
	B.set(0,1,15);
	B.set(0,2,2.5);
	B.set(0,3,0.3);
	B.set(0,4,0.2);
	B.set(1,0,-5);
	B.set(1,1,4.1);
	B.set(1,2,2.5);
	B.set(1,3,1.3);
	B.set(1,4,-0.2);
	B.set(2,0,-1);
	B.set(2,1,41);
	B.set(2,2,25);
	B.set(2,3,13);
	B.set(2,4,0.2);


	printMatrixCPU("%.1f ",A);
	printf("\n");
	printMatrixCPU("%.1f ",B);
	printf("\n");

	if(A.getColumns() != B.getRows()){
		throw "Nierówne wymiary macierzy!!";
	}


	void* d_A = NULL;
	void* d_B = NULL;
	void* result = NULL;

	
	cudaMallocManaged(&d_A, 2*sizeof(unsigned)+Ma*Na*sizeof(double));
	cudaMallocManaged(&d_B, 2*sizeof(unsigned)+Mb*Nb*sizeof(double));
	cudaMallocManaged(&result,  2*sizeof(unsigned)+Ma*Nb*sizeof(double));
	
	copyMatrixToGPU(d_A, A);
	copyMatrixToGPU(d_B, B);
	
	printf("Przed:\n");
	printMatrixGPU("%.1f ", d_A);
	printf("\n");
	printMatrixGPU("%.1f ", d_B);
	printf("\n");


	*(unsigned*)result = *(unsigned*)d_A;
	*((unsigned*)result+1) = *((unsigned*)d_B+1); 
	mul<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(d_A, d_B, result);
	cudaDeviceSynchronize();
	
	printf("Po:\n");
	

	printMatrixGPU("%.1f ", d_A);
	printf("\n");
	printMatrixGPU("%.1f ", d_B);
	printf("\n");
	printMatrixGPU("%.2f ", result);
	printf("\n");




	cudaFree(d_A);
	cudaFree(d_B);
	cudaFree(result);*/

 // startMenu();
 //  const unsigned M = 3; //test values
 //  const unsigned N = 3;
 //  Matrix matrix1(M, N);
 //  Matrix matrix2(M, N);
 //
	// //test values
 //  for(unsigned i=0; i<M; i++){
 //  	for(unsigned j=0; j<N; j++){
 //  		matrix1.set(i, j, (double)2*i+j);
 //  	}
 //  }
 //  for(unsigned i=0; i<M; i++){
 //  	for(unsigned j=0; j<N; j++){
 //  		matrix2.set(i, j, (double)2*i+j);
 //  	}
 //  }
 //
	// 	size_t sizeOfMatrix1 = sizeof(matrix1) + matrix1.getRows() * matrix1.getColumns() * sizeof(double);
	// 	size_t sizeOfMatrix2 = sizeof(matrix2) + matrix2.getRows() * matrix2.getColumns() * sizeof(double);
	//  //host pointers
 //   void* h_ptr_matrixInt = malloc(sizeOfMatrix1);
 //   void* h_ptr_matrixInt2 = malloc(sizeOfMatrix2);
 //   void* h_ptr_product = malloc(sizeOfMatrix1);
 //
 //

	//srand(time(NULL));
	//reserve matrix of ints
	// const unsigned M = 3; //test values
	// const unsigned N = 3;
	// Matrix A(M,N);
	// for(unsigned i=0; i<M; i++){
	// 	for(unsigned j=0; j<N; j++){
	// 		A.set(i, j, (double)2*i+j);
	// 	}
	// }
	// printMatrixCPU("%.1f ", A);
	// Matrix B(M,N);
	// for(unsigned i=0; i<M; i++){
	// 	for(unsigned j=0; j<N; j++){
	// 		B.set(i, j, (double)2*i+j);
	// 	}
	// }
	// printMatrixCPU("%.1f ", B);
	// Matrix C = A+B;
	// printMatrixCPU("%.1f ", C);
	// Matrix D = A*B;
	// printMatrixCPU("%.1f ", D);
	// printf("%.1f \n", D.det());
	// printf("%.1f \n", A.det());
	// Matrix E(M,N);
	// E.set(0,0,1.0);
	// E.set(0,1,-1.0);
	// E.set(0,2,2.0);
	// E.set(1,0,3.0);
	// E.set(1,1,0.0);
	// E.set(1,2,-4.0);
	// E.set(2,0,2.0);
	// E.set(2,1,3.0);
	// E.set(2,2,5.0);
	/*for(unsigned i=0; i<M; i++){
		for(unsigned j=0; j<N; j++){
			E.set(i, j, (double)(j*j*j+i));
		}
	}*/
	// printMatrixCPU("%.1f ", E);
	// printf("%.1f \n", E.det());
	// Matrix F = E.inverse();
	// printMatrixCPU("%.8f ", F);


	// copyMatrixToGPU("%.1f ", d_ptr_matrixInt, matrixInt);
	// printMatrixGPU("%.1f ", d_ptr_matrixInt);
	// double* w = new double; (*w)=0.0;
	// double* partial_sums = new double[THREADS_PER_BLOCK*BLOCKS_PER_GRID];
	// for(unsigned i=0; i<THREADS_PER_BLOCK*BLOCKS_PER_GRID; i++) partial_sums[i] = 0.0;
	// det(d_ptr_matrixInt, partial_sums, w);
	// printf("%.1f\n", *w);
	// delete[] partial_sums;
	// delete w;
	// cudaFree(d_ptr_matrixInt);

	//reserve matrix of ints
	
	const unsigned M = 3; //test values
	const unsigned N = 3;
	Matrix matrix1(M,N);
 	Matrix matrix2(M,N);
	for(int i = 0 ;i<M;i++){
		for(int j = 0 ; j<N ; j++){
			matrix1.set(i,j,i);
		}
	}
	matrix1.set(0,0,1.0);
	matrix1.set(0,1,-1.0);
	matrix1.set(0,2,2.0);
	matrix1.set(1,0,3.0);
	matrix1.set(1,1,0.0);
	matrix1.set(1,2,-4.0);
	matrix1.set(2,0,2.0);
	matrix1.set(2,1,3.0);
	matrix1.set(2,2,5.0);

	for(int i = 0 ;i<M;i++){
		for(int j = 0 ; j<N ; j++){
			matrix2.set(i,j,j);
		}
	}

	//testowe wywolanie w mainie	
	/*	
	void* d_test = NULL;
	cudaMallocManaged(&d_test, 2*sizeof(unsigned)+matrix1.getRows()*matrix1.getColumns()*sizeof(double));
	copyMatrixToGPU(d_test, matrix1);
	double* tab_test = new double[BLOCKS_PER_GRID * THREADS_PER_BLOCK];
	for(unsigned i=0;i<BLOCKS_PER_GRID * THREADS_PER_BLOCK; i++) tab_test[i] = 0.0;
	double* w_test = new double; *w_test = 0.0;
	detHelper<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(d_test, tab_test);
	cudaDeviceSynchronize();
	for(unsigned i=0;i<BLOCKS_PER_GRID*THREADS_PER_BLOCK; i++){
		*w_test += tab_test[i];	
	}
	printf("%f", *w_test);

	cudaFree(d_test);
	delete[] tab_test;
	delete w_test;
	*/

	int x = 11;
	while (x != 0)
	{
		printf("Macierz 1:\n");
		printMatrixCPU("%.2f ", matrix1);
		printf("\n");
		printf("Macierz 2:\n");
		printMatrixCPU("%.2f ", matrix2);
		printf("\n");

		printf( "*******************MENU********************\n");
		printf("1.Wczytaj macierz z pliku\n");
		printf("2.Zapisz macierz do pliku\n");
		printf("3.Wyznacznik ( det(m1) )\n");
		printf("4.Dodaj macierze (m1=m1+m2)\n");
		printf("5.Przemnóż macierze (m1=m1*m2)\n");
		printf("6.Macierz odwrotna ( m1=m1^(-1) )\n");
    		printf("0.Wyjdź\n");

		printf("Wybieram : \n");
		scanf("%d", &x);
		system("clear");

    		int subselect= 0;
		switch (x){
		case 1: //ladowanie
			printf("1 - Macierz 1\n2 - Macierz 2\n");
			scanf("%d", &subselect);
			subselect==1?loadMatrixFromFile(matrix1):loadMatrixFromFile(matrix2);
		break;
		case 2: //zapisywanie
			printf("1 - Macierz 1\n2 - Macierz 2\n");
			scanf("%d", &subselect);
			subselect==1?saveMatrixToFile(matrix1):saveMatrixToFile(matrix2);
		break;
		case 3: //wyznacznik
				
			if(matrix1.getRows() != matrix1.getColumns()){
				printf("Wymiary macierzy musza byc rowne!!!");	
			}
			else{
				//deklarcja wskaznikow na pamiec urzadzenia
				void* d_1 = NULL;

				//poczatek mierzenia czasu wliczajac alokacja pamieci
				high_resolution_clock::time_point tzk = high_resolution_clock::now();

				//alokacja pamieci
				cudaMallocManaged(&d_1, 2*sizeof(unsigned)+matrix1.getRows()*matrix1.getColumns()*sizeof(double));

				//kopia macierzy w pamieci GPU:
				copyMatrixToGPU(d_1, matrix1);

				//dynamiczne przystosowanie liczby watkow (ominiecie bledu BLOCKS_PER_GRID * THREADS_PER_BLOCK >= N!)
				unsigned tempBLOCKS_PER_GRID = BLOCKS_PER_GRID;
				unsigned tempTHREADS_PER_BLOCK = THREADS_PER_BLOCK;

				switch(matrix1.getRows()){
					case 1:
					tempBLOCKS_PER_GRID=1;
					tempTHREADS_PER_BLOCK=1;
					break;
					case 2:
					tempBLOCKS_PER_GRID=1;
					tempTHREADS_PER_BLOCK=1;
					break;
					case 3:
					tempBLOCKS_PER_GRID=2;
					tempTHREADS_PER_BLOCK=2;
					break;
					case 4:
					tempBLOCKS_PER_GRID=4;
					tempTHREADS_PER_BLOCK=4;
					break;
					case 5:
					tempBLOCKS_PER_GRID=4;
					tempTHREADS_PER_BLOCK=16;
					break;
					case 6:
					tempBLOCKS_PER_GRID=8;
					tempTHREADS_PER_BLOCK=64;
					break;
					default:
					tempBLOCKS_PER_GRID = BLOCKS_PER_GRID;
					tempTHREADS_PER_BLOCK = THREADS_PER_BLOCK;
				}

				//tablica i wskaznik pomocnicze
				double* tab = NULL;//new double[BLOCKS_PER_GRID * THREADS_PER_BLOCK];
				cudaMallocManaged(&tab, sizeof(double) * tempBLOCKS_PER_GRID * tempTHREADS_PER_BLOCK);
				for(unsigned i=0;i<tempBLOCKS_PER_GRID * tempTHREADS_PER_BLOCK; i++) tab[i] = 0.0;
				double w = 0.0;
				
				//poczatek mierzenia czasu wliczajac alokacja pamieci
				high_resolution_clock::time_point tbk = high_resolution_clock::now();

				//jezeli rozmiar macierzy wiekszy od jeden wywolaj kernel, oblicz wyznacznik
				if(matrix1.getRows()>1){

					detHelper<<<tempBLOCKS_PER_GRID, tempTHREADS_PER_BLOCK>>>(d_1, tab);
					cudaDeviceSynchronize();
					for(unsigned i=0;i<tempBLOCKS_PER_GRID*tempTHREADS_PER_BLOCK; i++){
						w += tab[i];	
					}

				}
				//w przeciwnym wypadku wyznacznikiem jest jedyny element macierzy
				else w = *((double*)((unsigned*)d_1+2));

				//koniec mierzenia czasu wliczajac alokacja pamieci
				high_resolution_clock::time_point tk = high_resolution_clock::now();

				duration<double> time_alloc  = duration_cast<duration<double>>(tk - tzk);
				duration<double> time_nalloc = duration_cast<duration<double>>(tk - tbk);

				printf("(Device) Wyznacznik obliczony przy uzyciu kernela:\n");
				printf("%f\n\n", w);
				
				printf("(Device) wykonano w ciagu: %fs\n", time_nalloc.count());
				printf("(Device) ... %fs, jesli uzwglednic zarzadzanie pamiecia.\n", time_alloc.count());
				printf("\n");
				
				cudaFree(d_1);
				cudaFree(tab);
				//delete[] tab;

				//poczatek pomiaru czasu na CPU
				high_resolution_clock::time_point tpc = high_resolution_clock::now();

				double d = matrix1.det();
				
				//koniec pomiaru czasu na CPU
				high_resolution_clock::time_point tpk = high_resolution_clock::now();

				duration<double> timeCPU = duration_cast<duration<double>>(tpk - tpc);

				printf("Wynik obliczen na CPU:\n");
				printf("%f\n\n", d);
				printf("(Host) wykonano w ciagu: %fs\n", timeCPU.count());
				printf("\n");


			}

		break;
		case 4: //dodawanie
			if(matrix1.getRows() != matrix2.getRows() || matrix1.getColumns()!=matrix2.getColumns()){
				printf("Macierze różnych rozmiarów!\n");
			}
			else{
				//realizacja na GPU:
				//deklaracja wskaznikow na pamiec urzadzenia
				void* d_1 = NULL;
				void* d_2 = NULL;
				void* d_r = NULL;
				
				//poczatek mierzenia czasu wliczajac alokacja pamieci
				high_resolution_clock::time_point tzk = high_resolution_clock::now();

				//alokacja pamieci
				cudaMallocManaged(&d_1, 2*sizeof(unsigned)+matrix1.getRows()*matrix1.getColumns()*sizeof(double));
				cudaMallocManaged(&d_2, 2*sizeof(unsigned)+matrix1.getRows()*matrix1.getColumns()*sizeof(double));
				cudaMallocManaged(&d_r, 2*sizeof(unsigned)+matrix1.getRows()*matrix1.getColumns()*sizeof(double));
		
				//kopie macierzy w pamieci GPU:
				copyMatrixToGPU(d_1, matrix1);
				copyMatrixToGPU(d_2, matrix2);

				//poczatek mierzenia czasu, jedynie obliczenia
				high_resolution_clock::time_point tbk = high_resolution_clock::now();

				//ustawienie rozmiaru macierzy wynikowej
				*(unsigned*)d_r = matrix1.getRows();
				*((unsigned*)d_r+1) = matrix1.getColumns();
				add<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(d_1, d_2, d_r);	
				cudaDeviceSynchronize();

				//koniec mierzenia czasu
				high_resolution_clock::time_point tk = high_resolution_clock::now();

				//obliczenie delt czasowych:
				duration<double> time_alloc  = duration_cast<duration<double>>(tk - tzk);
				duration<double> time_nalloc = duration_cast<duration<double>>(tk - tbk);

				printf("Wynik obliczen na GPU:\n");
				printMatrixGPU("%.2f ", d_r);
				printf("\n");
				printf("(Device) wykonano w ciagu: %fs\n", time_nalloc.count());
				printf("(Device) ... %fs, jesli uzwglednic zarzadzanie pamiecia.\n", time_alloc.count());
				printf("\n");

				cudaFree(d_1);
				cudaFree(d_2);
				cudaFree(d_r);
				
				//poczatek pomiaru czasu na CPU
				high_resolution_clock::time_point tpc = high_resolution_clock::now();

				Matrix result = matrix1 + matrix2;//(matrix1.getRows(),matrix.getColumns());
				
				//koniec pomiaru czasu na CPU
				high_resolution_clock::time_point tpk = high_resolution_clock::now();

				duration<double> timeCPU = duration_cast<duration<double>>(tpk - tpc);
				
					
				
				printf("Wynik obliczen na CPU:\n");
				printMatrixCPU("%.2f ", result);
				printf("\n");
				printf("(Host) wykonano w ciagu: %fs\n", timeCPU.count());
				printf("\n");

				for(unsigned i=0;i<matrix1.getRows();i++){
					for(unsigned j=0;j<matrix1.getColumns();j++){
						matrix1.set(i,j,result.get(i,j));
					}
				}

			}
			//else executeAdding(matrix1, matrix2);
		break;
		case 5: //mnozenie
			if(matrix1.getColumns() != matrix2.getRows()){
				printf("Nie zgadzaja sie rozmiary macierzy!!!\n");
			}
			else{
				//realizacja na GPU:
				//deklaracja wskaznikow na pamiec urzadzenia
				void* d_1 = NULL;
				void* d_2 = NULL;
				void* d_r = NULL;
				
				//poczatek mierzenia czasu wliczajac alokacja pamieci
				high_resolution_clock::time_point tzk = high_resolution_clock::now();

				//alokacja pamieci
				cudaMallocManaged(&d_1, 2*sizeof(unsigned)+matrix1.getRows()*matrix1.getColumns()*sizeof(double));
				cudaMallocManaged(&d_2, 2*sizeof(unsigned)+matrix2.getRows()*matrix2.getColumns()*sizeof(double));
				cudaMallocManaged(&d_r, 2*sizeof(unsigned)+matrix1.getRows()*matrix2.getColumns()*sizeof(double));

				//kopie macierzy w pamieci GPU:
				copyMatrixToGPU(d_1, matrix1);
				copyMatrixToGPU(d_2, matrix2);

				//poczatek mierzenia czasu, jedynie obliczenia
				high_resolution_clock::time_point tbk = high_resolution_clock::now();

				//ustawienie rozmiaru macierzy wynikowej
				*(unsigned*)d_r = matrix1.getRows();
				*((unsigned*)d_r+1) = matrix2.getColumns();
				mul<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(d_1, d_2, d_r);	
				cudaDeviceSynchronize();

				//koniec mierzenia czasu
				high_resolution_clock::time_point tk = high_resolution_clock::now();

				//obliczenie delt czasowych:
				duration<double> time_alloc  = duration_cast<duration<double>>(tk - tzk);
				duration<double> time_nalloc = duration_cast<duration<double>>(tk - tbk);

				printf("Wynik obliczen na GPU:\n");
				printMatrixGPU("%.2f ", d_r);
				printf("\n");
				printf("(Device) wykonano w ciagu: %fs\n", time_nalloc.count());
				printf("(Device) ... %fs, jesli uzwglednic zarzadzanie pamiecia.\n", time_alloc.count());
				printf("\n");

				cudaFree(d_1);
				cudaFree(d_2);
				cudaFree(d_r);
				
				//poczatek pomiaru czasu na CPU
				high_resolution_clock::time_point tpc = high_resolution_clock::now();

				Matrix result = matrix1 * matrix2;//(matrix1.getRows(),matrix.getColumns());
				
				//koniec pomiaru czasu na CPU
				high_resolution_clock::time_point tpk = high_resolution_clock::now();

				duration<double> timeCPU = duration_cast<duration<double>>(tpk - tpc);
					
				printf("Wynik obliczen na CPU:\n");
				printMatrixCPU("%.2f ", result);
				printf("\n");
				printf("(Host) wykonano w ciagu: %fs\n", timeCPU.count());
				printf("\n");

				matrix1.setSize(result.getRows(),result.getColumns());

				for(unsigned i=0;i<matrix1.getRows();i++){
					for(unsigned j=0;j<matrix1.getColumns();j++){
						matrix1.set(i,j,result.get(i,j));
					}
				}



			}
		break;
    		case 6: //macierz odwrotna
			if(matrix1.getColumns() != matrix1.getRows()){
				printf("Nie zgadzaja sie rozmiary macierzy!!!\n");
			}
			else{
				//deklarcja wskaznikow na pamiec urzadzenia
				void* d_1 = NULL;
				void* d_r = NULL;

				//poczatek mierzenia czasu wliczajac alokacja pamieci
				high_resolution_clock::time_point tzk = high_resolution_clock::now();

				//alokacja pamieci
				cudaMallocManaged(&d_1, 2*sizeof(unsigned)+matrix1.getRows()*matrix1.getColumns()*sizeof(double));
				cudaMallocManaged(&d_r, 2*sizeof(unsigned)+matrix1.getRows()*matrix1.getColumns()*sizeof(double));

				//kopia macierzy w pamieci GPU:
				copyMatrixToGPU(d_1, matrix1);

				//dynamiczne przystosowanie liczby watkow (ominiecie bledu BLOCKS_PER_GRID * THREADS_PER_BLOCK >= N!)
				unsigned tempBLOCKS_PER_GRID = BLOCKS_PER_GRID;
				unsigned tempTHREADS_PER_BLOCK = THREADS_PER_BLOCK;

				switch(matrix1.getRows()){
					case 1:
					tempBLOCKS_PER_GRID=1;
					tempTHREADS_PER_BLOCK=1;
					break;
					case 2:
					tempBLOCKS_PER_GRID=1;
					tempTHREADS_PER_BLOCK=1;
					break;
					case 3:
					tempBLOCKS_PER_GRID=2;
					tempTHREADS_PER_BLOCK=2;
					break;
					case 4:
					tempBLOCKS_PER_GRID=4;
					tempTHREADS_PER_BLOCK=4;
					break;
					case 5:
					tempBLOCKS_PER_GRID=4;
					tempTHREADS_PER_BLOCK=16;
					break;
					case 6:
					tempBLOCKS_PER_GRID=8;
					tempTHREADS_PER_BLOCK=64;
					break;
					default:
					tempBLOCKS_PER_GRID = BLOCKS_PER_GRID;
					tempTHREADS_PER_BLOCK = THREADS_PER_BLOCK;
				}


				//tablica i wskaznik pomocnicze
				double* tab = NULL;//new double[BLOCKS_PER_GRID * THREADS_PER_BLOCK];
				cudaMallocManaged(&tab, sizeof(double) * tempBLOCKS_PER_GRID * tempTHREADS_PER_BLOCK);
				for(unsigned i=0;i<tempBLOCKS_PER_GRID * tempTHREADS_PER_BLOCK; i++) tab[i] = 0.0;
				double w = 0.0;

				//poczatek mierzenia czasu wliczajac alokacja pamieci
				//high_resolution_clock::time_point tbk = high_resolution_clock::now();

				if(matrix1.getColumns()>1){
					detHelper<<<tempBLOCKS_PER_GRID, tempTHREADS_PER_BLOCK>>>(d_1, tab);
					cudaDeviceSynchronize();
					for(unsigned i=0;i<tempBLOCKS_PER_GRID*tempTHREADS_PER_BLOCK; i++){
						w += tab[i];	
					}
				}
				//w przeciwnym wypadku wyznacznikiem jest jedyny element macierzy
				else w = *((double*)((unsigned*)d_1+2));


				double mainDet = w;

				//debug
				printf("Wyznacznik glownej macierzy: %f\n", mainDet);
				
				
				
				
				
				if(mainDet == 0.0 || mainDet == -0.0){
					printf("Wyznacznik macierzy rowny 0.0 - macierz odwrotna nie istnieje!");
				}
				else{
					*(unsigned*)d_r = matrix1.getRows();
					*((unsigned*)d_r+1) = matrix1.getColumns();
					if(matrix1.getRows()==1){
						*((double*)((unsigned*)d_r+2)) = 1.0/mainDet;
					}
					else if(matrix1.getRows()>1){
						//Tworzymy bufor na minora:
						void* minor = NULL;
						cudaMallocManaged(&minor, 2*sizeof(unsigned)+(matrix1.getRows()-1)*(matrix1.getColumns()-1)*sizeof(double));

						//przekopiowanie	
						*(unsigned*)minor = matrix1.getRows()-1;
						*((unsigned*)minor+1) = matrix1.getColumns()-1;

						//debug
						/*
						printf("Udalo sie tu dojsc!\n");
						printMatrixGPU("%f ",minor);
						*/

						for(unsigned i=0; i<matrix1.getRows(); i++){
							for(unsigned j=0; j<matrix1.getColumns(); j++){
								//zerujemy tab i w
								for(unsigned q=0;q<tempBLOCKS_PER_GRID * tempTHREADS_PER_BLOCK; q++) tab[q] = 0.0;
								w = 0.0;
								for(unsigned q=0;q<matrix1.getRows(); q++){
									for(unsigned k=0;k<matrix1.getRows(); k++){
										if(q<i && k<j) *((double*)((unsigned*)minor+2)+q*(matrix1.getColumns()-1)+k) = *((double*)((unsigned*)d_1+2) + q*matrix1.getColumns() + k);
										if(q<i && k>j) *((double*)((unsigned*)minor+2)+q*(matrix1.getColumns()-1)+k-1) = *((double*)((unsigned*)d_1+2) + q*matrix1.getColumns() + k);
										if(q>i && k<j) *((double*)((unsigned*)minor+2)+(q-1)*(matrix1.getColumns()-1)+k) = *((double*)((unsigned*)d_1+2) + q*matrix1.getColumns() + k);
										if(q>i && k>j) *((double*)((unsigned*)minor+2)+(q-1)*(matrix1.getColumns()-1)+k-1) = *((double*)((unsigned*)d_1+2) + q*matrix1.getColumns() + k);

									}	
								}	

								unsigned temp2BLOCKS_PER_GRID = BLOCKS_PER_GRID;
								unsigned temp2THREADS_PER_BLOCK = THREADS_PER_BLOCK;

								switch(matrix1.getRows()-1){
									case 1:
									temp2BLOCKS_PER_GRID=1;
									temp2THREADS_PER_BLOCK=1;
									break;
									case 2:
									temp2BLOCKS_PER_GRID=1;
									temp2THREADS_PER_BLOCK=1;
									break;
									case 3:
									temp2BLOCKS_PER_GRID=2;
									temp2THREADS_PER_BLOCK=2;
									break;
									case 4:
									temp2BLOCKS_PER_GRID=4;
									temp2THREADS_PER_BLOCK=4;
									break;
									case 5:
									temp2BLOCKS_PER_GRID=4;
									temp2THREADS_PER_BLOCK=16;
									break;
									case 6:
									temp2BLOCKS_PER_GRID=8;
									temp2THREADS_PER_BLOCK=64;
									break;
									default:
									temp2BLOCKS_PER_GRID = BLOCKS_PER_GRID;
									temp2THREADS_PER_BLOCK = THREADS_PER_BLOCK;
								}

								if(matrix1.getColumns()-1>1){

									//kernel
									detHelper<<<temp2BLOCKS_PER_GRID, temp2THREADS_PER_BLOCK>>>(minor, tab);
									cudaDeviceSynchronize();

									for(unsigned q=0;q<temp2BLOCKS_PER_GRID*temp2THREADS_PER_BLOCK; q++){
										w += tab[q];	
									}
								}

								else w = *((double*)((unsigned*)minor+2));

								
								*((double*)((unsigned*)d_r+2)+j*matrix1.getColumns()+i) = (((i%2)+(j%2))%2==0?1.0:-1.0) * w / mainDet;

							}
						}
						
						cudaFree(minor);

					}
					//koniec mierzenia czasu
					high_resolution_clock::time_point tk = high_resolution_clock::now();

					//obliczenie delt czasowych:
					duration<double> time_alloc  = duration_cast<duration<double>>(tk - tzk);
					//duration<double> time_nalloc = duration_cast<duration<double>>(tk - tbk);

					printf("Wynik obliczen na GPU:\n");
					printMatrixGPU("%.2f ", d_r);
					printf("\n");
					printf("(Device) wykonano w ciagu: %fs\n", time_alloc.count());
					//printf("(Device) ... %fs, jesli uzwglednic zarzadzanie pamiecia.\n", time_alloc.count());
					printf("\n");


				}
				

				cudaFree(tab);
				cudaFree(d_r);
				cudaFree(d_1);

				if(mainDet != 0.0 && mainDet != -0.0){

					try{
						//poczatek pomiaru czasu na CPU
						high_resolution_clock::time_point tpc = high_resolution_clock::now();
						
						Matrix result = matrix1.inverse();

						//koniec pomiaru czasu na CPU
						high_resolution_clock::time_point tpk = high_resolution_clock::now();

						duration<double> timeCPU = duration_cast<duration<double>>(tpk - tpc);	
						
						printf("Wynik obliczen na CPU:\n");
						printMatrixCPU("%.2f ", result);
						printf("\n");
						printf("(Host) wykonano w ciagu: %fs\n", timeCPU.count());
						printf("\n");

						for(unsigned i=0;i<matrix1.getRows();i++){
							for(unsigned j=0;j<matrix1.getColumns();j++){
								matrix1.set(i,j,result.get(i,j));
							}
						}
					}
					catch(...){}
				}
			}
      		break;
    		case 0: //wyjscie
      			exit(0);
    		default:
			printf("Nie rozumiem - sprobuj ponownie!");
      		break;

		}
	}
	
	/*
	Matrix a(2,5);
	Matrix b = a;
	printMatrixCPU("%.1f ",a);
	printMatrixCPU("%.1f ",b);
	*/
	return 0;
}
