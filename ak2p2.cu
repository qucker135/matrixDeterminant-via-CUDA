#include <stdio.h> //printf
#include <stdlib.h>//srand, rand
#include <time.h>  //time

const unsigned BLOCKS_PER_GRID = 1;
const unsigned THREADS_PER_BLOCK = 2;

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
		for(unsigned i=0;i<M;i++){
			for(unsigned j=0;j<N;j++){
				double detElement = (((i%2)+(j%2))%2==0?1.0:-1.0)*minor(i,j).det();

				result.set(j,i,detElement/mainDet);
			}	
		}
		return result;
	}
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
	for(unsigned i=0; i<*(unsigned*)matrix; i++){  //wartosc *(unsigned*)matrix przechowuje liczbe wierszy
		for(unsigned j=0; j<*((unsigned*)matrix+1); j++){ //wartosc *(unsigned*)matrix+1 przechowuje liczbe kolumn 
			//duzo castowania
			printf(format, *((double*)((unsigned*)matrix+2)+i*(*(unsigned*)matrix)+j));
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
			*((double*)((unsigned*)ptr+2)+(i*matrix.getRows()+j)) = matrix.get(i,j);
		}
	}
	/*
	for(unsigned i=0; i<*(unsigned*)ptr; i++){  //wartosc *(unsigned*)matrix przechowuje liczbe wierszy
		for(unsigned j=0; j<*((unsigned*)ptr+1); j++){ //wartosc *(unsigned*)matrix+1 przechowuje liczbe kolumn 
			//duzo castowania
			printf(format, *((T*)((unsigned*)ptr+2)+i*(*(unsigned*)ptr)+j));
		}
		printf("\n");
	}
	*/
	//printMatrixGPU<T>(format, ptr);

}

__global__
void detHelper(void* matrix, double* tab_helper, unsigned N){
	unsigned idx = threadIdx.x + blockDim.x * blockIdx.x;
	unsigned stride = blockDim.x * gridDim.x;
	
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

			//printf("%d\n",tab_helper[idx]);

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


int main(){
	//srand(time(NULL));
	//reserve matrix of ints
	const unsigned M = 3; //test values
	const unsigned N = 3;
	Matrix A(M,N);
	for(unsigned i=0; i<M; i++){
		for(unsigned j=0; j<N; j++){
			A.set(i, j, (double)2*i+j);				
		}
	}
	printMatrixCPU("%.1f ", A);
	Matrix B(M,N);
	for(unsigned i=0; i<M; i++){
		for(unsigned j=0; j<N; j++){
			B.set(i, j, (double)2*i+j);				
		}
	}
	printMatrixCPU("%.1f ", B);
	Matrix C = A+B;
	printMatrixCPU("%.1f ", C);
	Matrix D = A*B;
	printMatrixCPU("%.1f ", D);
	printf("%.1f \n", D.det());
	printf("%.1f \n", A.det());
	Matrix E(M,N);
	E.set(0,0,1.0);	
	E.set(0,1,-1.0);	
	E.set(0,2,2.0);	
	E.set(1,0,3.0);	
	E.set(1,1,0.0);	
	E.set(1,2,-4.0);	
	E.set(2,0,2.0);	
	E.set(2,1,3.0);	
	E.set(2,2,5.0);	
	/*for(unsigned i=0; i<M; i++){
		for(unsigned j=0; j<N; j++){
			E.set(i, j, (double)(j*j*j+i));				
		}
	}*/
	printMatrixCPU("%.1f ", E);
	printf("%.1f \n", E.det());
	Matrix F = E.inverse();
	printMatrixCPU("%.8f ", F);
	//Kopia macierzy w pamięci GPU
	/*
	void* d_ptr_matrixInt = NULL;
		
	cudaMallocManaged(&d_ptr_matrixInt, sizeof(matrixInt));
	
	
	//Serializacja	
	*((unsigned*)d_ptr_matrixInt) = matrixInt.getRows();
	*((unsigned*)d_ptr_matrixInt+1) = matrixInt.getColumns();
	for(unsigned i=0;i<matrixInt.getRows();i++){
		for(unsigned j=0;j<matrixInt.getColumns();j++){
			*((int*)((unsigned*)d_ptr_matrixInt+2)+(i*matrixInt.getRows()+j)) = matrixInt.get(i,j);
		}
	}
	

	copyMatrixToGPU("%.1f ", d_ptr_matrixInt, matrixInt);
	printMatrixGPU("%.1f ", d_ptr_matrixInt);
	double* w = new double; (*w)=0.0; 
	double* partial_sums = new double[THREADS_PER_BLOCK*BLOCKS_PER_GRID];
	for(unsigned i=0; i<THREADS_PER_BLOCK*BLOCKS_PER_GRID; i++) partial_sums[i] = 0.0;
	det(d_ptr_matrixInt, partial_sums, w);
	printf("%.1f\n", *w);
	delete[] partial_sums;
	delete w;
	cudaFree(d_ptr_matrixInt);
	*/
	return 0;
}
