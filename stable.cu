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
	//printf("%d\n",*(unsigned*)matrix);
	//printf("%d\n",*((unsigned*)matrix+1));
	
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

			tab_helper[idx] += product;

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




using namespace std::chrono;

int main(){
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
	
	return 0;
}
