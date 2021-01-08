#include <locale>
#include <stdlib.h>
#include <iostream>
#include <malloc.h>
#include "cmath"
#include<string>
#include <fstream>

using namespace std;

#define N_MIN	3			// минимальное количество вершин
#define N_MAX	30			// максимальное количество вершин

#define ALPHA	1			// вес фермента
#define BETA	3			// коэффициент эвристики

#define T_MAX	20			// время жизни колонии
#define M		50			// количество муравьев в колонии
#define Q		100			// количество феромона
#define RHO		0.2			// коэффициент испарения феромона
#define MIN_Q	20          //минимальное количество феромона
#define MAX_Q   100         //максимальное количество феромона

// структура ПУТЬ (длинна, массив вершин, количество вершин)
struct WAY_TYPE {
	int itabu;
	int length;
	int* tabu;
	int days;
};
//проверка выполнения предшествующих задач
bool all_vertex(int* passed_vertex, int* all_vertex, int itabu, int N)
{
	bool f1=true, f2 = true;
	for (int i = 0; i < N; ++i)
	{
		if (all_vertex[i] == 1)
		{
			f1 = false;
			//поиск задачи из списка необходимых в списке сделанных
			for (int j = 0; j < itabu; ++j)
			{
				if (passed_vertex[j] == i)
					f1 = true;
			}
		}
		if (!f1)
		{
			f2 = false;
			break;
		}
	}
	return f2;
}
// вероятность перехода муравья ant в вершину to
double probability(int to, WAY_TYPE ant, double** pheromone, double** distance, int vertex) {
	// если вершина уже посещена, возвращаем 0
	for (int i = 0; i < ant.itabu; ++i) if (to == ant.tabu[i]) return 0;

	double sum = 0.0;
	int from = ant.tabu[ant.itabu - 1];
	// считаем сумму в знаминателе
	for (int j = 0; j < vertex; ++j) {
		int flag = 1;
		// проверяем, посещал ли муравей j вершину
		for (int i = 0; i < ant.itabu; ++i) if (j == ant.tabu[i]) flag = 0;
		// если нет, тогда прибавляем к общей сумме
		if (flag) sum += pow(pheromone[from][j], ALPHA) * pow(distance[from][j], BETA);
	}
	// возвращаем значение вероятности
	return pow(pheromone[from][to], ALPHA) * pow(distance[from][to], BETA) / sum;
}

// основная функция алгоритма поиска
WAY_TYPE AntColonyOptimization(double** distance0, int vertex, int start, double h, int z, int**List) {
	// инициализация данных о лучшем маршруте
	WAY_TYPE way;
	way.itabu = 0;
	way.length = -1;
	way.tabu = (int*)malloc(sizeof(int) * vertex);
	way.days = -1;
	// инициализация данных о расстоянии и количестве феромона
	double** distance = NULL, ** pheromone = NULL;
	distance = (double**)malloc(sizeof(double*) * vertex);
	pheromone = (double**)malloc(sizeof(double*) * vertex);
	for (int i = 0; i < vertex; ++i) {
		distance[i] = (double*)malloc(sizeof(double) * vertex);
		pheromone[i] = (double*)malloc(sizeof(double) * vertex);
		for (int j = 0; j < vertex; ++j) {
			pheromone[i][j] = 200.0 / vertex;
			if (i != j) distance[i][j] = 1.0 / distance0[i][j];
		}
	}
	// инициализация муравьев
	WAY_TYPE ants[M];
	for (int k = 0; k < M; ++k) {
		ants[k].itabu = 0;
		ants[k].days = 0;
		ants[k].length = 0.0;
		ants[k].tabu = (int*)malloc(sizeof(int) * vertex);
		ants[k].tabu[ants[k].itabu++] = start;
	}
	int w = 0;
	int pr[100000];
	int hours;
	// основной цикл
	for (int t = 0; t < T_MAX; ++t) {
		// цикл по муравьям
		hours = z;
		for (int k = 0; k < M; ++k) {
			// поиск маршрута для k-го муравья

			do {
				int j_max = -1;
				double p_max = 0.0;
				for (int j = 0; j < vertex; ++j) {
					// Проверка вероятности перехода в вершину j
					if (ants[k].tabu[ants[k].itabu - 1] != j) {
						double p = probability(j, ants[k], pheromone, distance, vertex);
						if (rand() % 101 < p * 100 && p > 0 && all_vertex(ants[k].tabu,List[j], ants[k].itabu, vertex))
						{
							pr[w] = j;
							w++;
						}
					}
				}
				if (w) {
					j_max = pr[rand() % w];
					for (int i = 0; i < w + 1; i++)
					{
						pr[i] = 0;
					}
					w = 0;
					ants[k].length += distance0[ants[k].tabu[ants[k].itabu - 1]][j_max];
					ants[k].tabu[ants[k].itabu++] = j_max;


					if (hours + distance0[ants[k].tabu[ants[k].itabu - 1]][j_max] > h)
					{
						ants[k].days += 1;
						hours = distance0[ants[k].tabu[ants[k].itabu - 1]][j_max];
					}
					else
						hours += distance0[ants[k].tabu[ants[k].itabu - 1]][j_max];
				}
				if (ants[k].itabu == vertex && hours != 0)
					ants[k].days += 1;
			} while (ants[k].itabu != vertex);
			// проверка на лучшее решение
			if (ants[k].days < way.days || way.days < 0) {
				way.itabu = ants[k].itabu;
				way.days = ants[k].days;
				way.length = ants[k].length;
				for (int i = 0; i < way.itabu; ++i) way.tabu[i] = ants[k].tabu[i];
			}
			// обновление муравьев
			ants[k].itabu = 1;
			ants[k].length = 0.0;
			ants[k].days = 0;
		}
		// оставляем феромон на пути муравья
		for (int i = 0; i < way.itabu - 1; ++i) {
			int from = way.tabu[i % way.itabu];
			int to = way.tabu[(i + 1) % way.itabu];
			if (pheromone[from][to] + Q / way.length < MAX_Q)
			{
				pheromone[from][to] += Q / way.length;
				pheromone[to][from] = pheromone[from][to];
			}
		}
		// цикл по ребрам
		for (int i = 0; i < vertex; ++i)
		{
			for (int j = 0; j < vertex; ++j)
				// обновление феромона для ребра (i, j)
				if (i != j)
				{
					if ((pheromone[i][j] * (1 - RHO)) > MIN_Q)
					{
						pheromone[i][j] *= (1 - RHO);

					}
					else
						pheromone[i][j] = MIN_Q;

				}
		}



	}
	// возвращаем кратчайший маршрут
	return way;
}

// точка входа в программу
int main() {
	setlocale(LC_ALL, "Russian");

	double** D = NULL;
	char** L = NULL;
	int** List = NULL;
	double* F = NULL;
	int N = 0, A = 0, B = 0, h;
	double h1;
	// инициализация матрицы расстояний
	//////////////////////////////////////////////
	fstream infile;
	string a;
	cout << "Введите имя файла:\n";
	cin >> a;

	infile.open(a);
	if (infile.fail())
	{
		cout << "Ошибка открытия файла";
	}
	else {
		infile >> N;
	L = (char**)malloc(sizeof(char*) * N);
	F = (double*)malloc(sizeof(double) * N);
	for (int i = 0; i < N; ++i)
	{
		L[i] = (char*)malloc(sizeof(char) * 256);
		infile >> L[i];
		infile >> F[i];
	}	
	infile >>A >> h >> h1;//считывание точки начала, количества людей, которых можно задействовать, количество рабочих часов в день
	List = (int**)malloc(sizeof(int*) * N);

	for (int i = 0; i < N; ++i) {
		List[i] = (int*)malloc(sizeof(int) * N);
		for (int j = 0; j < N; ++j)
			infile >> List[i][j];
	}
	D = (double**)malloc(sizeof(double*) * N);

	for (int i = 0; i < N; ++i) {
		D[i] = (double*)malloc(sizeof(double) * N);
		for (int j = 0; j < N; ++j)
			D[i][j] = F[j];

	}
	// инициализация начальной точеки
	int z = F[A - 1];
	// запускаем алгоритм
	WAY_TYPE way = AntColonyOptimization(D, N, --A, h * h1, z,List);

	// выводим результат
	cout << "Время выполнения: " << way.days <<" дней" <<endl;
	
	int p = 1;
	cout << "Последовательность выполнения работ:\nДень 1 :\n" ;
	int hours = 0;
	for (int i = 0; i < way.itabu; ++i) 
	{
		
		
		hours += F[++way.tabu[i] - 1];	
		if (hours > h1*h) 
		{
			cout << "День: " <<++p <<'\n';
			hours = F[way.tabu[i] - 1];
		}
		
		cout << L[way.tabu[i] - 1]<<endl;
		
	}
	}
	free(L);
	free(F);
	for (int i = 0; i < N; ++i) {
		free(D[i]);
		free(List[i]);
		
	}
	return 0;
}
