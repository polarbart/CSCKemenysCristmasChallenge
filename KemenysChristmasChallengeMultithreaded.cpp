// CSC.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <intrin.h>
#include <iostream>
#include <algorithm>
#include <ctime>
#include <thread>
#include <vector>
#include <atomic>

#define N 11
#define SOL 20
#define INTCEIL(x, y) (((x) + (y) - 1) / (y))

#define SET(i, j) do {mb[(i)] |= (1ui16 << (j)); mb[(j)] &= ~(1ui16 << (i));} while(false) // "1ui16" may only work in visual studio. Try to use 1 if not working
#define CLEAR(i, j) do {mb[(i)] &= ~(1ui16 << (j)); mb[(j)] |= (1ui16 << (i));} while(false)
#define NOE(i) (__popcnt16(mb[(i)])) // popcount may also only work with visual studio
#define NOEU(i, j) (__popcnt16(mb[(i)] & ((1ui16 << (j)) - 1)))
#define GETB(i, j) (mb[(i)] & (1ui16 << (j))) // 0 if false else != 0
#define GETI(i, j) ((mb[(i)] & (1ui16 << (j))) >> (j)) // 0 if false else 1


std::atomic_int64_t calls = 0;


bool swap_it_like_its_hot(const uint16_t mb[], int idc[]) {
	for (int i = N - 1; i >= 1; i--)
	{
		int v = 0;
		for (int k = i - 1; k >= 0; k--)
		{
			v += GETB(idc[i], idc[k]) ? -1 : 1; //make sense?
			if (v > 0) {
				int t = idc[i];
				for (int x = i; x >= k + 1; x--)
				{
					idc[x] = idc[x - 1];
				}
				idc[k] = t;
				return true;
			}
		}
	}
	return false;
}

void check_one(const uint16_t mb[], int idc[]) {
	calls++;
	int aks = 0;
	for (int k = 1; k < N; k++)
		for (int l = 0; l < k; l++)
			aks += GETI(idc[l], idc[k]);
	while (aks > SOL && swap_it_like_its_hot(mb, idc)) {
		aks--;
	}

	int cnt = 0;
	while (aks > SOL && cnt < 10000) {
		std::random_shuffle(idc, idc + N);
		aks = 0;
		for (int k = 1; k < N; k++)
			for (int l = 0; l < k; l++)
				aks += GETI(idc[l], idc[k]);
		while (aks > SOL && swap_it_like_its_hot(mb, idc)) {
			aks--;
		}
		cnt++;
	}
	if (cnt > 5000) {
		std::cout << "Error!" << std::endl;
	}
}

void setup_next_col_perm(const int i, const int max, uint16_t mb[], int indices[]) {
	for (int j = i + 1; j < N; j++) {
		if ((j - i) <= max)
			SET(i, j);
		else
			CLEAR(i, j);
	}
	indices[i] = i + 1;
}

bool next_col_perm(const int i, const int min, uint16_t mb[], int indices[]) {
	if (!GETB(i, indices[i])) {
		return false;
	}
	int k = 0;
	while (indices[i] + k + 1 < N && GETB(i, indices[i] + k + 1)) {
		CLEAR(i, indices[i] + k);
		SET(i, k + i + 1);
		k++;
	}
	if (indices[i] + k == N - 1) {
		if (k < min)
			return false;
		CLEAR(i, indices[i] + k);
		indices[i] = i + 1;
		return true;
	}
	CLEAR(i, indices[i] + k);
	SET(i, indices[i] + k + 1);
	if (k > 0)
		indices[i] = i + 1;
	else
		indices[i] += 1;
	return true;
}

int get_min_noe(const int i, const uint16_t mb[]) {
	int max = 0;
	int sum = NOEU(i, i);
	for (int j = i + 1; j < N; j++) {
		sum += NOEU(j, i);
		if (NOEU(j, i) > max)
			max = NOEU(j, i);
	}
	int ret = INTCEIL((N - i) * (N - i - 1) / 2 + sum, N - i);
	ret = (max > ret ? max : ret) - NOEU(i, i);
	return ret < 0 ? 0 : ret;
}

int get_max_noe(const int i, const uint16_t mb[]) {
	int a = NOE(i - 1) - NOEU(i, i) - GETI(i, i - 1);
	int b = N - i - 1;
	return (a < b ? a : b);
}

void check_graphs(const uint16_t start) {
	uint16_t mb[N] = {};
	int idc[N];
	int indices[N] = {1};
	int min_noe[N] = {N/2};

	for (int k = 0; k < N; k++)
		idc[k] = k;

	mb[0] = start;
	for (int i = 1; i < N; i++)
		if (!GETB(0, i))
			mb[i] = (1ui16);
	
	int j = 0;
	do {
		for (int i = j + 1; i < N - 1; i++) {
			min_noe[i] = get_min_noe(i, mb);
			int max_noe = get_max_noe(i, mb);
			if (min_noe[i] > max_noe) {
				j = i - 1;
				goto abc;
			}
			setup_next_col_perm(i, max_noe, mb, indices);
		}

		if ((NOE(N - 2) - GETI(N - 1, N - 2)) >= NOE(N - 1))
			check_one(mb, idc);

		j = N - 2;
	abc:
		while (j >= 0 && !next_col_perm(j, min_noe[j], mb, indices))
			j--;
	} while (j > 0);
}

void check_all_graphs() {
	uint16_t mb[N] = {};
	int indices[N] = {};
	setup_next_col_perm(0, N - 1, mb, indices);
	std::vector<std::thread*> v;
	do v.push_back(new std::thread(check_graphs, mb[0]));
	while (next_col_perm(0, N / 2, mb, indices));
	for (auto t = v.begin(); t != v.end(); ++t)
		(*t)->join();
		
}

int main() {
	clock_t begin = std::clock();
	check_all_graphs();
	std::cout << double(std::clock() - begin) / CLOCKS_PER_SEC << "s" << std::endl;
	std::cout << calls << std::endl << std::endl;
}