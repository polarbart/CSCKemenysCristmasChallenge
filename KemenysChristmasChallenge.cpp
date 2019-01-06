#include <intrin.h>
#include <iostream>
#include <algorithm>
#include <ctime>

// number of nodes
#define N 12 
// number of edges that are allowed to be removed
#define SOL 22

// ceil(x/y)
#define INTCEIL(x, y) (((x) + (y) - 1) / (y)) 

// m_{i, j} = 1 and m_{j, i} = 0
#define SET(i, j) do {mb[(i)] |= (1ui16 << (j)); mb[(j)] &= ~(1ui16 << (i));} while(false) // "1ui16" may only work with visual studio. Try to use 1_u or ((uint_16) 1) instead

// m_{i, j} = 0 and m_{j, i} = 1
#define CLEAR(i, j) do {mb[(i)] &= ~(1ui16 << (j)); mb[(j)] |= (1ui16 << (i));} while(false)

// c^N(i) (noe == Number Outgoing Edges)
#define NOE(i) (__popcnt16(mb[(i)]))  // __popcnt16 may also only work with visual studio

// c^j(i)
#define NOEU(i, j) (__popcnt16(mb[(i)] & ((1ui16 << (j)) - 1)))

// m_{i, j}
#define GET(i, j) ((mb[(i)] >> (j)) & 1ui16)

// the adjacency martix M
uint16_t mb[N];

// current ranking of vertices
int idc[N];

// index of the current leftmost one for the respective column
int indices[N];

// number of graphs we've checked
int calls = 0;

//improves the current ranking by moving one vertex up the ranking such that the anti-kemeny-score decreases by 1, returns false if not possible
bool swap_idcs() {
	for (int i = N - 1; i >= 1; i--) { //over all vertices
		int v = 0;
		for (int k = i - 1; k >= 0; k--) { //over all higher ranked vertices
			v += GET(idc[i], idc[k]) ? -1 : 1; //add up total change of aks
			if (v > 0) { // as soon as v = 1, perform the swap and return
				int t = idc[i];
				for (int x = i; x >= k + 1; x--)
					idc[x] = idc[x - 1];
				idc[k] = t;
				return true;
			}
		}
	}
	return false;
}

//checks whether the current graph has a permutation with less than "SOL" backwards edges using the described algorithm
void check_one() {
	calls++;
	for (int cnt = 0; cnt < 5000; cnt++) { //arbitrary number to terminate in case of impossibility (jumps to error print then)
	//compute the anti-kemeny-score
		int aks = 0;
		for (int k = 1; k < N; k++)
			for (int l = 0; l < k; l++)
				aks += GET(idc[l], idc[k]);
		//perform iterative improvements to the ranking as long as necessary and possible
		while (aks > SOL && swap_idcs())
			aks--; //the aks always decreases by exactly one
		if (aks <= SOL)
			return;
		//if aks still > SOL do random shuffle and try again
		std::random_shuffle(idc, idc + N);
	}
	std::cout << "Error!" << std::endl;
}

// initialize column i of lower triangular martix with max ones e.g. setup_next_col_perm(3, 4) -> m[3] = ????111100... 
void setup_next_col_perm(const int i, const int max) {
	for (int j = i + 1; j - i <= max; j++)
		SET(i, j);
	for (int j = i + max + 1; j < N; j++)
		CLEAR(i, j);
	indices[i] = i + 1;
}

/*
This function generates the next "permutation" if there is any i.e. if the number of ones >= min
it returns true if it generated a valid "permutation" else false

it generates the next "permutation" in the following fashion:
if the element right from the leftmost one is a zero then just shift the leftmost one once to the right
   e.g. 0100110 -> 0010110
else shift the one from the right once to the right and shift the leftmost one to the leftmost position
   e.g. 0011010 -> 1000110
   or   0011100 -> 1100010 -> 1010010 -> 0110010 -> 1001010
also if all ones are on the right and the number of ones is > min shift all ones to the leftmost posision and clear the rightmost one
   e.g. 000111 -> 110000

one full example with min = 2 and max = 4
   111100 -> 111010 -> 110110 -> 101110 -> 011110 -> 111001 -> 110101 
-> 101101 -> 011101 -> 110011 -> 101011 -> 011011 -> 100111 -> 010111 
-> 001111 -> 111000 -> 110100 -> 101100 -> 011100 -> 110010 -> 101010 
-> 011010 -> 100110 -> 010110 -> 001110 -> 110001 -> 101001 -> 011001 
-> 100101 -> 010101 -> 001101 -> 100011 -> 010011 -> 001011 -> 000111 
-> 110000 -> 101000 -> 011000 -> 100100 -> 010100 -> 001100 -> 100010 
-> 010010 -> 001010 -> 000110 -> 100001 -> 010001 -> 001001 -> 000101 
-> 000011
*/
bool next_col_perm(const int i, const int min) {
	// if there are no ones in the column there isn't a next "permutation"
	if (!GET(i, indices[i]))
		return false;

	int k = 0;
	// while there is a one to the right of the current one shift the current one to the leftmost position 
	while (indices[i] + k + 1 < N && GET(i, indices[i] + k + 1)) {
		CLEAR(i, indices[i] + k);
		SET(i, k + i + 1);
		k++;
	}

	// if all ones were on the right
	if (indices[i] + k == N - 1) {
		// if the number of ones == min there there isn't a next "permutation"
		if (k < min)
			return false;
		// else clear the rightmost one. The other ones have already been shifted to the left
		CLEAR(i, indices[i] + k);
		indices[i] = i + 1;
		return true;
	}

	// else shift the current one to the right
	CLEAR(i, indices[i] + k);
	SET(i, indices[i] + k + 1);

	// set indices[i] to the index of the leftmost one
	if (k > 0)
		indices[i] = i + 1;
	else
		indices[i] += 1;
	return true;
}

// compute the minimum number of outgoing edges i.e. min_i
int get_min_noe(const int i) {
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

// compute the maximum number of outgoing edges i.e. max_i
int get_max_noe(const int i) {
	int a = NOE(i - 1) - NOEU(i, i) - GET(i, i - 1);
	int b = N - i - 1;
	return (a < b ? a : b);
}

// check all graphs
void check_all() {
	int min_noe[N];

	min_noe[0] = N / 2;
	setup_next_col_perm(0, N - 1);

	for (int k = 0; k < N; k++)
		idc[k] = k;

	int j = 0; // current column
	while (j >= 0) {

		// initialize all columns that are uninitialized
		for (int i = j + 1; i < N - 1; i++) {
			min_noe[i] = get_min_noe(i);
			int max_noe = get_max_noe(i);
			// if min_k > max_k there are no "permutations"
			if (min_noe[i] > max_noe) {
				j = i - 1;
				goto abc;
			}
			setup_next_col_perm(i, max_noe);
		}

		// if the last column has <= ones than the second last column i.e. the graph is valid, check it 
		if ((NOE(N - 2) - GET(N - 1, N - 2)) >= NOE(N - 1))
			check_one();

		j = N - 2;
	abc:
		// generate the next valid graph
		while (j >= 0 && !next_col_perm(j, min_noe[j]))
			j--;
	}
}

int main() {
	clock_t begin = std::clock();

	check_all();

	std::cout << double(std::clock() - begin) / CLOCKS_PER_SEC << "s" << std::endl;
	std::cout << calls << std::endl;
}