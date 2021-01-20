#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

///////////////////////////////////////////////////////////////////////////////

// Define squared difference macro
#define sqdiff(a, b) (int) pow((double)a - (double)b, 2);

// DECLARATIONS

typedef struct
{
	size_t nRows, nCols;
	int** value;
} Matrix;

// Affinity Propagation declarations

#define MAX_ITER 200
#define DAMPING_FACTOR 0.5f

Matrix* sim; // similarity matrix
Matrix* res; // responsibility matrix
Matrix* ava; // availability matrix
Matrix* cri; // criterion matrix

void computeSimilarity(Matrix* data);
void computeResponsibility();
void computeAvailability();
void computeCriterion();

void damping(Matrix* new, Matrix* old);

void affinityPropagation(Matrix* data);

void affinityPropagationDebug();

///////////////////////////////////////////////////////////////////////////////

// Utilities

Matrix* createMatrix(size_t nRows, size_t nCols);

Matrix* copyMatrix(Matrix* src);

void deleteMatrix(Matrix* m);

int equalsMatrix(Matrix* a, Matrix* b);

void displayMatrix(FILE* out, Matrix* m);

void error(char* msg);

///////////////////////////////////////////////////////////////////////////////

// ENTRY POINT

int main(int argc, char* argv[])
{
	affinityPropagationDebug();
	return 0;
}

///////////////////////////////////////////////////////////////////////////////

// IMPLEMENTATIONS

void computeSimilarity(Matrix* data)
{
	
	// Similarity of x and y = negative euclidean distance between x and y
	for (size_t i = 0; i < sim->nRows; i++)
	{
		for (size_t j = i + 1; j < sim->nCols; j++)
		{
			for (size_t k = 0; k < data->nCols; k++)
			{
				int dif = sqdiff(data->value[i][k], data->value[j][k]);
				sim->value[i][j] -= dif;
				sim->value[j][i] -= dif;
			}
		}
	}

	// Find minimum in matrix 
	// (distance between the two most distant points from each other)
	int m = INT_MAX;
	for (size_t i = 0; i < sim->nRows; i++)
	{
		for (size_t j = 0; j < sim->nCols; j++)
		{
			m = min(sim->value[i][j], m);
		}
	}

	// Set self-similarity to the minimal value found beforehand
	// (So any point would never choose itself as an exemplar)
	for (size_t i = 0; i < sim->nRows; i++)
	{
		sim->value[i][i] = m;
	}
}

void computeResponsibility()
{
	for (size_t i = 0; i < res->nRows; i++)
	{
		for (size_t j = 0; j < res->nCols; j++)
		{
			/*
			* The responsibility of B to A = r(A,B) is the similarity of A and
			* B minus the maximum between, for each K which is not B, of the
			* similarity between A and K + the availability of K to A.
			*/
			res->value[i][j] = sim->value[i][j];
			int m = INT_MIN;
			for (size_t k = 0; k < res->nCols; k++)
			{
				if (j != k)
				{
					m = max(sim->value[i][k] + ava->value[i][k], m);
				}
			}
			res->value[i][j] -= m;
		}
	}
}

void computeAvailability()
{
	for (size_t i = 0; i < ava->nRows; i++)
	{
		for (size_t j = 0; j < ava->nCols; j++)
		{
			if (i == j)
			{
				/*
				* The self availability of A is the sum, for each K which is
				* not A, of the responsibility of A to K only if it is greater
				* than zero.
				*/
				ava->value[i][i] = 0;
				for (size_t k = 0; k < ava->nRows; k++)
				{
					if (i != k)
					{
						ava->value[i][i] += max(0, res->value[k][i]);
					}
				}
			}
			else
			{
				/*
				* The availability of B to A = a(A,B) is the self
				* responsibility of B + the sum, for each K which is not A
				* nor B, of the responsibility of B to K only if it is greater
				* than 0. If the computed availability is greater than zero, it
				* is reduced to zero.
				*/
				ava->value[i][j] = res->value[j][j];
				for (size_t k = 0; k < ava->nRows; k++)
				{
					if (i != k && j != k)
					{
						ava->value[i][j] += max(0, res->value[k][j]);
					}
				}
				ava->value[i][j] = min(0, ava->value[i][j]);
			}
		}
	}
}

void computeCriterion()
{
	for (size_t i = 0; i < cri->nRows; i++)
	{
		for (size_t j = 0; j < cri->nCols; j++)
		{
			cri->value[i][j] = res->value[i][j] + ava->value[i][j];
		}
	}
}

void damping(Matrix* new, Matrix* old)
{
	// Damping keeps the new values relative to the old ones to prevent 
	// numerical oscillations between iterations
	for (size_t i = 0; i < new->nRows; i++)
	{
		for (size_t j = 0; j < new->nCols; j++)
		{
			new->value[i][j] =
				(int)(new->value[i][j] * DAMPING_FACTOR) +
				(int)(old->value[i][j] * (1 - DAMPING_FACTOR));
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void affinityPropagation(Matrix* data)
{
	// Initialize all matrices.
	sim = createMatrix(data->nRows, data->nRows);
	res = createMatrix(data->nRows, data->nRows);
	ava = createMatrix(data->nRows, data->nRows);
	cri = createMatrix(data->nRows, data->nRows);

	Matrix* oldRes = createMatrix(data->nRows, data->nRows);
	Matrix* oldAva = createMatrix(data->nRows, data->nRows);

	computeSimilarity(data);

	int nIter = 0;
	int isStable = 0;
	do
	{
		nIter++;

		oldRes = copyMatrix(res);
		computeResponsibility();
		damping(res, oldRes);

		oldAva = copyMatrix(ava);
		computeAvailability();
		damping(ava, oldAva);

		isStable = equalsMatrix(res, oldRes) && equalsMatrix(ava, oldAva);

		deleteMatrix(oldRes);
		deleteMatrix(oldAva);

	} while (!isStable && nIter < MAX_ITER);

	deleteMatrix(oldRes);
	deleteMatrix(oldAva);

	computeCriterion();

	deleteMatrix(sim);
	deleteMatrix(res);
	deleteMatrix(ava);
}

void affinityPropagationDebug()
{
	Matrix* data = createMatrix(5, 5);
	int staticData[5][5] =
	{
		{3, 4, 3, 2, 1},
		{4, 3, 5, 1, 1},
		{3, 5, 3, 3, 3},
		{2, 1, 3, 3, 2},
		{1, 1, 3, 2, 3}
	};

	for (size_t i = 0; i < 5; i++)
	{
		memcpy(data->value[i], staticData[i], 5 * sizeof(int));
	}
	
	printf("Data :\n");
	displayMatrix(stdout, data);

	// Initialize all matrices.
	sim = createMatrix(data->nRows, data->nRows);
	res = createMatrix(data->nRows, data->nRows);
	ava = createMatrix(data->nRows, data->nRows);
	cri = createMatrix(data->nRows, data->nRows);

	Matrix* oldRes;
	Matrix* oldAva;

	computeSimilarity(data);

	printf("\nSimilarity :\n");
	displayMatrix(stdout, sim);

	int isStable = 0;
	do
	{
		system("PAUSE");

		oldRes = copyMatrix(res);
		computeResponsibility();
		damping(res, oldRes);

		printf("\nResponsibility :\n");
		displayMatrix(stdout, res);

		oldAva = copyMatrix(ava);
		computeAvailability();
		damping(ava, oldAva);

		printf("\nAvailability :\n");
		displayMatrix(stdout, ava);

		isStable = equalsMatrix(res, oldRes) && equalsMatrix(ava, oldAva);

		deleteMatrix(oldRes);
		deleteMatrix(oldAva);

	} while (!isStable);

	computeCriterion();

	deleteMatrix(sim);
	deleteMatrix(res);
	deleteMatrix(ava);

	printf("\nCriterion :\n");
	displayMatrix(stdout, cri);
}

///////////////////////////////////////////////////////////////////////////////

Matrix* createMatrix(size_t nRows, size_t nCols)
{
	Matrix* m = (Matrix*)malloc(sizeof(Matrix));
	if (m == NULL)
	{
		error("Failed to create Matrix.");
	}
	else
	{
		m->nRows = nRows;
		m->nCols = nCols;

		// calloc fills array with zeroes (unlike malloc)
		m->value = calloc(nRows, sizeof(int*));

		if (m->value == NULL)
		{
			error("Failed to create Matrix.value");
		}
		else
		{
			for (size_t i = 0; i < nRows; i++)
			{
				m->value[i] = calloc(nCols, sizeof(int));
				if (m->value[i] == NULL)
				{
					error("Failed to create Matrix.value col");
				}
			}
		}

	}
	return m;
}

Matrix* copyMatrix(Matrix* src)
{
	Matrix* dst = (Matrix*)malloc(sizeof(Matrix));

	if (dst == NULL)
	{
		error("Failed to create Matrix");
	}
	else
	{
		memcpy(dst, src, sizeof(Matrix));

		dst->value = calloc(dst->nRows, sizeof(int*));

		if (dst->value == NULL)
		{
			error("Failed to copy Matrix.value row");
		}
		else
		{
			for (size_t i = 0; i < dst->nRows; i++)
			{
				dst->value[i] = calloc(dst->nCols, sizeof(int));
				if (dst->value[i] == NULL)
				{
					error("Failed to copy Matrix.value col");
				}
				else
				{
					memcpy(dst->value[i], src->value[i], dst->nCols * sizeof(int));
				}
			}
		}
	}

	return dst;
}

void deleteMatrix(Matrix* m)
{
	free(m->value);
	free(m);
}

int equalsMatrix(Matrix* a, Matrix* b)
{
	if (a->nRows != b->nRows || a->nCols != b->nCols)
	{
		return 0; // false
	}
	
	for (size_t i = 0; i < a->nRows; i++)
	{
		for (size_t j = 0; j < a->nCols; j++)
		{
			if (a->value[i][j] != b->value[i][j])
			{
				return 0; // false
			}
		}
	}

	return 1; // true
}

void displayMatrix(FILE* out, Matrix* m)
{
	for (size_t i = 0; i < m->nRows; i++)
	{
		for (size_t j = 0; j < m->nCols; j++)
		{
			fprintf(out, "%d ", m->value[i][j]);
		}
		fprintf(out, "\n");
	}
}

void error(char* msg)
{
	fprintf(stderr, msg);
}