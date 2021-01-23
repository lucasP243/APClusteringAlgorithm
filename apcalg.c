#pragma warning(disable:4996)

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

#define MAX_ITER INT_MAX
#define DAMPING_FACTOR 0.5f

Matrix* sim; // similarity matrix
Matrix* res; // responsibility matrix
Matrix* ava; // availability matrix
Matrix* cri; // criterion matrix

void damping(Matrix* new, Matrix* old);

void computeSimilarity(Matrix* data);
void computeResponsibility();
void computeAvailability();
void computeCriterion();

Matrix* extractExemplars();

Matrix* affinityPropagation(Matrix* data);

void affinityPropagationDebug();

///////////////////////////////////////////////////////////////////////////////

// Utilities

Matrix* createMatrix(size_t nRows, size_t nCols);

Matrix* copyMatrix(Matrix* src);

void deleteMatrix(Matrix* m);

int equalsMatrix(Matrix* a, Matrix* b);

void displayMatrix(FILE* out, Matrix* m);

void error(char* msg);

struct size_file findSize(FILE* fdata);

Matrix* loadData(char* path);

///////////////////////////////////////////////////////////////////////////////

// ENTRY POINT

int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		error("Usage : >> apcalg path_to_data_file");
		return -1;
	}
	else
	{
		char* path = argv[1];
		if (strcmp(path, "debug_test") == 0)
		{
			affinityPropagationDebug();
		}
		else
		{
			Matrix* data = loadData(path);
			Matrix* result = affinityPropagation(data);

			FILE* out = fopen("result.txt", "w");
			if (out == NULL)
			{
				error("Failed to create file.");
			}
			else
			{
				displayMatrix(out, result);
			}

			deleteMatrix(data);
			deleteMatrix(result);
		}
	}
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
			// C(i,j) = R(i,j) + A(i,j)
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

Matrix* extractExemplars()
{
	// For each individual (row), the column with the highest value defines
	// this individual's best choice for exemplar
	Matrix* exe = createMatrix(cri->nRows, cri->nCols);
	for (size_t i = 0; i < exe->nRows; i++)
	{
		int m = INT_MIN;
		for (size_t j = 0; j < cri->nCols; j++)
		{
			m = max(cri->value[i][j], m);
		}
		for (size_t j = 0; j < exe->nCols; j++)
		{
			exe->value[i][j] = (cri->value[i][j] == m);
		}
	}
	return exe;
}

///////////////////////////////////////////////////////////////////////////////

Matrix* affinityPropagation(Matrix* data)
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

	computeCriterion();
	Matrix* exemplars = extractExemplars();

	deleteMatrix(sim);
	deleteMatrix(res);
	deleteMatrix(ava);
	deleteMatrix(cri);

	return exemplars;
}

void affinityPropagationDebug()
{
	Matrix* data = createMatrix(5, 5);
	Matrix* expected = createMatrix(5, 5);
	int staticData[5][5] =
	{
		{3, 4, 3, 2, 1},
		{4, 3, 5, 1, 1},
		{3, 5, 3, 3, 3},
		{2, 1, 3, 3, 2},
		{1, 1, 3, 2, 3}
	}; 
	int staticResult[5][5] =
	{
		{1, 0, 0, 0, 0},
		{1, 0, 0, 0, 0},
		{1, 0, 0, 0, 0},
		{0, 0, 0, 1, 0},
		{0, 0, 0, 1, 0}
	};

	for (size_t i = 0; i < 5; i++)
	{
		memcpy(data->value[i], staticData[i], 5 * sizeof(int));
		memcpy(expected->value[i], staticResult[i], 5 * sizeof(int));
	}

	Matrix* result = affinityPropagation(data);
	printf(equalsMatrix(expected, result) ? "test ok\n" : "test failed\n");

	deleteMatrix(data);
	deleteMatrix(expected);
	deleteMatrix(result);
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
					error("Failed to create Matrix.value[i]");
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
	for (size_t i = 0; i < m->nRows; i++)
	{
		free(m->value[i]);
	}
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

struct size_file
{
	int row;
	int col;
};

struct size_file findSize(FILE* fdata)
{
	rewind(fdata);
	struct size_file res = { .row = 0, .col = 0};
	char line[1024];
	char* token;

	while (fgets(line, sizeof(line), fdata))
	{
		res.row++;
		res.col = 0;
		token = strtok(line, ";");

		while (token != NULL)
		{
			res.col++;
			token = strtok(NULL, ";");
		}
	}
	rewind(fdata);
	return res;
}

Matrix* loadData(char* path)
{
	FILE* fdata = fopen(path, "r");
	if (fdata != NULL)
	{
		struct size_file len = findSize(fdata);

		Matrix* data = createMatrix(len.row, len.col);
		char line[1024];
		char* token;
		int i = 0;
		int j = 0;

		while (fgets(line, sizeof(line), fdata))
		{
			token = strtok(line, ";");
			j = 0;

			do
			{
				data->value[i][j++] = (int)strtol(token, NULL, 10);
			} while ((token = strtok(NULL, ";")) != NULL);
			i++;
		}
		fclose(fdata);
		return data;
	}
	else
	{
		error("Failed to open file.");
		return NULL;
	}
}