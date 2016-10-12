#include "readFromFile.h"
using namespace std;

/* ****************************************************************** */
int readInfo(Grid* grid, std::string file_path) {
	std::ifstream infofile, supportfile, typefile;
	if(openFlatFile(file_path + "input/info.txt", &infofile)){return 1;};
	if(openFlatFile(file_path + "input/support.txt", &supportfile)){return 1;};
	if(openFlatFile(file_path + "input/types.txt", &typefile)){return 1;};

	int nf, nc, n_sup;

	// Read sizes
	infofile >> nf;
	infofile >> nc;

	// Read offsets
	(*grid).offsets = new int[nc + 1];
	int i = 0;
	while (!infofile.eof() && i < nc + 1) {
		infofile >> (*grid).offsets[i];
		i++;
	}

	n_sup = (*grid).offsets[nc];
	// Read support
	(*grid).support = new int[n_sup];
	(*grid).boundary = new bool[n_sup]();
	int* celltypes = new int[n_sup];
	i = 0; int j =1; int n_Two = 0;
	while (!supportfile.eof() && !typefile.eof() && i < n_sup) {
		supportfile >> (*grid).support[i];
		typefile >> celltypes[i];
		i++;
	}

	for (int i = 0; i < n_sup; i++ ){
		if ( celltypes[i] == 2 ){
			(*grid).boundary[i] = true;
		}
	}

  	(*grid).N = nf;
	(*grid).M = nc;
  	(*grid).n_sup = n_sup;
	(*grid).basis = new double[n_sup];

	if(readBasisOperator(grid, (*grid).basis, file_path)) return 1;

	delete[] celltypes;
	return 0;
};

/* ****************************************************************** */
int readBasisOperator(Grid * grid, double* basis, std::string file_path) {
	int n_el = (*grid).n_sup;
	std::ifstream f;
	openFlatFile(file_path + "input/operator.txt", &f);
	int i = 0;
	while (!f.eof() && i < n_el) {
		f >> basis[i];
		i++;
	}
	return 0;
};

/* ****************************************************************** */
int readMatrix(Grid * grid, Matrix * mat, std::string file_path) {
	std::ifstream f, fs;
	if(openFlatFile(file_path + "input/matrix.txt", &f)){return 1;};
	if(openFlatFile(file_path + "input/sparsity.txt", &fs)){return 1;};
	int maxRow, n_el;
	n_el = (*grid).N;
	f >> maxRow;

	(*mat).mat_coef    = new double[n_el*maxRow];
	(*mat).mat_index = new int[n_el*maxRow];
	(*mat).N = n_el;
	(*mat).maxRow = maxRow;

	int i = 0;
	int ix;
	while (!f.eof() && i < n_el) {
		for (int j = 0; j < maxRow; j++) {
			ix = i*maxRow + j;
			f >> (mat->mat_coef)[ix];
			fs >> (mat->mat_index)[ix];
		}
		i++;
	}

	for (int i = 0; i < n_el; i++) {
		for (int j = 0; j < maxRow; j++) ix = i*maxRow + j;
	}

	return 0;
};

// ********************************************************************************************** //
int writeBasis(Grid& g_grid, std::string file_path){
	std::ofstream outfile;
	if(openFlatFile(file_path + "output/basis.txt", &outfile)){return 1;};

    for (int i = 0; i < g_grid.n_sup; i++){
        outfile<<g_grid.basis[i]<<" ";
    }
	outfile.close();
}

/* ****************************************************************** */
int openFlatFile(std::string file_path, std::ifstream * file) {
	(*file).open(file_path.c_str());
	if (!(*file).is_open()) {
		std::cout << "Unable to open file: " << file_path << "\r\n";
		return 1;
	}
	return 0;
}

/* ****************************************************************** */
int openFlatFile(std::string file_path, std::ofstream * file) {
	(*file).open(file_path.c_str());
	if (!(*file).is_open()) {
		std::cout << "Unable to open file: " << file_path << "\r\n";
		return 1;
	}
	return 0;
}
