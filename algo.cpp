#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <set>
#include <chrono>

/*
COMMAND LINE ARGUMENTS FORMAT EXPECTED:
	./prgm_name <#rows of mat1> <#columns> <sparsity> <#rows of mat2> <#columns> <sparsity>
PLS FOLLOW FORMAT
*/

typedef std::vector<std::pair<int, int>> sparse; //pair has (position in matrix, data value in that position)
//position of matrix[i][j] (matrix is dense) is defined as i*C+j where C is #columns in matrix 

void uniqueRand(sparse &v, int n, int tot) {
	std::set<int> s; //to sort according to position
	while(s.size() < n) {
		s.insert(std::rand() % tot);  //randomly assigning positions for the non-zero data
	}
	for(auto &i: s) {
		v.push_back(std::make_pair(i, rand() % 5 + 1));  //assigning the non-zero data with random values
	}
}

void makeTranspose(sparse &v, int R, int C) {
	//my multiplication algo needs the mat2 transposed and entries sorted by position
	for(int i = 0; i < v.size(); i++) {
		int row = v[i].first/C;
		int col = v[i].first%C;
		v[i].first = col * R + row;
	}
	std::sort(v.begin(), v.end());
}


int main(int argc, char* argv[]) {
	std::srand(static_cast<unsigned int>(std::time(nullptr)));

	if(argc < 7) {
		std::cout<<"\nERR: Enter all arguments";
		return 0;
	}


	int R1 = atoi(argv[1]), C1 = atoi(argv[2]);
	int R2 = atoi(argv[4]), C2 = atoi(argv[5]);
	
	if(R2 != C1) {
		std::cout<<"\nERR: Cannot multiply with these dimensions";
		return 0;
	}

	int N1 = R1 * C1, N2 = R2 * C2;
	float sp1 = atof(argv[3]), sp2 = atof(argv[6]);
	sparse f1, f2, f3;

	uniqueRand(f1, sp1 * N1, N1);
	uniqueRand(f2, sp2 * N2, N2);

	void normalMatMul(sparse, int, int, sparse, int, int); //ONLY for timing comparison, NOT part of my algo 
	normalMatMul(f1, R1, C1, f2, R2, C2); //I repeat, the matrices are still sparse only

	// UNCOMMENT TO PRINT INPUT SPARSE MATRICES
	// std::cout<<"\nMatrix 1:\n";
	// for(int k = 0; k < f1.size(); k++) {
	// 	std::cout<<f1[k].first<<":"<<f1[k].second<<"\n";
	// }
	// std::cout<<"Matrix 2:\n";
	// for(int k = 0; k < f2.size(); k++) {
	// 	std::cout<<f2[k].first<<":"<<f2[k].second<<"\n";
	// }

	int i = 0, j = 0;  //indices for f1 and f2
	int r = 0, c = 0;  //keep track of which row of mat1 and column of mat2 currently
	int n1 = f1.size(), n2 = f2.size();

	auto start = std::chrono::high_resolution_clock::now();

	makeTranspose(f2, R2, C2);  //we need the second matrix to be transposed

	auto transp = std::chrono::high_resolution_clock::now();

	while(r < R1) {
	    int s = 0, startRow = i;   //keeping track of where this row starts
	    while(i < n1 && j < n2 && f1[i].first/C1 == r && f2[j].first/R2 == c) {
            //std::cout<<"i:"<<i<<" j:"<<j<<"\n";
	        if(f1[i].first % C1 == f2[j].first % R2) {
	            s += f1[i].second * f2[j].second;   //multiply!
                j++;
	        } else if(f1[i].first % C1 > f2[j].first % R2) {
	            j++;
	        } else {
	            i++;
	        }
	    }
	    if(s != 0) {
    	    f3.push_back(std::make_pair(C2 * r + c, s));    //we have a non-zero entry in the product!
	    }
        if (i == n1){ //last row, any more new columns?
	        while(j < n2 && f2[j].first/R2 == c) j++;  //look for entry of next non-zero column of mat2
			if(j == n2) r = R1; //done, no more new columns
	        else { //prep for last row of mat1 with the next non-zero column of mat2
				c = f2[j].first/R2;
	        	i = startRow;
			}
	    } else if(j == n2) {  //last column, any more new rows?
	        while(i < n1 && f1[i].first/C1 == r) i++;  //look for entry of next non-zero row of mat1
	        if(i == n1) r = R1; //done, no more new rows
	        else { //prep for a new row of mat1 to be multiplied with all columns from start of mat2
				r = f1[i].first/C1;
	        	c = 0; j = 0;
			}
	    } else if(f1[i].first/C1 != r) {
	        r = f1[i].first/C1;  //new non-zero row, go to start of columns of mat2
	        c = 0; j = 0;
	    } else if(f2[j].first/R2 != c) {
	        i = startRow;        //new non-zero column, go to start of current row of mat1
	        c = f2[j].first/R2;
	    }
	}

	auto end = std::chrono::high_resolution_clock::now();

	std::chrono::duration<float> dur_with_transp = end - start;
	std::chrono::duration<float> dur_without_transp = end - transp;
	
	// UNCOMMENT TO PRINT RESULT
	// std::cout<<"Resultant Sparse Matrix:\n";
	// for(int k = 0; k < f3.size(); k++) {
	//     std::cout<<f3[k].first<<":"<<f3[k].second<<"\n";
	// }
	std::cout<<"\nMemory Ratio:"<<2.0*(f1.size()+f2.size()+f3.size())/float(N1+N2+R1*C2);
	std::cout<<"\nTime with transpose = "<<dur_with_transp.count()*1e3;
	std::cout<<"\nTime without transpose = "<<dur_without_transp.count()*1e3;

	return 0;
}

void normalMatMul(sparse v1, int R1, int C1, sparse v2, int R2, int C2) {
	std::vector<std::vector<int>> a(R1, std::vector<int>(C1));
	std::vector<std::vector<int>> b(R2, std::vector<int>(C2));
	std::vector<std::vector<int>> c(R1, std::vector<int>(C2));
	for(int i = 0; i < v1.size(); i++) {
		int row = v1[i].first/C1;
		int col = v1[i].first%C1;
		a[row][col] = v1[i].second;
	}
	for(int i = 0; i < v2.size(); i++) {
		int row = v2[i].first/C2;
		int col = v2[i].first%C2;
		b[row][col] = v2[i].second;
	}

	auto start = std::chrono::high_resolution_clock::now();
	for(int i = 0; i < R1; i++) {
		for(int j = 0; j < C2; j++) {
			c[i][j] = 0;
			for(int k = 0; k < C1; k++) {
				c[i][j] += a[i][k] * b[k][j];
			}
		}
	}
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> dur = end - start; 
	std::cout<<"\nTime normalMatMul = "<<dur.count()*1e3;

	// UNCOMMENT THIS IF U WANT TO PRINT THE DENSE MATRIX AND CHECK RESULTS
	// std::cout<<"\nDense Product Matrix:";
	// for(int i = 0; i < R1; i++) {
	// 	std::cout<<"\n";
	// 	for(int j = 0; j < C2; j++) {
	// 		std::cout<<c[i][j]<<" ";
	// 	}
	// }
}