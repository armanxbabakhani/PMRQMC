#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <complex>
#include <algorithm>
#include <numeric>
#include <cmath>

using namespace std;


template<typename T>
void Print_matrix(const vector<vector<T>>& matrix) {
    int m = matrix.size();
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

// This function mod p multiplies the matrix A by a scalar c!
vector<vector<int>> Modp_scalmult(const int c , const vector<vector<int>> A , int p){
    int rowsA = A.size() , colsA = A[0].size();
    vector<vector<int>> cA(rowsA , vector<int> (colsA , 0));
    for(int i = 0; i < rowsA ; i++){
        for(int j = 0; j < colsA ; j++){
            if(c >= 0){
                cA[i][j] = c*A[i][j] % p;
            }
            else{
                cA[i][j] = (p + (c*A[i][j] % p)) % p;
            }
        }
    }
    return cA;
}

// This function computes matrix multiplication of A and b mod p.
vector<vector<int>> Modp_mult(const vector<vector<int>> A , const vector<vector<int>> B , int p){
    int rowsA = A.size(), colsA = A[0].size() , rowsB = B.size(), colsB = B[0].size();

    if (colsA != rowsB) {
        cerr << "Matrices are not compatible for multiplication." << std::endl;
        return {};
    }

    vector<vector<int>> AB(rowsA, vector<int>(colsB, 0));
    
    // Perform element-wise multiplication and summation
    for (int i = 0; i < rowsA; ++i) {
        for (int j = 0; j < colsB; ++j) {
            for (int k = 0; k < colsA; ++k) {
                AB[i][j] += A[i][k] * B[k][j];
            }
            AB[i][j] = AB[i][j] % p;
        }
    }

    return AB;
}

// This function divides num by divider (mod p)
int Modp_divide(int num , int divider , int p){
    int answer;
    for (int i = 0; i < p; i++){
        if( (divider * i) % p == num ){
            answer = i;
        }
    }
    return answer;
}

// This function transposes matrix A
vector<vector<int>> Transpose(vector<vector<int>> A){
    vector<vector<int>> A_transpose;
    int rowsA = A.size(), colsA = A[0].size();
    for(int j=0; j < colsA ; j++){
        vector<int> A_transpose_j(rowsA , 0);
        for(int i = 0; i < rowsA ; i++){
            A_transpose_j[i] = A[i][j];
        }
        A_transpose.push_back(A_transpose_j);
    }
    return A_transpose;
}

// Divide_vector divides the given matrix in 2 vertically!
vector<vector<vector<int>>> Divide_vector(const vector<vector<int>> Vec , int vecDim){
    int rowsVec = Vec.size();
    vector<vector<int>> Vec1, Vec2;
    vector<vector<vector<int>>> VecDiv;

    for (int i = 0; i < rowsVec ; i++){
        if(i < vecDim){
            Vec1.push_back(Vec[i]);
        }
        else{
            Vec2.push_back(Vec[i]);
        }
    }

    VecDiv.push_back(Vec1);
    VecDiv.push_back(Vec2);
    
    return VecDiv;
}

vector<vector<vector<int>>> Divide_blocks(const vector<vector<int>> A , int blockDim){
    int rowsA = A.size(), colsA = A[0].size();

    vector<vector<vector<int>>> BlockMatrix;
    vector<vector<int>> TopLeft(blockDim , vector<int> (blockDim , 0)), BottomRight(rowsA - blockDim , vector<int> (rowsA - blockDim , 0));
    vector<vector<int>> TopRight(blockDim , vector<int> (rowsA - blockDim , 0));

    for (int i = 0; i < blockDim; i++) {
        for (int j = i; j < blockDim; j++) {
            TopLeft[i][j] = A[i][j];
        }
    }
    for (int i = blockDim; i < rowsA; i++) {
        for (int j = i; j < rowsA; j++) {
            BottomRight[i - blockDim][j - blockDim] = A[i][j];
        }
    }
    for (int i = 0; i < blockDim; i++) {
        for (int j = blockDim; j < rowsA; j++) {
            TopRight[i][j - blockDim] = A[i][j];
        }
    }

    BlockMatrix.push_back(TopLeft);
    BlockMatrix.push_back(TopRight);
    BlockMatrix.push_back(BottomRight);

    return BlockMatrix;
}

vector<vector<vector<int>>> Divide_blocks_lower(const vector<vector<int>> A , int blockDim){
    int rowsA = A.size(), colsA = A[0].size();

    vector<vector<vector<int>>> BlockMatrix;
    vector<vector<int>> TopLeft(blockDim , vector<int> (blockDim , 0)), BottomRight(rowsA - blockDim , vector<int> (rowsA - blockDim , 0));
    vector<vector<int>> BottomLeft(rowsA-blockDim , vector<int> (blockDim , 0));

    for (int i = 0; i < blockDim; i++) {
        for (int j = 0; j <= i; j++) {
            TopLeft[i][j] = A[i][j];
        }
    }
    for (int i = blockDim; i < rowsA; i++) {
        for (int j = 0; j <= i; j++) {
            BottomRight[i - blockDim][j - blockDim] = A[i][j];
        }
    }
    for (int i = blockDim; i < rowsA; i++) {
        for (int j = 0; j < blockDim; j++) {
            BottomLeft[i-blockDim][j] = A[i][j];
        }
    }

    BlockMatrix.push_back(TopLeft);
    BlockMatrix.push_back(BottomLeft);
    BlockMatrix.push_back(BottomRight);

    return BlockMatrix;
}

// This function vertically concatenates to matrices X1 and X2
vector<vector<int>> Vert_conc(vector<vector<int>> X1 , vector<vector<int>> X2){
    vector<vector<int>> X = X1;
    for(int i = 0; i < X2.size(); i++){
        X.push_back(X2[i]);
    }
    return X;
}

// This function horizontally concatenates matrices X1 and X2
vector<vector<int>> Horz_conc(vector<vector<int>> X1 , vector<vector<int>> X2){
    vector<vector<int>> X = X1;
    if(X1.size() == 0){
        return X2;
    }
    else if(X2.size() == 0){
        return X1;
    }
    else{
        for(int i = 0; i < X.size(); i++){
            for(int j = 0; j < X2[0].size(); j++){
                X[i].push_back(X2[i][j]);
            }
        }
        return X;
    }
}

vector<vector<int>> Vec_sub(vector<vector<int>> Vec1, vector<vector<int>> Vec2 , int p){
    vector<vector<int>> Vecdiff;
    int diff;
    for(int i =0; i < Vec1.size() ; i++){
        vector<int> Vecdiff_i = {};
        for(int j = 0; j < Vec1[0].size(); j++){
            diff = (Vec1[i][j] - Vec2[i][j])%p;
            if(diff < 0){
                diff += p;
            }
            Vecdiff_i.push_back(diff);
        }
        Vecdiff.push_back(Vecdiff_i);
    }
    return Vecdiff;
}

vector<vector<int>> Vec_add(vector<vector<int>> Vec1, vector<vector<int>> Vec2 , int p){
    vector<vector<int>> Vecsum;
    int sum;
    for(int i =0; i < Vec1.size() ; i++){
        vector<int> Vecsum_i = {};
        for(int j = 0; j < Vec1[0].size(); j++){
            sum = (Vec1[i][j] + Vec2[i][j])%p;
            Vecsum_i.push_back(sum);
        }
        Vecsum.push_back(Vecsum_i);
    }
    return Vecsum;
}

bool All_zeros(vector<vector<int>> A){
    for(int i = 0; i < A.size(); i++){
        for(int j = 0; j < A[0].size(); j++){
            if(A[i][j] != 0){
                return false;
            }
        }
    }
    return true;
}

bool Num_found_inVec(int num , vector<int> Vec){
    bool found = false;
    for (int i=0; i < Vec.size(); i++){
        if(num == Vec[i]){
            found = true;
            break;
        }
    }
    return found;
}

typedef vector<vector<int>> Matrix;

// Mod p Gaussian Elimination:
void Modp_GE(Matrix &A , int p) {
    int rowsA = A.size();
    int colsA = A[0].size();
    int lead = 0;

    for (int r = 0; r < rowsA; r++) {
        if (lead >= colsA) {
            break;
        }
        int i = r;
        while (A[i][lead] == 0) {
            i++;
            if (i == rowsA) {
                i = r;
                lead++;
                if (lead == colsA) {
                    return;
                }
            }
        }
        
        swap(A[i], A[r]);

        int lv = A[r][lead] , multiplier = Modp_divide(1, A[r][lead], p); // Multiplier is the number that multiplies the leading term to 1.
        for (int j = 0; j < colsA; j++) {
            A[r][j] = (A[r][j] * multiplier) % p;
        }
        for (int i = 0; i < rowsA; i++) {
            if (i != r) {
                int lv2 = A[i][lead];
                for (int j = 0; j < colsA; j++) {
                    A[i][j] = (A[i][j] - (A[r][j] * lv2) % p + p) % p;
                }
            }
        }
        lead++;
    }
}

int Sum(vector<int> A, int p){
    int sum = 0;
    for(int i=0; i < A.size() ; i++){
        sum =+ A[i];
    }
    return sum%p;
}

// The following function decomposes an integer into its prime factors
vector<pair<int , int>> Prime_decomp(int n){
    vector<pair<int, int>> primes;

    int power2 = 0;
    if(n%2 == 0){
        while(n % 2 == 0){
            n /= 2;
            power2++;
        }
        primes.push_back({2 , power2});
    }

    for(int i = 3; i*i < n ; i += 2){
        int poweri = 0;
        if(n%i == 0){
            while(n % i == 0){
                n /= i;
                poweri++;
            }
            primes.push_back({i , poweri});
        }
    }

    // In case n is itself a prime:
    if(n > 1){
        primes.push_back({n , 1});
    }

    return primes;
}

// Nullspace for finite fields (mod p, where p is a prime)!
vector<vector<int>> Modp_Nullspace(vector<vector<int>> A , int p){
    Modp_GE(A , p);
    vector<vector<int>> nullspace;

    int colsA = A[0].size() , rowsA = A.size() , rank=0 , nullspaceDim;

    vector<int> PivotCols , PivotRows;
    //Determine pivot columns:
    for(int i = 0; i < rowsA; i++){
        for(int j = 0; j < colsA; j++){
            if(A[i][j] == 1){
                PivotCols.push_back(j);
                PivotRows.push_back(i);
                rank++;
                break;
            }
        }
    }

    nullspaceDim = colsA - rank;

    for(int i=0; i < colsA; i++){
        if(!Num_found_inVec(i, PivotCols)){ 
            vector<int> BasisVec(A[0].size() , 0);
            BasisVec[i] = 1;
            // Mark the rows that have a non-zero entry in this column:
            for(int k = 0; k < rank; k++){
                BasisVec[PivotCols[k]] = (p - A[PivotRows[k]][i])%p;
            }
            nullspace.push_back(BasisVec);
        }
    }
    return Transpose(nullspace);
}

// The following function solves the system of equations Ax = b (mod p)
// Even though b is a vector<vector<int>>, it must be only a single column vector!
vector<vector<int>> Modp_solver(vector<vector<int>> A , vector<vector<int>> b , int p){
    vector<vector<int>> Ab = Horz_conc(A , Modp_scalmult(-1 , b , p));
    Modp_GE(Ab , p);
    vector<vector<int>> x = Modp_Nullspace(Ab , p) , x_final;
    int n = x[0].size() , m = x.size();
    for(int i = 0; i < m; i++){
        if(x[i][n-1] != 0){
            int mult = Modp_divide(x[i][n-1] , 1 , p);
            vector<int> x_i = Modp_scalmult(mult , {x[i]} , p)[0];
            x_i.erase(x_i.begin() + n - 1);
            x_final.push_back(x_i);
        }
    }
    return x_final;
}

// This function computes nullspace basis of A mod p^r! when r = 1, it simply makes a single call to Modp_Nullspace(A , p)
vector<vector<int>> Modp_Nullspace_r(vector<vector<int>> A , int p , int r){
    if(r == 1){
        return Modp_Nullspace(A , p);
    }
    else{
        vector<vector<int>> x = Modp_Nullspace_r(A , p , r-1);
        vector<vector<int>> b2 = Modp_scalmult(-1 , Modp_mult(A , x , pow(p , r)), pow(p , r));

        x = Modp_scalmult(p , x , pow(p,r));
        vector<vector<int>> b2trans = Transpose(b2);

        for(int i = 0; i < b2trans.size(); i++){
            vector<vector<int>> xr = Transpose(Modp_solver(A , Transpose({b2trans[i]}) , r-1));
            if(!All_zeros(xr)){
                x = Horz_conc(x , xr);
            }
        }
        return x;
    }
}


// This function computes the mod n nullspace basis of A, where n is any integer!
vector<vector<int>> Nullspace_n(const vector<vector<int>> A , int n){
    vector<pair<int, int>> ps = Prime_decomp(n);
    vector<vector<int>> Null;

    if(A.size() > 0){
        for(int i=0 ; i < ps.size() ; i++){
            int prime = ps[i].first , power = ps[i].second;
            vector<vector<int>> Nullspacei = Modp_Nullspace_r(A , prime , power);
            Nullspacei = Modp_scalmult(n/pow(prime , power) , Nullspacei , n);
            Null = Horz_conc(Null , Nullspacei);
        }
    }
    return Null;
}

int main(){
    // Testing Nullspace finder!
    vector<vector<int>> A = Transpose({{1 , 2 , 1} , {2 , 1 , 0} , {0 , 2  , 1} , {1 , 0 , 1}}), Arr;
    
    Arr = A;
    int p = 7;
    Modp_GE(Arr , p); // Gaussian Elimination performed on A!
    cout << "A is " << endl;
    Print_matrix(A);
    cout << endl;

    cout << "A in row echelon form " << endl;
    Print_matrix(Arr);

    cout << endl;
    vector<vector<int>> EmptyMat;
    vector<vector<int>> x = Nullspace_n(EmptyMat , p);
    cout << "The nullspace basis of A " << endl;
    Print_matrix(x);

    return 0;
}