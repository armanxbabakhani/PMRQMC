#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <complex>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <bitset>

using namespace std;

int NumOfParticles = 0; // number of particles is a global variable!

typedef vector<vector<complex<double>>> Coeffs;
typedef vector<vector<int>> ZVecs;
typedef pair<int , int> Diag; // A diagonal matrix is specified by two numbers z and k!
typedef vector<vector<vector<Diag>>> DVecs;
typedef pair<vector<int> , pair<vector<Diag> , vector<complex<double>>>> PauliCDPs; // This typedef is to summarize the Paulis as sum of DPs (with corresponding coefficients);
struct PZData {
    vector<vector<int>> Ps;
    Coeffs coeffs;
    // ZVecs Dtrack;
    DVecs DMatrices;
};

template<typename T>
void printMatrix(const vector<vector<T>>& matrix) {
    int m = matrix.size();
    int n = matrix[0].size();

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

void Print_diagonals(vector<vector<Diag>> Diags , vector<complex<double>> Cs){
    for(int j = 0; j < Diags.size(); j++){
        cout <<  Cs[j];
        for(int k = 0; k < Diags[j].size() ; k++){
            if(Diags[j][k].first == 1){
                cout << " D(z," << Diags[j][k].second << ")";
            }
            else{
                cout << " D(" << Diags[j][k].second << ")";
            }
        }
        if(j < Diags.size()-1){
            cout << " + ";
        }
        else{
            cout << endl;
        }
    }
    cout << endl;
}

void Print_data(vector<vector<complex<double>>> C , vector<vector<vector<Diag>>> D , vector<int> P){
    for(int i = 0; i < P.size() ; i++){
        cout << "Permutation index: " << P[i] << endl;
        for(int j = 0; j < D[i].size(); j++){
            cout <<  C[i][j];
            for(int k = 0; k < D[i][j].size() ; k++){
                if(D[i][j][k].first == 1){
                    cout << " D(z," << D[i][j][k].second << ")";
                }
                else{
                    cout << " D(" << D[i][j][k].second << ")";
                }
            }
            if(j < D[i].size()-1){
                cout << " + ";
            }
            else{
                cout << endl;
            }
        }
        cout << endl;
    }
}

// This function finds an instance of an integer in a vector of integers!
pair<bool , int> Find_number(int Num , vector<int> NumVec){
    for(int i = 0; i < NumVec.size() ; i++){
        if(NumVec[i] == Num){
            return {true , i};
        }
    }
    return {false , 0};
}

void Permutation_append(vector<int>& AllPermsOnParticle , vector<vector<vector<Diag>>>& AllDiagsOnParticle , vector<vector<complex<double>>>& AllCoeffsOnParticle , PauliCDPs CurrentPnDnCs , int twosp1){
    vector<int> CurrentPs = CurrentPnDnCs.first;
    vector<Diag> CurrentDs = CurrentPnDnCs.second.first;
    vector<complex<double>> CurrentCs = CurrentPnDnCs.second.second;

    // We have to go through AllPermsOnParticle and generate a new term by adding the c
    vector<int> AllPInit = AllPermsOnParticle;
    AllPermsOnParticle.clear();
    vector<vector<vector<Diag>>> AllDInit = AllDiagsOnParticle;
    AllDiagsOnParticle.clear();
    vector<vector<complex<double>>> AllCoes = AllCoeffsOnParticle;
    AllCoeffsOnParticle.clear();

    // Initializing the diagonal element to be moved to the left of the existing permutations
    Diag CurrDj = {CurrentDs[0].first , 0};

    for(int j=0; j < CurrentPs.size(); j++){
        vector<int> AllP2 = AllPInit;
        vector<vector<vector<Diag>>> AllDiags2 = AllDInit;
        vector<vector<complex<double>>> AllCoes2 = AllCoes;
        for(int i=0; i < AllPInit.size(); i++){
            AllP2[i] = (AllP2[i] + CurrentPs[j]) % twosp1;

            CurrDj.second = (CurrentDs[j].second + AllPInit[i]) % twosp1; // Moving the new D matrix to the right of the existing permutation matrices!
            // Adding the new D matrix to the list of the existing D matrices!
            for(int k = 0; k < AllDiags2[i].size() ; k++){
                AllDiags2[i][k].push_back(CurrDj);
            }
            
            // If the new permutation already exists, only add the new diagonal to the list of existing ones
            pair<bool , int> PermFound = Find_number(AllP2[i] , AllPermsOnParticle); 
            if(PermFound.first){
                for(int k=0; k < AllCoes2[i].size(); k++){
                    AllCoes2[i][k] = AllCoes2[i][k]*CurrentCs[j];
                }
                AllDiagsOnParticle[PermFound.second].insert(AllDiagsOnParticle[PermFound.second].end() , AllDiags2[i].begin() , AllDiags2[i].end());
                AllCoeffsOnParticle[PermFound.second].insert(AllCoeffsOnParticle[PermFound.second].end() , AllCoes2[i].begin() , AllCoes2[i].end());

                // Remove the permutation, and the added diagonals with corresponding coefficients
                AllP2.erase(AllP2.begin() + i);
                AllDiags2.erase(AllDiags2.begin() + i);
                AllCoes2.erase(AllCoes2.begin() + i);
            }
        }
        AllPermsOnParticle.insert(AllPermsOnParticle.end() , AllP2.begin() , AllP2.end());
        AllDiagsOnParticle.insert(AllDiagsOnParticle.end() , AllDiags2.begin() , AllDiags2.end());

        // Multiplying all the existing coefficients with the new coefficients!
        for(int i = 0; i < AllCoes2.size(); i++){
            for(int k = 0; k < AllCoes2.size() ; k++){
                AllCoes2[i][k] = AllCoes2[i][k] * CurrentCs[j];
            }
        }
        AllCoeffsOnParticle.insert(AllCoeffsOnParticle.end() , AllCoes2.begin() , AllCoes2.end());
    }
}

template<typename T>
vector<T> Concat_one(vector<T> A , vector<T> B){
    vector<T> output = A;
    for(int i = 0; i < B.size(); i++){
        output.push_back(B[i]);
    }
    return output;
}

template<typename T>
vector<vector<T>> Concat_two(vector<vector<T>> A , vector<vector<T>> B){
    vector<vector<T>> output = A;
    for(int i = 0; i < B.size(); i++){
        output.push_back(B[i]);
    }
    return output;
}

void Diag_Multiply(vector<vector<Diag>>& ExistingDiags , vector<complex<double>>& ExistingCoeffs , vector<vector<Diag>> NewDiags , vector<complex<double>> NewCoeffs){
    vector<vector<Diag>> ExistingDiags0 = ExistingDiags;
    vector<complex<double>> ExistingCoeffs0 = ExistingCoeffs;
    int NumExisting = ExistingDiags.size() , NumNew = NewDiags.size(); 
    ExistingDiags.clear();
    ExistingCoeffs.clear();
    for(int i = 0; i < NumExisting; i++){
        for(int j = 0; j < NumNew; j++){
            ExistingDiags.push_back(Concat_one(ExistingDiags0[i] , NewDiags[j]));
            ExistingCoeffs.push_back(ExistingCoeffs0[i]*NewCoeffs[j]);
        }
    }
}

void PMR_otimes(vector<vector<pair<int,int>>>& PermutationSet , DVecs& DiagonalSet , Coeffs& CoeffSet , vector<int> NewPermutations , DVecs NewDiagonals, Coeffs NewCoeffs, int NewParticleNo , int twoSplusOne){
    vector<vector<pair<int, int>>> PermutationSet0 = PermutationSet;
    int len = PermutationSet.size();
    for(int j=0; j < NewPermutations.size(); j++){
        vector<vector<pair<int, int>>> PermutationSetj = PermutationSet0;
        for(int k=0; k < len; k++){
            // Maybe: Check for uniqueness here OR do it during the concatenation of two consecutive lines!
            PermutationSetj[k].push_back({NewParticleNo , NewPermutations[j]});
            Diag_Multiply(DiagonalSet[k] , CoeffSet[k] , NewDiagonals[j] , NewCoeffs[j]); // Here DiagonalSet and CoeffSet are appended by the New Diagonal terms!
        }
        // Multiply the corresponding Diagonals and coefficients.
        // Then vertically concatenate PermutationSet with PermutationSetj and the corresponding Diagonals and Coefficients.

    }
} 

struct PDdata {
    vector<vector<int>> Permutations;
    Coeffs coeffs;
    // ZVecs Dtrack;
    DVecs DMatrices;
};

PDdata CDPconvert(const vector<pair<complex<double>,vector<int>>>& data) {
    PDdata pdData;
    int NumLines = data.size();
    extern int NumOfParticles;
    vector<pair<vector<int>,vector<int>>> PMatrices;  // The first in the pair is the set of particle numbers and the second vector is the set of powers of permutations
    Coeffs Coefficients;
    //ZVecs Dtrack;
    DVecs DMatrices; //This vector maps the Zs to Ps it is a many to one mapping!

    // Definining the diagonal matrices and their complex coefficients
    complex<double> one(1,0) , plusi(0,1) , minusi(0,1);
    Diag Dplus = {0 , 0} , Dminus = {0 , -1} , Dz = {1 , 0};

    // Defining the diagonal matrices with the pair convection. If the first integer is zero, we have D^{(k)}, and if the
    //      first integer is one, we have D^{(z,k)}. The second integer specifies k!
    for (int l = 0; l < NumLines; l++){
        complex<double> coeffl = data[l].first;
        vector<int> datal = data[l].second; // Extracts the array of spins and paulis for every line of input!
        vector<int> ParticlNos;
        vector<vector<int>> Permsl; // This is a vector (for each particle number) of (vector<int> The set of permutations acting on that particle)
        vector<DVecs> Diagsl; // The first vector is for each particle, the second vector is for the set of permutations, and the third vector is for the set of diagonal vectors for the specified permutation
        vector<Coeffs> Coeffsl; // Since we can get multiple permutation matrices per line of data, we need to keep track of the coefficients for each 

        for (size_t i = 0; i < datal.size() / 4; i++) {
            int Particlei = datal[4 * i] , Pauli = datal[4 * i + 1] , Poweri = datal[4 * i + 2] , twoSpinPlusOnei = datal[4 * i + 3];
            if (Particlei > NumOfParticles)
                NumOfParticles = Particlei;
            PauliCDPs Operator;
            if (Pauli == 1){
                Operator = {{1 , twoSpinPlusOnei - 1} , {{Dplus , Dminus} , {one , one}}};  // Operator = X!
            } else if (Pauli == 2){
                Operator = {{1 , twoSpinPlusOnei - 1} , {{Dplus , Dminus} , {minusi , plusi}}}; // Operator = Y!
            } else if (Pauli == 3){
                Operator = {{0} , {{Dz} , {one}}}; // Operator = Z!
            }

            pair<bool, int> PartFound = Find_number(Particlei , ParticlNos);
            if(PartFound.first){
                while(Poweri > 0){
                    Permutation_append(Permsl[PartFound.second] , Diagsl[PartFound.second], Coeffsl[PartFound.second] , Operator , twoSpinPlusOnei);
                    Poweri--;
                }
            }
            else{
                ParticlNos.push_back(Particlei);
                Permsl.push_back(Operator.first);
                DVecs NewD;
                Coeffs NewC;
                for(int p = 0; p < Operator.second.first.size() ; p++){
                    NewD.push_back({{Operator.second.first[p]}});
                    NewC.push_back({{Operator.second.second[p]}});
                }
                Diagsl.push_back(NewD);
                Coeffsl.push_back(NewC);
                while(Poweri > 1){
                    Permutation_append(Permsl[PartFound.second] , Diagsl[PartFound.second], Coeffsl[PartFound.second] , Operator , twoSpinPlusOnei);
                    Poweri--;
                }
            }
        }



    }

    return pdData;
}


int main(){
    int twosp1 = 3;
    Diag Dplus = {0 , 0} , Dminus = {0 , -1} , Dz = {1 , 0}; 
    complex<double> one(1.0 , 0.0) , plusi(0.0 , 1.0) , minusi(0.0 , -1.0);
    PauliCDPs X = {{1 , twosp1 - 1} , {{Dplus , Dminus} , {one , one}}} , Y = {{1 , twosp1 - 1} , {{Dplus , Dminus} , {minusi , plusi}}} , Z = {{0} , {{Dz} , {one}}};
    vector<int> P = {1 , 2 , 0};
    vector<vector<vector<Diag>>> D = {{{Dplus}} , {{Dminus}} , {{Dz} , {Dminus , Dplus}}};
    vector<vector<complex<double>>> C = {{plusi} , {minusi} , {one , one}};

    cout << "The data before multiplication" << endl;
    Print_data(C , D , P);

    Permutation_append(P , D , C , X , twosp1);

    return 0;
}