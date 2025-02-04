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

void Print_data(Coeffs C , DVecs D , vector<vector<pair<int, int>>> P){
    for(int i = 0; i < P.size() ; i++){
        cout << "The permutation is ";
        for(int l = 0; l < P[i].size(); l++){
            cout << "< spin # " << P[i][l].first << " , P^" << P[i][l].second << " > ";
        }
        cout << endl;
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

vector<pair<complex<double>, vector<int>>> data_extract(const string& fileName){
    vector<pair<complex<double>, vector<int>>> data;

    ifstream inputFile(fileName);
    if (!inputFile) {
        cout << "Failed to open the input file!" << endl;
        return data;
    }

    string line;
    while (getline(inputFile, line)) {
        // The first non-empty element is the coefficient
        istringstream iss(line);
        pair<complex<double>,vector<int>> linedata;

        // Extracting the complex coefficient:
        double realpart, imagpart=0;
        char sign;
        string complexPart;
        iss >> complexPart; 

        istringstream complexIss(complexPart);
        complexIss >> realpart >> imagpart;
        linedata.first = complex<double> (realpart , imagpart);

        //Extracting the integer vectors of qubits and paulis:
        string token;
        vector<int> integers;
        while (iss >> token){
            // int num = std::stoi(token);
            int num = token == "X" || token == "x" ? 1 : 
                      token == "Y" || token == "y" ? 2 : 
                      token == "Z" || token == "z" ? 3 : std::stoi(token);
            integers.push_back(num);
        }
        linedata.second = integers;
        data.push_back(linedata);
        }
        inputFile.close();
    return data;
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
    DVecs DiagonalSet0 = DiagonalSet;
    Coeffs CoeffSet0 = CoeffSet;
    int len = PermutationSet.size();
    PermutationSet.clear();
    DiagonalSet.clear();
    CoeffSet.clear();

    for(int j=0; j < NewPermutations.size(); j++){
        vector<vector<pair<int, int>>> PermutationSetj = PermutationSet0;
        DVecs DiagonalSetj = DiagonalSet0;
        Coeffs CoeffSetj = CoeffSet0;
        for(int k=0; k < len; k++){
            // Maybe: Check for uniqueness here OR do it during the concatenation of two consecutive lines!
            PermutationSetj[k].push_back({NewParticleNo , NewPermutations[j]}); // Here, we are multiplying the permutation matrices
            Diag_Multiply(DiagonalSetj[k] , CoeffSetj[k] , NewDiagonals[j] , NewCoeffs[j]); // Here, we are multiplying the diagonal matrices with corresponding coefficients!
        }
        PermutationSet.insert(PermutationSet.end() , PermutationSetj.begin() , PermutationSetj.end());
        DiagonalSet.insert(DiagonalSet.end() , DiagonalSetj.begin() , DiagonalSetj.end());
        CoeffSet.insert(CoeffSet.end() , CoeffSetj.begin() , CoeffSetj.end());
    }
} 

struct PDdata {
    vector<vector<pair<int, int>>> Permutations;
    DVecs Diagonals;
    Coeffs Coefficients;
};

bool Perm_compare(vector<pair<int,int>> A , vector<pair<int,int>> B){
    if(A.size() != B.size())
        cerr << "The sizes of the two input vectors are not equal! Comparison failed." << endl;
    for(int i = 0; i < A.size(); i++){
        if(A[i] != B[i]){
            return false;
        }
    }
    return true;
}

pair<bool , int> Find_permutation(vector<pair<int,int>> LinePerms , vector<vector<pair<int,int>>> AllPerms){
    for(int i=0; i<AllPerms.size(); i++){
        if(Perm_compare(LinePerms , AllPerms[i])){
            return {true, i};
        }
    }
    return {0,0};
}
void PMR_append(PDdata& pdData, PDdata pdDataLine){
    vector<vector<pair<int,int>>> AllPerms = pdData.Permutations , LinePerms = pdDataLine.Permutations;
    DVecs AllDiags = pdData.Diagonals , LineDiags = pdDataLine.Diagonals;
    Coeffs AllCoes = pdData.Coefficients , LineCoes = pdDataLine.Coefficients;

    for(int i=0; i < LinePerms.size(); i++){
        pair<bool , int> PermFound = Find_permutation(LinePerms[i] , AllPerms);
        if(PermFound.first){
            AllDiags[PermFound.second].insert(AllDiags[PermFound.second].end() , LineDiags[i].begin() , LineDiags[i].end());
            AllCoes[PermFound.second].insert(AllCoes[PermFound.second].end() , LineCoes[i].begin() , LineCoes[i].end());
        }
        else{
            AllPerms.push_back(LinePerms[i]);
            AllDiags.push_back(LineDiags[i]);
            AllCoes.push_back(LineCoes[i]);
        }
    }
    pdData.Permutations = AllPerms;
    pdData.Coefficients = AllCoes;
    pdData.Diagonals = AllDiags;
}

PDdata CDPconvert(const vector<pair<complex<double>,vector<int>>> data) {
    int NumLines = data.size();
    extern int NumOfParticles;
    PDdata pdData;

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

        vector<vector<pair<int , int>>> PMatricesLine;  // The first in the pair is the set of particle numbers and the second vector is the set of powers of permutations
        DVecs DMatricesLine; //This vector maps the Zs to Ps it is a many to one mapping!
        Coeffs CsLine;
        PDdata pdDataLine;

        int twoSpinplus1l; 
        for (size_t i = 0; i < datal.size()/4; i++) {
            int Particlei = datal[4 * i] , Pauli = datal[4 * i + 1] , Poweri = datal[4 * i + 2] , twoSpinPlusOnei = datal[4 * i + 3];
            twoSpinplus1l = twoSpinPlusOnei; // Assuming only a single spin species for now. In future, if there are multiple spin species, additional matrices should be created to separate these permutation operators!
            if (Particlei > NumOfParticles)
                NumOfParticles = Particlei;
            PauliCDPs Operator;
            if (Pauli == 1){
                Operator = {{1 , twoSpinPlusOnei - 1} , {{Dplus , Dminus} , {one , one}}};  // Operator = S_x!
            } else if (Pauli == 2){
                Operator = {{1 , twoSpinPlusOnei - 1} , {{Dplus , Dminus} , {minusi , plusi}}}; // Operator = S_y!
            } else if (Pauli == 3){
                Operator = {{0} , {{Dz} , {one}}}; // Operator = S_z!
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
                    NewC.push_back({{Operator.second.second[p]}}); // This takes care of the front coefficient multiplication!
                }
                Diagsl.push_back(NewD);
                Coeffsl.push_back(NewC);
                while(Poweri > 1){
                    Permutation_append(Permsl[PartFound.second] , Diagsl[PartFound.second], Coeffsl[PartFound.second] , Operator , twoSpinPlusOnei);
                    Poweri--;
                }
            }
        }
        // Combining all the P and D matrices from a single line into corresponding PRM forms (there could be multiple PMR terms from each line)
        for(int k = 0; k < Permsl.size(); k++){
            if(k == 0){
                for(int kk = 0; kk < Permsl[k].size(); kk++){
                    PMatricesLine.push_back({{ParticlNos[k] , Permsl[k][kk]}});
                    DMatricesLine.push_back(Diagsl[k][kk]);
                    CsLine.push_back(Coeffsl[k][kk]);
                }
            }
            else{
                PMR_otimes(PMatricesLine , DMatricesLine , CsLine , Permsl[k] , Diagsl[k] , Coeffsl[k] , ParticlNos[k] , twoSpinplus1l);
            }
        }
        // Multiply all coefficients by the coefficient at front:
        for(int c = 0; c < CsLine.size(); c++){
            for(int cc = 0; cc < CsLine[c].size(); cc++){
                CsLine[c][cc] = CsLine[c][cc]*coeffl;
            }
        }

        // Adding the new CDP from the line to the entire set of total CDPs
        pdDataLine.Coefficients = CsLine;
        pdDataLine.Permutations = PMatricesLine;
        pdDataLine.Diagonals = DMatricesLine;

        PMR_append(pdData, pdDataLine);
    }

    return pdData;
}

vector<vector<int>> Convert_perms(vector<vector<pair<int,int>>> PMatrices){
    extern int NumOfParticles;
    vector<vector<int>> ColPermMatrix(NumOfParticles , vector<int>(PMatrices.size() , 0));
    for(int i=0; i< PMatrices.size(); i++){
        for(int j = 0; j < PMatrices[i].size(); j++){
            int particle = PMatrices[i][j].first , power = PMatrices[i][j].second;
            ColPermMatrix[particle-1][i] = power;
        }
    }
    return ColPermMatrix;
}


int main(int argc , char* argv[]){
    string fileName(argv[1]);  // Reading the name of the input .txt file describing the Hamiltonian
    vector<pair<complex<double>, vector<int>>> data = data_extract(fileName);
    cout << "Now processing ... " << endl;
    cout << endl;
    PDdata CDPdata = CDPconvert(data);
    Coeffs Cs = CDPdata.Coefficients;
    DVecs DMatrices = CDPdata.Diagonals;
    vector< vector<pair<int,int>>> PMatrices = CDPdata.Permutations;

    cout << "The following is the breakdown of the data " << endl;
    Print_data(Cs , DMatrices , PMatrices);

    // Converting the PMatrices into vector<vector<int>> to make matrix of column permutations

    vector<vector<int>> PermMatrixColumn = Convert_perms(PMatrices);
    cout << endl;
    printMatrix(PermMatrixColumn);

    return 0;
}