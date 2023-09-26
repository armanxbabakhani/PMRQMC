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

// This function converts the array of integers into a corresponding binary string (used to convert the indices of Z into string of bitsets)
string int_to_str(vector<int> Z){
    string Z_string = "";
    if(Z.size() < 1){
        return "0";
    }
    std::sort(Z.begin() , Z.end());
    int Z_size = Z.size();
    int count = 1, indZcount = 0, max_z = Z[Z_size-1];

    while(indZcount < Z_size){
	if(Z[indZcount] < count){
		cout << endl << "Error: repeating spin indices are detected in a Pauli string" << endl; exit(1);
	} else if(Z[indZcount] == count){
            Z_string = "1" + Z_string;
            indZcount++;
        }
        else{
            Z_string = "0" + Z_string;
        }
        count++;
    }
    return Z_string;
}

// This function extracts the input file information into a vector of pairs!
// The assumption about the input file:
//      coeff_1  (particle-number_1)_1  (pauli-number_1)_1  (power-of-pauli_1)_1  ((2*spin+1)_1)_1 ... (particle-number_1)_k  (pauli-number_1)_k  (power-of-pauli_1)_k  ((2*spin+1)_1)_k

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


typedef vector<complex<double>> Coeffs;
typedef vector<vector<int>> ZVecs;
typedef pair<int , int> Diag; // A diagonal matrix is specified by two numbers z and k!
typedef vector<vector<pair<int , vector<pair<int , int>>>>> DVecs;
typedef pair<vector<int> , pair<vector<Diag> , vector<complex<double>>>> PDCs; // This typedef is to summarize the Paulis as sum of DPs (with corresponding coefficients);
struct PZData {
    vector<set<int>> Ps;
    Coeffs coeffs;
    ZVecs Dtrack;
    DVecs DMatrices;
};


void Permutation_append(vector<int>& AllPermsOnParticle , vector<vector<vector<Diag>>>& AllDiagsOnParticle , vector<vector<complex<double>>>& AllCoeffsOnParticle , PDCs CurrentPnDnCs , int spin){
    vector<int> CurrentPs = CurrentPnDnCs.first , AllPInit = AllPermsOnParticle;
    vector<Diag> CurrentDs = CurrentPnDnCs.second.first;
    vector<Complex<double>> CurrentCs = CurrentPnDnCs.second.second;

    // We have to go through AllPermsOnParticle and generate a new term by adding the c
    vector<int> AllPInit = AllPermsOnParticle;
    AllPermsOnParticle.clear();
    vector<vector<vector<Diag>>> AllDInit = AllDiagsOnParticle;
    AllDiagsOnParticle.clear();
    vector<vector<complex<double>>> AllCoes = AllCoeffsOnParticle;
    AllCoeffsOnParticle.clear();

    Diag CurrDj = {CurrentDs[0].first , 0};

    for(int j=0; j < CurrentPs.size(); j++){
        vector<int> AllP2 = AllPInit;
        vector<vector<vector<Diag>>> AllDiags2 = AllDInit;
        vector<complex<double>> AllCoes2 = AllCoes;
        for(int i=0; i < AllPermsOnParticle.size() ; i++){
            AllP2[i] = (AllP2[i] + CurrentPs[j]) % spin;
            CurrDj.second = (CurrrentDs[j].second + AllPerms[i]) % spin; // Moving the new D matrix to the right of the existing permutation matrices!
            
            // Adding the new D matrix to the list of the existing D matrices!
            for(int k = 0; k < AllDiags2[i].size() ; k++){
                AllDiags[i][k].push_back(CurrDj);
            }
        }
        AllPermsOnParticle.insert(AllCoeffsOnParticle.end() , AllP2);
        AllDiagsOnParticle.insert(AllDiagsOnParticle.end() , AllDiags2);

        // Multiplying all the existing coefficients with the new coefficients!
        for(int i = 0; i < AllCoes2.size(); i++){
            for(int k = 0; k < AllCoes2.size() ; k++){
                AllCoes2[i][k] *= CurrentCs[j];
            }
        }
        AllCoeffsOnParticle.insert(AllCoeffsOnParticle.end() , AllCoes2);
    }
}

PZdata PZcomp(const vector<pair<complex<double>,vector<int>>>& data) {
    PZdata PZ_data;
    int l = data.size(), zcount = 0;
    extern int no_qubit;
    vector<vector<int>> Ps;
    Coeffs coeffs;
    ZVecs Dtrack;
    DVecs DMatrices; //This vector maps the Zs to Ps it is a many to one mapping!
    pair<int,int> Dplus = {0 , 0} , Dminus = {0 , -1} , Dz = {1 , 0}; 
    complex<double> one(1,0) , plusi(0,1) , minusi(0,1);

    // Defining the diagonal matrices with the pair convection. If the first integer is zero, we have D^{(k)}, and if the
    //      first integer is one, we have D^{(z,k)}. The second integer specifies k!
    for (int i = 0; i < l; i++){
        complex<double> coeff_i = data[i].first;
        vector<int> zsi; // For every line zs extracts the spins on which a pauli Z acts! 
        vector<int> datai = data[i].second; // Extracts the array of spins and paulis for every line of input!
        // The following data structure is to 
        vector<int> ParticlNos;
        vector<vector<int>> Permsi; // This is a vector (for each particle number) of (vector<int> The set of permutations acting on that particle)
        vector<vector<vector<Diag>>> Diagsi; // The first vector is for each particle, the second vector is for the set of permutations, and the third vector is for the set of diagonal vectors for the specified permutation
        vector<vector<complex<double>>> coeffsi; // Since we can get multiple permutation matrices per line of data, we need to keep track of the coefficients for each 

        for (size_t j = 0; j < datai.size() / 4; j++) {
            int Particlej = datai[4 * j] , Paulij = datai[4 * j + 1] , Powerj = datai[4 * j + 2] , twoSpinPlusOnej = datai[4 * j + 3];
            if (Particlej > NumOfParticles)
                NumOfParticles = Particlej;
            if (Paulij == 1) {
                // If Pauli X, then we have a P and a P^{-1} = P^{2s}. Appending both into the bitnum vector
                Permsi.insert(Permsi.end() , {1 , twoSpinPlusOnej - 1});
                // Setting the diagonals D^+ and D^- for P and P^- respectively
                Diagsi.insert(Diagsi.end() , {Dplus , Dminus});
                // The coefficient for each term is +1
                coeffsi.insert(coeffsi.end() , {one , one});
            } else if (Paulij == 2) {
                Permsi.insert(Permsi.end() , {1 , twoSpinPlusOnej - 1});
                Diagsi.insert(Diagsi.end() , {Dplus , Dminus});
                // The coefficient for the Pauli y is -i and i for D^+ and D^- respectively.
                coeffsi.insert(coeffsi.end() , {minusi , plusi});
            } else if (pauli_j == 3) {
                Permsi.push_back(0);
                Diagsi.push_back(Dz);
                coeffsi.push_back(one);
            }
        }
        coeffs.push_back(coeff_i);
        Zs.push_back(zs_i); // If zs is empty then no non-trivial diagonal components! (All identity operators)
            
        // Look for num in the previous list of Ps (permutations)
        pair<bool , int> bit_in_set = Bit_is_in_set(bit_num , Ps);
        if(!bit_in_set.first){ 
            // num was not found in Ps, thus a new permutation matrix!
            Ps.push_back(bit_num);
            // We will add i-th element of coeffs and Zs to be associated with the current Ps!
            Z_track.push_back(vector<int> {i});
        }
        else{
            // num was found in Ps, and P_index will be the index of Ps that matched num.
            int P_index = bit_in_set.second;

            vector<int> z_indices = Z_track[P_index];

            bool z_found = false;
            for (int k = 0; k < z_indices.size(); k++){
                if (Zs[z_indices[k]] == zs_i){
                    coeffs[P_index] += coeff_i;
                    z_found = true;
                    break;
                }
            }
            // If the z array is new, we add it to the Z_track for the Ps associated with num! 
            if (!z_found){
                Z_track[P_index].push_back(i);
            }
        }
    }

    // Throw away the zero coefficients:
    vector<bitset<5000>> Ps_kept;
    ZVecs Z_track_kept;

    for(int k = 0; k < Z_track.size(); k++){
        vector<int> ztrack_k;
        for(int l = 0; l < Z_track[k].size(); l++){
            int k_l_index = Z_track[k][l];
            complex<double> coeff_k_l = coeffs[k_l_index];
            if (abs(coeff_k_l) > 1e-8){
                //zero_coeffs_for_k.push_back(l);
                ztrack_k.push_back(k_l_index);
            }
        }
        if(ztrack_k.size() > 0){
            Z_track_kept.push_back(ztrack_k);
            Ps_kept.push_back(Ps[k]);
        }
    }

    // Sorting everything based on indices of Ps:
    vector<int> indices;
    for(int i = 0; i < Ps.size(); i++){
        indices.push_back(i);
    }
    sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
        return BitsetComparator()(Ps_kept[a], Ps_kept[b]);
    });


    vector<bitset<5000>> Ps_sorted;
    ZVecs Z_track_sorted;
    for(int i = 0; i < Ps_kept.size(); i++){
        Ps_sorted.push_back(Ps_kept[indices[i]]);
        Z_track_sorted.push_back(Z_track_kept[indices[i]]);
    }

    PZ_data.coeffs = coeffs; //Coeffs and Zs are kept as the original data
    PZ_data.Ps = Ps_sorted;
    PZ_data.Zs = Zs;
    PZ_data.Z_track = Z_track_sorted;

    return PZ_data;
}