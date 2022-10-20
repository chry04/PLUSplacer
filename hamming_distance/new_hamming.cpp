#include <iostream>
#include <fstream>
#include <omp.h>


int main( int argc, char **argv ){
    if( argc <= 6 ){
        std::cerr << "Usage: "<<argv[0]<<" [ref infile] [nbr ref records] [query infile] [nbr query records] [outfile] [nbr of leaves returned]" << std::endl;
        return -1;
    }
 
    std::ifstream input_q(argv[1]);
    if(!input_q.good()){
        std::cerr << "Error opening '"<<argv[1]<<"'. Bailing out." << std::endl;
        return -1;
    }
    std::cout << argv[2] << std::endl;
    std::string name_arr[std::stoi(argv[2])];
    std::string seq_arr[std::stoi(argv[2])];
    std::string line, name, content;
    int count1 = 0;

    while( std::getline( input_q, line ).good() ){
        if( line.empty() || line[0] == '>' ){ // Identifier marker
            if( !name.empty() ){ // Print out what we read from the last entry
                //std::cout << name << " : " << content << std::endl;
                name_arr[count1] = name.c_str();
                seq_arr[count1] = content.c_str();
                name.clear();
                //std::cout << count1<< " : " << name_arr[count1] <<std::endl;
                count1++;
            }
            if( !line.empty() ){
                name = line.substr(1);
            }
            content.clear();
        } else if( !name.empty() ){
            if( line.find(' ') != std::string::npos ){ // Invalid sequence--no spaces allowed
                name.clear();
                content.clear();
            } else {
                content += line;
            }
        }

    }
    
    
    if( !name.empty() ){ // Print out what we read from the last entry
        //std::cout << name << " : " << content << std::endl;        
        name_arr[count1] = name;
        seq_arr[count1] = content;
        count1++;
    }
  
  
    std::ifstream input(argv[3]);
    if(!input.good()){
        std::cerr << "Error opening '"<<argv[3]<<"'. Bailing out." << std::endl;
        return -1;
    }
    std::cout << argv[4] << std::endl;
    std::string q_name_arr[std::stoi(argv[4])+3];
    std::string q_seq_arr[std::stoi(argv[4])+3];
    int count2 = 0;
    name = "";
    while( std::getline( input, line ).good() ){
        if( line.empty() || line[0] == '>' ){ // Identifier marker
            if( !name.empty() ){ // Print out what we read from the last entry
                //std::cout << name << " : " << content << std::endl;
                q_name_arr[count2] = name.c_str();
                q_seq_arr[count2] = content.c_str();
                name.clear();
                //std::cout << count2 << " : " << q_name_arr[count2] <<std::endl;
                count2++;
            }
            if( !line.empty() ){
                name = line.substr(1);
            }
            content.clear();
        } else if( !name.empty() ){
            if( line.find(' ') != std::string::npos ){ // Invalid sequence--no spaces allowed
                name.clear();
                content.clear();
            } else {
                content += line;
            }
        }

    }
    

    if( !name.empty() ){ // Print out what we read from the last entry
        //std::cout << name << " : " << content << std::endl;        
        q_name_arr[count2] = name;
        q_seq_arr[count2] = content;
        count2++;
    }

    std::cout << "ref count: "<< count1 <<" query count2: " <<count2 << std::endl;

    std::ofstream outFile(argv[5]);

    #pragma omp parallel for 
    for (int c2=0; c2<count2; c2++){ //query seq array
    
        int size = std::stoi(argv[6]);
        int best_hamming[size];
        int best_index[size];
        int furthest_index = 0;

        for (int i=0; i<size; i++){
            best_index[i] = 0;
            best_hamming[i] = 999999999;
        }

        for (int c1=0; c1<count1 ; c1++) { //ref seq array
            int count = 0;
            int len = seq_arr[c1].length();
            for(int i=0; i < len; i++) {
                if(seq_arr[c1][i] != q_seq_arr[c2][i]) {
                    count++;
                    if (count > best_hamming[furthest_index]) {
                        break;
                    }
                }
            }
            //std::cout << "here" << std::endl;
            if (count <= best_hamming[furthest_index]) {
                best_hamming[furthest_index] = count;
                best_index[furthest_index] = c1;
                int high_hamming = 0;
                int high_index = 0;
                for (int i=0; i<size; i++){
                    if (best_hamming[i] > high_hamming){
                        high_hamming = best_hamming[i];
                        high_index = i;
                    }
                furthest_index = high_index;
                }
            }
        }
        #pragma omp critical
        {
            outFile << q_name_arr[c2];
            for (int i=0; i<size; i++){
                outFile << "," << name_arr[best_index[i]] << ":" << best_hamming[i];
            }
            outFile << std::endl;
        }
    }
    outFile.close();
    return 0;
}
