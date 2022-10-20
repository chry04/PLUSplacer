#include <iostream>
#include <fstream>
#include <omp.h>


int main( int argc, char **argv ){
    if( argc <= 7 ){
        std::cerr << "Usage: "<<argv[0]<<" [ref infile] [nbr ref records] [query infile] [nbr query records] [outfile] [nbr of subtrees] [nbr of votes]" << std::endl;
        return -1;
    }

  
    // read in query sequences first
    std::ifstream input(argv[3]);
    if(!input.good()){
        std::cerr << "Error opening '"<<argv[3]<<"'. Bailing out." << std::endl;
        return -1;
    }
    std::string line, name, content;
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
    


    // read in the reference sequences second
    double smallest_total_hamming[count2];
    int best_tree_index[count2];
    for (int c=0; c<count2; c++)
    {
        smallest_total_hamming[c] = 999999999999999;
        best_tree_index[c] = 0;
    }
    std::cout << argv[2] << std::endl;
    int subtree_size = std::stoi(argv[2]);
    int nbr_subtrees = std::stoi(argv[6]);
    
    std::cout << nbr_subtrees << std::endl;
    
    
    // std::cout << "here!!" << std::endl;
    
    std::string subtree_path = argv[1];
    int count1 = 0;
    
    std::string name_arr[subtree_size+1];
    std::string seq_arr[subtree_size+1];
    
    for (int subtree_idx = 0; subtree_idx < nbr_subtrees; subtree_idx++){
    
        
        std::string subtree_full_path = subtree_path + std::to_string(subtree_idx);
        
        std::cout << subtree_full_path << std::endl;
        
        std::ifstream input_q(subtree_full_path); 

        if(!input_q.good()){
            std::cerr << "Error opening '"<< subtree_full_path <<"'. Bailing out." << std::endl;
            return -1;
        }

        count1 = 0;

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
        //std::cout << "finished all but last reference sequences" << std::endl;
    
        if( !name.empty() ){ // Print out what we read from the last entry
            //std::cout << name << " : " << content << std::endl;
            
            name_arr[count1] = name.c_str();                 
            seq_arr[count1] = content.c_str();
            count1++;
        }    
        std::cout << "ref count: "<< count1 <<" query count2: " <<count2 << std::endl;
        //std::cout << "finished reading reference sequences" << std::endl;
    
        // find n (size) closest reference sequences by Hamming distance to each query 
        // print this to outfile with the one query per line followed by the n closest
        // reference sequences separated by a comma, with their requisite Hamming
        // distances separated with a semicolon.
    
        #pragma omp parallel for 
        for (int c2=0; c2<count2; c2++){ //query seq array
    
        int size = std::stoi(argv[7]);
        int best_hamming[size];
        int best_index[size];
        int furthest_index = 0;

        for (int i=0; i<size; i++){
            best_index[i] = 0;
            best_hamming[i] = 999999999;
        }

            int q_len = q_seq_arr[c2].length();
            int start_idx = q_len;
            int end_idx = q_len;

            for (int i=0; i < q_len; i++) {
                if (q_seq_arr[c2][i] != '-' && start_idx >= q_len) {
                    start_idx = i;
                }
                if (q_seq_arr[c2][q_len-i-1] != '-' && end_idx >= q_len) {
                    end_idx = q_len - i - 1;
                }
            }


            for (int c1=0; c1<count1 ; c1++) { //ref seq array
                int count = 0;
                for(int i=start_idx; i < end_idx+1; i++) { //individual seq hamming distances
                    if(seq_arr[c1][i] != q_seq_arr[c2][i]) {
                        count++;
                        if (count > best_hamming[furthest_index]) {
                            break;
                        }
                    }
                }
                
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
            int count = 0;
            for (int i= 0; i<size; i++){
                count += best_hamming[i];
            }
            //std::cout << "here" << std::endl;
            if (count <= smallest_total_hamming[c2]) {
                smallest_total_hamming[c2] = count;
                best_tree_index[c2] = subtree_idx;
            }
        }
        std::cout << "finished hamming distance " << std::endl;
        
    }
    
    std::ofstream outFile(argv[5]);

    for (int c2=0; c2<count2; c2++){ //query seq array
        outFile << q_name_arr[c2] << "," << best_tree_index[c2] << ":" << smallest_total_hamming[c2] <<std::endl;
    }
    
    outFile.close();

    return 0;
}
