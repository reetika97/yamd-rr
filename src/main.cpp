//
// Created by raor on 19.07.24.
//
#include "header_files/xyz.h"
#include "header_files/functions.h"
#include <iostream>

int main(int argc, char* argv[]) {

    std::string pgm_selection = "Default";
    if(argc==1) {
        std::cout << "Select_program: " && std::cin >> pgm_selection;
        std::cout << pgm_selection << " is being executed." << std::endl;
    }
    else if(argc==2 ){
        pgm_selection = argv[1];
        std::cout <<pgm_selection<< " is being executed." << std::endl;
    }
    if(pgm_selection=="Energy_conservation"){

        Energy_conservation();

    }
    else {
        std::cout << "Invalid argument(s)" << std::endl;
        exit(0);
    }

    return 0;
}