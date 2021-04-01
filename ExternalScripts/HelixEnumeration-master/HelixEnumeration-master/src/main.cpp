#include "Enumerate.h"
#include <iostream>
#include <string>
#include <vector>
#include <unistd.h>

using namespace Enumerate;

void PrintHelp()
{
    std::cout << "HelixEnumeration version x.x.x" << std::endl << std::endl;
    
    std::cout << "\t-h: Print this help message and exit" << std::endl;
    std::cout << "\t-s [sequence]: Provide the sequence as a series of capital letters. ";
    std::cout << "If this argument is not set then the program will look for the ";
    std::cout << "sequence through cin." << std::endl;
    std::cout << "\t-k [int]: The minimum length of a helix class to be included in ";
    std::cout << "the enumeration. Defaults to 1." << std::endl;
    std::cout << "\t-l [int]: The length of a hairpin bend. This is the minimum gap ";
    std::cout << "between sides of a helix class. Defaults to 3." << std::endl;

    std::cout << std::endl;

    std::cout << "Usage: ./HelixEnumeration -s AAAAUUUUGGGGCCCC" << std::endl;
    std::cout << "Usage: echo AAAAUUUUGGGGCCCC | ./HelixEnumeration" << std::endl;
    std::cout << "Usage: cat sequence.txt | ./HelixEnumeration" << std::endl;

    std::cout << std::endl;
}

void ReadArguments(
        const int argc, 
        char **argv,
        std::string& sequence,
        int& min_k,
        int& hairpin_length)
{
    sequence = "";
    min_k = 1;
    hairpin_length = 3;

    int arg_letter;
    while((arg_letter = getopt(argc,argv,"hs:k:l:")) != -1) 
    {
        switch(arg_letter)
        {
        case 'h':
            PrintHelp();
            exit(0);
        case 's':
            sequence = optarg;
            break;
        case 'k':
            min_k = std::stoi(optarg);
            break;
        case 'l':
            hairpin_length = std::stoi(optarg);
            break;
        case '?':

            break;
        }
    }

    if (sequence == "")
    {
        std::getline(std::cin,sequence);
    }
}

int main(int argc, char **argv) 
{
    std::string sequence;
    int min_k, hairpin_length;

    ReadArguments(
        argc, argv,
        sequence, min_k, hairpin_length);
    
    std::vector<HelixClass> classes = EnumerateHelices(
                                            sequence,
                                            min_k,
                                            hairpin_length);

    for (std::size_t idx = 0; idx < classes.size(); idx++)
    {
        std::cout << classes.at(idx).i << " " << 
            classes.at(idx).j << " " <<
            classes.at(idx).k << std::endl;
    }

    return 0;
}

