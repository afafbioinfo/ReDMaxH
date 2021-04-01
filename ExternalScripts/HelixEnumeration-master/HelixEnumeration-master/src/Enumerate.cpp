#include "Enumerate.h"
#include <set>

using namespace Enumerate;

namespace {

const std::set<char>* getPairedLetters(char letter)
{
    static const std::set<char> a {'U'};
    static const std::set<char> u {'A','G'};
    static const std::set<char> g {'C','U'};
    static const std::set<char> c {'G'};
    static const std::set<char> blank {};

    switch (letter) {
        case 'A':
            return &a;
        case 'U':
            return &u;
        case 'G':
            return &g;
        case 'C':
            return &c;
        default:
            fprintf(stderr, "ERROR: sequence contained unrecognized character\n");
            return &blank;
    }
}

//calculates an approximation of -energy of any given base pair
int getBasePairEnergy(char a, char b)
{
    //flip the order of a and b to reduce the number of if statements necessary
    if (b < a)
    {
        const char temp = a;
        a = b;
        b = temp;
    }

    if (a == 'A')
    {
        if (b == 'U')
            return 2;
    }
    else if (a == 'C')
    {
        if (b == 'G')
            return 3;
    }
    else if (a == 'G')
    {
        if (b == 'U')
            return 1;
    }

    return 0;
}

//calculates an approximation of -energy for a helix class
int getHelixClassEnergy(int start, int end, int length, const std::string& sequence)
{
    int weight = 0;
    for (int i = 0; i < length; i++)
    {
        weight += getBasePairEnergy(
            sequence.at(start + i),
            sequence.at(end - i));
    }
    return weight;
}

} //anonymous namespace

std::vector<HelixClass> Enumerate::EnumerateHelices(
        const std::string& sequence, 
        const int min_k, 
        const int hairpin_length)
{

    //Weight, start, end, length vector
    std::vector<HelixClass> classes;
    std::vector<std::set<int> > skip_array(sequence.size());

    for (int idx = 0; idx < (int)sequence.size(); idx++)
    {
        const char letter = sequence.at(idx);
        const std::set<char>* paired_letters = getPairedLetters(letter);

        int end_idx = sequence.size() - 1;

        while (end_idx > idx + min_k * 2 - 1 + hairpin_length)
        {
            if (skip_array.at(idx).count(end_idx))
            {
                end_idx--;
                continue;
            }

            int length = 0;
            const std::set<char>* paired_temp = paired_letters;

            int moving_idx = end_idx;
            while (
                (idx + length + hairpin_length < moving_idx) &&
                (paired_temp->count(sequence.at(moving_idx))))
            {
                skip_array.at(idx + length).insert(moving_idx);

                length++;
                paired_temp = getPairedLetters(sequence.at(idx+length));
                moving_idx--;
            }

            const int end_pair = moving_idx + length;
            const int true_length = std::min(length,(int)(moving_idx+length-idx+1)/2);

            if (true_length >= min_k)
                classes.push_back(
                    HelixClass(idx,end_pair,true_length));

            end_idx--;
        }
    }

    return classes;
}
