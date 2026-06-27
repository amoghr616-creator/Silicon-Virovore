// engine.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "virovore.h" // Linked to your structural definitions

int amino_acid_charge[256] = {
    ['A'] = 0,  ['B'] = 0,  ['C'] = 0,  ['D'] = -1, ['E'] = -1,
    ['F'] = 0,  ['G'] = 0,  ['H'] = 0,  ['I'] = 0,  ['J'] = 0,
    ['K'] = 1,  ['L'] = 0,  ['M'] = 0,  ['N'] = 0,  ['O'] = 0,
    ['P'] = 0,  ['Q'] = 0,  ['R'] = 1,  ['S'] = 0,  ['T'] = 0,
    ['U'] = 0,  ['V'] = 0,  ['W'] = 0,  ['X'] = 0,  ['Y'] = 0, ['Z'] = 0
};

float calculate_charge_penalty(int net_charge, int threshold)
{
    int abs_charge = (net_charge < 0) ? -net_charge : net_charge;
    if (abs_charge > threshold)
    {
        return 0.0; 
    }
    else
    {
        return (float)(threshold + 1 - abs_charge);
    }
}

int compute_charge_profile(const char *sequence, int charge_threshold, charge_summary *summary_out)
{
    if (sequence == NULL || summary_out == NULL)
    {
        fprintf(stderr, "ERROR: NULL pointer in compute_charge_profile\n");
        return -1;
    }
    
    int num_positive = 0;
    int num_negative = 0;
    
    for (int i = 0; sequence[i] != '\0'; i++)
    {
        char aa = toupper((unsigned char)sequence[i]);
        
        if (aa == 'K' || aa == 'R')
        {
            num_positive++;
        }
        else if (aa == 'D' || aa == 'E')
        {
            num_negative++;
        }
        else if (aa < 'A' || aa > 'Z')
        {
            fprintf(stderr, "ERROR: Invalid amino acid '%c' at position %d\n", sequence[i], i);
            return -1;
        }
    }
    
    summary_out->num_positive = num_positive;
    summary_out->num_negative = num_negative;
    summary_out->net_charge = num_positive - num_negative;
    summary_out->penalty = calculate_charge_penalty(summary_out->net_charge, charge_threshold);
    
    return 0;  
}