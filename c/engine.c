// engine.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "virovore.h" 

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
    size_t len = strlen(sequence);
    double max_patch = 0.0;
    
    // Window-based Localized Charge Density Tracking (Debye-Hückel precursor)
    int window_size = 6;
    if (len >= (size_t)window_size) {
        for (size_t i = 0; i <= len - window_size; i++) {
            int local_charge = 0;
            for (int j = 0; j < window_size; j++) {
                char local_aa = toupper((unsigned char)sequence[i + j]);
                if (local_aa == 'K' || local_aa == 'R') local_charge++;
                if (local_aa == 'D' || local_aa == 'E') local_charge--;
            }
            double current_density = abs(local_charge) / (double)window_size;            
            if (current_density > max_patch) {
                max_patch = current_density;
            }
        }
    }
    
    for (size_t i = 0; i < len; i++)
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
            fprintf(stderr, "ERROR: Invalid amino acid '%c' at position %zu\n", sequence[i], i);
            return -1;
        }
    }
    
    summary_out->num_positive = num_positive;
    summary_out->num_negative = num_negative;
    summary_out->net_charge = num_positive - num_negative;
    summary_out->penalty = calculate_charge_penalty(summary_out->net_charge, charge_threshold);
    summary_out->max_patch_density = max_patch;

    return 0;  
}

// Chou-Fasman alpha-helix propensity lookup array for amino acids
double get_chou_fasman_score(char aa) {
    switch (toupper((unsigned char)aa)) {
        case 'E': case 'A': case 'L': return 1.42; 
        case 'H': case 'M': case 'Q': case 'W': case 'V': case 'F': return 1.10;
        case 'K': case 'I': return 1.00; 
        case 'D': case 'T': case 'S': case 'R': case 'C': return 0.80; 
        case 'N': case 'Y': return 0.60; 
        case 'G': case 'P': return 0.50; 
        default: return 0.00;
    }
}

double calculate_chou_fasman(const char *sequence) {
    if (!sequence || *sequence == '\0') return 0.0;
    
    double total_score = 0.0;
    size_t len = strlen(sequence);
    
    for (size_t i = 0; i < len; i++) {
        total_score += get_chou_fasman_score(sequence[i]);
    }
    
    return total_score / (double)len;
}

// Atomic Solvation Parameter lookup (Simplified SASA free energy in kcal/mol per residue)
double get_solvation_parameter(char aa) {
    switch (toupper((unsigned char)aa)) {
        case 'I': case 'L': case 'V': case 'F': return -2.5; 
        case 'W': case 'Y': case 'M': return -2.0;
        case 'A': case 'C': return -0.5;
        case 'G': case 'P': case 'S': case 'T': return 0.5;  
        case 'H': case 'Q': case 'N': return 1.5;
        case 'D': case 'E': case 'K': case 'R': return 3.0;  
        default: return 0.0;
    }
}

double calculate_solvation_energy(const char *sequence) {
    if (!sequence || *sequence == '\0') return 0.0;
    
    double delta_g = 0.0;
    size_t len = strlen(sequence);
    
    for (size_t i = 0; i < len; i++) {
        delta_g += get_solvation_parameter(sequence[i]);
    }
    
    return delta_g; 
}

/* ============================================================================
 * LINKING UTILITIES: Required by ga_loop.c initialization steps
 * ============================================================================
 */

// Dynamic memory allocation constructor for a new Variant sequence string
Variant init_variant(size_t initial_capacity) {
    Variant v;
    v.sequence = (char *)malloc(initial_capacity * sizeof(char));
    if (!v.sequence) {
        fprintf(stderr, "FATAL ERROR: Out of memory during init_variant allocation\n");
        exit(EXIT_FAILURE);
    }
    v.sequence[0] = '\0'; 
    v.length = 0;
    v.capacity = initial_capacity;
    v.fitness_score = 0.0;
    v.solvation_energy = 0.0;
    v.hydrophobic_moment = 0.0;
    v.helix_propensity = 0.0;
    return v;
}

// Appends string chunks dynamically, auto-resizing the block via realloc if needed
void append_sequence(Variant *v, const char *chunk) {
    if (!v || !chunk) return;
    
    size_t chunk_len = strlen(chunk);
    size_t required_capacity = v->length + chunk_len + 1; 
    
    if (required_capacity > v->capacity) {
        size_t new_capacity = v->capacity * 2;
        if (new_capacity < required_capacity) {
            new_capacity = required_capacity;
        }
        
        char *new_seq = (char *)realloc(v->sequence, new_capacity * sizeof(char));
        if (!new_seq) {
            fprintf(stderr, "FATAL ERROR: Out of memory during append_sequence realloc\n");
            exit(EXIT_FAILURE);
        }
        v->sequence = new_seq;
        v->capacity = new_capacity;
    }
    
    strcpy(v->sequence + v->length, chunk);
    v->length += chunk_len;
}