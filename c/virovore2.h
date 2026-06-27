// virovore.h
#ifndef VIROVORE_H
#define VIROVORE_H

#include <stddef.h>

// Population Dynamic Structures
typedef struct {
    char *sequence;
    size_t length;
    size_t capacity; 
    double fitness_score;
} Variant;

typedef struct {
    Variant *variants;
    size_t count;
    size_t capacity;
} Population;

// Engine Gate 2 Structures: Net Charge Profile
typedef struct {
    int net_charge;
    float penalty;
    int num_positive;  
    int num_negative;  
} charge_summary;

// Function Declarations
Variant init_variant(size_t initial_capacity);
void append_sequence(Variant *v, const char *chunk);
double calculate_max_hydropathy(const char *sequence, int window_size);
int compute_charge_profile(const char *sequence, int charge_threshold, charge_summary *summary_out);
float calculate_charge_penalty(int net_charge, int threshold);

// virovore.h snippet
typedef struct {
    float max_peak;
    int peak_detected;         
    int total_windows_scanned;
    float solvation_penalty;   
} hydropathy_summary;

// Declare the optimized function signature
int scan_hydropathy_profile(const char *sequence, int window_size, hydropathy_summary *summary_out);

#endif