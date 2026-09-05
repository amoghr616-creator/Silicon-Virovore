// virovore.h
#ifndef VIROVORE_H
#define VIROVORE_H

#include <stddef.h>

#ifndef SEQ_LEN
#define SEQ_LEN 21
#endif

// Core Structures upgraded for High-Fidelity Biophysics
typedef struct {
    char *sequence;
    size_t length;
    size_t capacity; 
    double fitness_score;
    double solvation_energy;    // Delta G solvation parameter
    double hydrophobic_moment;   // Eisenberg structural vector magnitude
    double helix_propensity;    // Chou-Fasman folding stability score
    double target_alignment;    // Alignment score with target sequence
    double decoy_penalty;       // Exponential penalty for decoy binding
    double charge_penalty;      // Charge-based penalty
    double charge_density;      // Maximum localized charge density
} Variant;

typedef struct {
    Variant *variants;
    size_t count;
    size_t capacity;
} Population;

typedef struct {
    int net_charge;
    float penalty;      
    int num_positive;  
    int num_negative;
    double max_patch_density;   // Track localized coulombic patches
} charge_summary;

typedef struct {
    float max_peak;
    int peak_detected;         
    int total_windows_scanned;
    float solvation_penalty;   
} hydropathy_summary;

// Lookup Table Sizes
#define CHOU_FASMAN_SIZE 256

// Function Prototypes
Variant init_variant(size_t initial_capacity);
void append_sequence(Variant *v, const char *chunk);
double calculate_alignment_score(const char *variant, const char *profile);
void c_generate_adaptive_population(
    const char *seed,
    double importance[],
    char output[][SEQ_LEN + 1],
    int pop_size
);


// Advanced Physics Math Engines
int compute_charge_profile(const char *sequence, int charge_threshold, charge_summary *summary_out);
float calculate_charge_penalty(int net_charge, int threshold);
int scan_hydropathy_profile(const char *sequence, int window_size, hydropathy_summary *summary_out);

// New Biophysical Math Engines for July
double calculate_chou_fasman(const char *sequence);
double calculate_eisenberg_moment(const char *sequence);
double calculate_solvation_energy(const char *sequence);

#endif