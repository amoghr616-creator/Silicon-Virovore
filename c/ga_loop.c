/*
 * ga_loop.c
 * Silicon Virovore Project — May 2026
 *
 * Implements a minimal, seeded Genetic Algorithm (GA) loop run.
 * Evolves a sequence across generations to minimize total biophysical penalties.
 *
 * Compile: gcc -o ga_demo ga_loop.c -lm
 * Run:     ./ga_demo
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>

const char AA_POOL[] = "ACDEFGHIKLMNPQRSTVWY";

typedef struct {
    int net_charge;
    float charge_penalty;
} charge_summary;

typedef struct {
    float max_peak;
    int peak_detected;
    float solvation_penalty;
} hydropathy_summary;

typedef struct {
    float total_fitness_penalty;
    int passed_all_gates;
    charge_summary charge;
    hydropathy_summary hydro;
} pipeline_score;

int amino_acid_charge[26] = {
    0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0
};

float amino_acid_hydropathy[26] = {
    1.8,  0.0,  2.5, -3.5, -3.5,  2.8, -0.4, -3.2,  4.5,  0.0, 
   -3.9,  3.8,  1.9, -3.5,  0.0, -1.6, -3.5, -4.5, -0.8, -0.7, 
    0.0,  4.2, -0.9,  0.0, -1.3,  0.0
};

float get_amino_acid_hydropathy(char amino_acid) {
    amino_acid = toupper(amino_acid);
    if (amino_acid < 'A' || amino_acid > 'Z') return 0.0;
    return amino_acid_hydropathy[amino_acid - 'A'];
}

float calculate_charge_penalty(int net_charge, int threshold) {
    int abs_charge = (net_charge < 0) ? -net_charge : net_charge;
    if (abs_charge > threshold) return 0.0;
    return (float)(threshold + 1 - abs_charge);
}

int evaluate_sequence(const char *sequence, int charge_threshold, pipeline_score *eval_out) {
    int len = strlen(sequence);
    
    // Charge calculation
    int num_positive = 0, num_negative = 0;
    for (int i = 0; i < len; i++) {
        char aa = toupper(sequence[i]);
        if (aa == 'K' || aa == 'R') num_positive++;
        else if (aa == 'D' || aa == 'E') num_negative++;
    }
    eval_out->charge.net_charge = num_positive - num_negative;
    eval_out->charge.charge_penalty = calculate_charge_penalty(eval_out->charge.net_charge, charge_threshold);

    // Hydropathy calculation
    float max_peak = -99.0;
    for (int i = 0; i <= len - 19; i++) {
        float window_sum = 0.0;
        for (int j = 0; j < 19; j++) {
            window_sum += get_amino_acid_hydropathy(sequence[i + j]);
        }
        float current_avg = window_sum / 19.0;
        if (current_avg > max_peak) max_peak = current_avg;
    }
    eval_out->hydro.max_peak = max_peak;
    eval_out->hydro.peak_detected = (max_peak > 1.6) ? 1 : 0;
    eval_out->hydro.solvation_penalty = eval_out->hydro.peak_detected ? 0.0f : (1.6f - max_peak) * 5.0f;

    // Combined Score
    eval_out->total_fitness_penalty = eval_out->charge.charge_penalty + eval_out->hydro.solvation_penalty;
    eval_out->passed_all_gates = (eval_out->charge.charge_penalty == 0.0f && eval_out->hydro.peak_detected) ? 1 : 0;
    
    return 0;
}

void mutate_sequence(char *sequence) {
    int len = strlen(sequence);
    int target_index = rand() % len;
    int random_aa_index = rand() % 20;
    sequence[target_index] = AA_POOL[random_aa_index];
}

int main(void) {
    // Seed random number generator with a fixed value for reproducible demo runs
    srand(42); 

    printf("=============================================================\n");
    printf("Silicon Virovore — Minimal Evolutionary GA Loop Run\n");
    printf("May 2026\n");
    printf("=============================================================\n\n");

    // Seed sequence: 25 residues, balanced near charge 0 and highly soluble (Terrible starting fitness)
    char current_sequence[64] = "MKLAVDALLVTFAGSSDKKRRKKRR"; 
    char candidate_sequence[64];
    
    pipeline_score current_score;
    evaluate_sequence(current_sequence, 1, &current_score);

    printf("INITIAL SEED CANDIDATE:\n");
    printf("Sequence: %s\n", current_sequence);
    printf("Initial Penalty Score: %.2f\n", current_score.total_fitness_penalty);
    printf("-------------------------------------------------------------\n");
    printf("Starting optimization iterations...\n\n");

    int generation = 0;
    int max_generations = 2000;

    while (generation < max_generations && current_score.total_fitness_penalty > 0.0f) {
        generation++;

        // Duplicate current sequence to modify
        strcpy(candidate_sequence, current_sequence);
        mutate_sequence(candidate_sequence);

        pipeline_score candidate_score;
        evaluate_sequence(candidate_sequence, 1, &candidate_score);

        // Natural Selection: Accept mutation only if it decreases the penalty score
        if (candidate_score.total_fitness_penalty < current_score.total_fitness_penalty) {
            strcpy(current_sequence, candidate_sequence);
            current_score = candidate_score;
            
            printf("Gen %4d | New Fittest Penalty: %5.2f | Sequence: %s\n", 
                   generation, current_score.total_fitness_penalty, current_sequence);
        }
    }

    printf("\n-------------------------------------------------------------\n");
    printf("EVOLUTIONARY SEARCH TERMINATED AT GENERATION: %d\n", generation);
    printf("Final Optimized Sequence: %s\n", current_sequence);
    printf("Final Net Charge        : %d\n", current_score.charge.net_charge);
    printf("Final Hydropathy Peak   : %.2f\n", current_score.hydro.max_peak);
    printf("FINAL PIPELINE STATUS   : %s\n", current_score.passed_all_gates ? "PASS (Zero Penalty Achieved)" : "FAIL");
    printf("=============================================================\n");

    return 0;
}