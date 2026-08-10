/*
 * ga_loop.c
 * Silicon Virovore Project — 10,000 Generation Master Sync Engine
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "virovore.h"

// --- LOCALISED BACKUP MACROS (Ensures compilation if header paths are unsynced) ---
#ifndef HERV_K_ENV
#define HERV_K_ENV "MKLAVDALLVTFAGSSDKKRR"
#endif
#ifndef ALBUMIN
#define ALBUMIN "DAHKSEVAHRFKDLGEENFKALVL"
#endif
#ifndef NCAM1
#define NCAM1 "MLQTKDLIWTLFFLGTAVS"
#endif

// Forward declare external biophysical fitness evaluator
void evaluate_fitness(Variant *v, const char *target, const char *decoy1, const char *decoy2, void *params);

// Standard 20 natural L-amino acids pool
const char AA_POOL[] = "ACDEFGHIKLMNPQRSTVWY";

// Thread-safe quicksort comparison utility for rankings
int compare_variants(const void *a, const void *b) {
    double fit_a = ((const Variant *)a)->fitness_score;
    double fit_b = ((const Variant *)b)->fitness_score;
    if (fit_a > fit_b) return -1;
    if (fit_a < fit_b) return 1;
    return 0;
}

// In-place mutation with structural boundaries
void mutate_variant(Variant *v, double rate) {
    if (!v || !v->sequence) return;
    for (size_t i = 0; i < v->length; i++) {
        double r = (double)rand() / RAND_MAX;
        if (r < rate) {
            v->sequence[i] = AA_POOL[rand() % 20];
        }
    }
}

// Single-point crossover producing safe structural progeny
void crossover(const Variant *p1, const Variant *p2, Variant *child) {
    if (!p1 || !p2 || !child) return;
    
    size_t min_len = (p1->length < p2->length) ? p1->length : p2->length;
    if (min_len <= 1) {
        strncpy(child->sequence, p1->sequence, child->capacity - 1);
        child->sequence[child->capacity - 1] = '\0';
        child->length = strlen(child->sequence);
        return;
    }

    size_t split = 1 + (rand() % (min_len - 1));
    
    // Copy first half from parent 1
    memcpy(child->sequence, p1->sequence, split);
    
    // Copy second half from parent 2
    size_t second_half_len = p2->length - split;
    if (split + second_half_len >= child->capacity) {
        second_half_len = child->capacity - split - 1;
    }
    memcpy(child->sequence + split, p2->sequence + split, second_half_len);
    
    child->sequence[split + second_half_len] = '\0';
    child->length = strlen(child->sequence);
}

// High-performance tournament selection to avoid genetic drift
int tournament_selection(const Population *pop, int tournament_size) {
    int pop_limit = (int)pop->count; 
    if (pop_limit <= 0) {
        pop_limit = tournament_size * 2; // Dynamic safety fallback
    }

    int best_idx = rand() % pop_limit;
    double best_fitness = pop->variants[best_idx].fitness_score;

    for (int i = 1; i < tournament_size; i++) {
        int idx = rand() % pop_limit;
        if (pop->variants[idx].fitness_score > best_fitness) {
            best_idx = idx;
            best_fitness = pop->variants[idx].fitness_score;
        }
    }
    return best_idx;
}

// -------------------------------------------------------------------------
// MASTER EVOLUTIONARY LOOP EXECUTION PIPELINE
// -------------------------------------------------------------------------
void run_evolutionary_loop(int generations, int pop_size, double mutation_rate, int tournament_size) {
    printf("[*] Instantiating Silicon Virovore population vectors...\n");
    printf("    - Population Size: %d | Target Generations: %d\n", pop_size, generations);
    
    // Seed random number generator
    srand((unsigned int)time(NULL));

    // Initialize population
    Population pop;
    pop.count = (size_t)pop_size;
    pop.capacity = (size_t)pop_size;
    pop.variants = (Variant *)malloc(pop_size * sizeof(Variant));
    if (!pop.variants) {
        fprintf(stderr, "FATAL ERROR: Failed to allocate memory for population array.\n");
        exit(EXIT_FAILURE);
    }

    // Assign initial variant states safely
    for (int i = 0; i < pop_size; i++) {
        pop.variants[i] = init_variant(32); // Pre-allocate 32 bytes to prevent fragmentation
        for (size_t j = 0; j < 21; j++) {
            pop.variants[i].sequence[j] = AA_POOL[rand() % 20];
        }
        pop.variants[i].sequence[21] = '\0';
        pop.variants[i].length = 21;
    }

    // Open high-performance CSV streams for visualization dashboard mapping
    FILE *metrics_file = fopen("metrics.csv", "w");
    FILE *heatmap_file = fopen("evolution_consensus.csv", "w");
    if (!metrics_file || !heatmap_file) {
        fprintf(stderr, "FATAL ERROR: Unable to construct downstream CSV file sync handles.\n");
        exit(EXIT_FAILURE);
    }

    // Print headers for mathematical modeling scripts
    fprintf(metrics_file, "Generation,FitnessScore,SolvationEnergy,Eisenberg,HelixPropensity\n");
    fprintf(heatmap_file, "Generation,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y\n");

    // Pre-allocated structural buffer for offspring creation
    Variant child = init_variant(32);

    // 10,000 Generation Optimization Routine
    for (int gen = 0; gen < generations; gen++) {
        
        // --- STEP 1: DYNAMIC ADAPTIVE MUTATION RATE ---
        // Decays the exploration rate linearly from the user-defined base rate to 0.5%
        double progress = (double)gen / (double)generations;
        double dynamic_mutation_rate = mutation_rate * (1.0 - progress) + 0.005 * progress;

        // Step 2: Calculate fitness for all active variants
        for (int i = 0; i < pop_size; i++) {
            evaluate_fitness(&pop.variants[i], HERV_K_ENV, ALBUMIN, NCAM1, NULL);
        }

        // Step 3: Sort population in-place (Elitist selection)
        qsort(pop.variants, pop_size, sizeof(Variant), compare_variants);

        // Track mathematical consensus values for the heatmap dashboard
        int aa_counts[256] = {0};
        for (int i = 0; i < pop_size; i++) {
            for (size_t j = 0; j < pop.variants[i].length; j++) {
                aa_counts[(unsigned char)pop.variants[i].sequence[j]]++;
            }
        }

        // Export metrics for current generation
        fprintf(metrics_file, "%d,%.4f,%.4f,%.4f,%.4f\n", 
                gen, 
                pop.variants[0].fitness_score, 
                pop.variants[0].solvation_energy, 
                pop.variants[0].hydrophobic_moment, 
                pop.variants[0].helix_propensity);

        // Export consensus distribution profile
        fprintf(heatmap_file, "%d", gen);
        for (int a = 0; a < 20; a++) {
            fprintf(heatmap_file, ",%d", aa_counts[(unsigned char)AA_POOL[a]]);
        }
        fprintf(heatmap_file, "\n");

        // Print telemetry summary every 1000 generations
        if (gen % 1000 == 0 || gen == generations - 1) {
            printf("[Gen %d] Current Champion: %s | Fitness: %.4f | ΔG Solv: %.2f | Mut-Rate: %.3f%%\n", 
                   gen, pop.variants[0].sequence, pop.variants[0].fitness_score, pop.variants[0].solvation_energy, dynamic_mutation_rate * 100.0);
        }

        // Step 4: Breed next generation with safe pointer reallocations
        Population next_pop;
        next_pop.count = (size_t)pop_size;
        next_pop.capacity = (size_t)pop_size;
        next_pop.variants = (Variant *)malloc(pop_size * sizeof(Variant));
        if (!next_pop.variants) {
            fprintf(stderr, "FATAL ERROR: Memory allocation failure during generation breeding.\n");
            exit(EXIT_FAILURE);
        }

        // Preserve top 10% elite candidates unchanged (Elitism)
        int elite_count = pop_size / 10;
        for (int i = 0; i < elite_count; i++) {
            next_pop.variants[i] = init_variant(32);
            strcpy(next_pop.variants[i].sequence, pop.variants[i].sequence);
            next_pop.variants[i].length = pop.variants[i].length;
            next_pop.variants[i].fitness_score = pop.variants[i].fitness_score;
            next_pop.variants[i].solvation_energy = pop.variants[i].solvation_energy;
            next_pop.variants[i].hydrophobic_moment = pop.variants[i].hydrophobic_moment;
            next_pop.variants[i].helix_propensity = pop.variants[i].helix_propensity;
        }

        // Fill remaining population with bred offspring using adaptive mutation
        for (int i = elite_count; i < pop_size; i++) {
            int p1_idx = tournament_selection(&pop, tournament_size);
            int p2_idx = tournament_selection(&pop, tournament_size);

            crossover(&pop.variants[p1_idx], &pop.variants[p2_idx], &child);
            mutate_variant(&child, dynamic_mutation_rate);

            next_pop.variants[i] = init_variant(32);
            strcpy(next_pop.variants[i].sequence, child.sequence);
            next_pop.variants[i].length = child.length;
        }

        // Deep free previous population allocation safely
        for (int i = 0; i < pop_size; i++) {
            free(pop.variants[i].sequence);
        }
        free(pop.variants);

        // Advance to next generation
        pop = next_pop;
    }

    // Final file output flushes
    fclose(heatmap_file);
    fclose(metrics_file);
    printf("[SUCCESS] File sync pipelines complete. Simulation data exported.\n");

    // Sort final generation to isolate the true biological winner
    qsort(pop.variants, pop_size, sizeof(Variant), compare_variants);

    // ------------------------------------------------------------
    // TERMINAL CHAMPION REPORT (The Final Verdict)
    // ------------------------------------------------------------
    printf("\n=============================================================\n");
    printf("DESIGN COMPLETE: THE WINNING SILICON VIROVORE PEPTIDE\n");
    printf("=============================================================\n");
    printf("Final Peak Fitness: %.4f\n", pop.variants[0].fitness_score);
    printf("Sequence (1-21):    %s\n", pop.variants[0].sequence);
    printf("α-Helix Propensity: %.4f\n", pop.variants[0].helix_propensity);
    printf("Hydrophobic Moment: %.4f\n", pop.variants[0].hydrophobic_moment);
    printf("=============================================================\n\n");

    // Clean up temporary child buffer
    free(child.sequence);

    // Final memory sanitation step (Leak profile: 0 bytes)
    for (int i = 0; i < pop_size; i++) {
        free(pop.variants[i].sequence);
    }
    free(pop.variants);
    printf("[CLEANUP] Memory sanitation complete. System stabilized.\n");
}