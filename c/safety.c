#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "virovore.h"

// Forward declaration from engine.c
extern void evaluate_fitness(Variant *v, const char *target, const char *decoy1, const char *decoy2, void *params);

#define HERV_K_ENV "MKLAVDALLVTFAGSSDKKRR"
#define ALBUMIN "DAHKSEVAHRFKDLGEENFKALVL"
#define NCAM1 "MLQTKDLIWTLFFLGTAVS"

#define SEQ_LEN 21
static const char AA_ALPHABET[] = "ACDEFGHIKLMNPQRSTVWY";

static char random_amino_acid_different_from(char parental, char current) {
    char allowed[20];
    int count = 0;
    for (int i = 0; i < 20; i++) {
        char candidate = AA_ALPHABET[i];
        if (candidate == current) continue;
        /* If this site was already mutated, don't silently revert to parent. */
        if (current != parental && candidate == parental) continue;
        allowed[count++] = candidate;
    }
    if (count == 0) return current;
    return allowed[rand() % count];
}

static int is_valid_amino_acid(char c) {
    return (strchr(AA_ALPHABET, toupper((unsigned char)c)) != NULL);
}

/**
 * Executes full biophysical evaluation from engine.c on the given candidate.
 */
int c_evaluate_sequence(
    const char *sequence,
    FitnessMetrics *out
) {
    if (sequence == NULL || out == NULL) {
        return -1;
    }

    int len = strlen(sequence);
    if (len < 9) {
        return -2;
    }

    for (int i = 0; i < len; i++) {
        if (!is_valid_amino_acid(sequence[i])) {
            return -3;
        }
    }

    Variant v = init_variant(len + 1);
    strcpy(v.sequence, sequence);
    v.length = len;

    /*
     * Target alignment is currently not a biologically valid
     * peptide-Env interaction metric. Keep it neutral until a
     * proper target-specific objective is implemented.
     */
    evaluate_fitness(&v, NULL, ALBUMIN, NCAM1, NULL);

    out->fitness_score = v.fitness_score;
    out->solvation_energy = v.solvation_energy;
    out->hydrophobic_moment = v.hydrophobic_moment;
    out->helix_propensity = v.helix_propensity;
    out->target_alignment = v.target_alignment;
    out->decoy_penalty = v.decoy_penalty;
    out->charge_penalty = v.charge_penalty;
    out->charge_density = v.charge_density;

    free(v.sequence);
    return 0;
}

double c_check_sequence_fitness(const char *sequence) {
    FitnessMetrics metrics;

    if (c_evaluate_sequence(sequence, &metrics) != 0) {
        return -1.0;
    }

    return metrics.fitness_score;
}

void c_seed_random(unsigned int seed) {
    srand(seed);
}

/**
 * C dynamic population generator
 */
void c_generate_mutated_population(const char* seed_sequence, char output_population[][SEQ_LEN + 1], int pop_size, double mutation_rate) {
    if (seed_sequence == NULL || strlen(seed_sequence) != SEQ_LEN) return;

    for (int p = 0; p < pop_size; p++) {
        for (int i = 0; i < SEQ_LEN; i++) {
            double r = (double)rand() / RAND_MAX;
            if (r < mutation_rate) {
                output_population[p][i] = AA_ALPHABET[rand() % 20];
            } else {
                output_population[p][i] = seed_sequence[i];
            }
        }
        output_population[p][SEQ_LEN] = '\0';
    }
}

void c_generate_adaptive_population(
    const char *seed,
    double importance[],
    char output[][SEQ_LEN + 1],
    int pop_size)
{
    if (seed == NULL || strlen(seed) != SEQ_LEN) return;

    for (int p = 0; p < pop_size; p++) {
        strcpy(output[p], seed);

        for (int i = 0; i < SEQ_LEN; i++) {
            /*
             * High importance
             * -> mutate less
             *
             * Low importance
             * -> mutate more
             */

            double mutation_rate = 0.05 + (1.0 - importance[i]) * 0.20;
            double r = (double)rand() / RAND_MAX;

            if (r < mutation_rate) {
                output[p][i] = AA_ALPHABET[rand() % 20];
            }
        }

        output[p][SEQ_LEN] = '\0';
    }
}

/*
 * Generate a population with an identical per-candidate mutation-event
 * budget across uniform and hotspot-guided modes. Mutation positions are
 * sampled without replacement so repeated draws cannot silently reduce the
 * number of distinct changed positions.
 * guidance_mode: 0=uniform, 1=hotspot-only, 2/3=probability-biased hotspot.
 */
void c_generate_policy_population(
    const char *seed,
    double mutation_rate,
    int hotspot_start,
    int hotspot_end,
    int guidance_mode,
    double hotspot_probability,
    char output[][SEQ_LEN + 1],
    int pop_size)
{
    if (seed == NULL || strlen(seed) != SEQ_LEN || pop_size <= 0) return;

    int hotspot_valid = (
        hotspot_start >= 0
        && hotspot_start < hotspot_end
        && hotspot_end <= SEQ_LEN
    );
    if (!hotspot_valid || guidance_mode == 0) {
        guidance_mode = 0;
    }
    if (hotspot_probability < 0.0) hotspot_probability = 0.0;
    if (hotspot_probability > 1.0) hotspot_probability = 1.0;

    int mutation_events = (int)(mutation_rate * SEQ_LEN + 0.5);
    if (mutation_rate > 0.0 && mutation_events < 1) mutation_events = 1;

    for (int p = 0; p < pop_size; p++) {
        strcpy(output[p], seed);
        int used_positions[SEQ_LEN] = {0};
        for (int event = 0; event < mutation_events; event++) {
            int position = -1;
            int use_hotspot = 0;
            if (guidance_mode == 1) {
                use_hotspot = 1;
            } else if ((guidance_mode == 2 || guidance_mode == 3)
                       && (double)rand() / RAND_MAX < hotspot_probability) {
                use_hotspot = 1;
            }

            /* Select without replacement so the event budget equals the number
             * of distinct mutation positions whenever valid positions exist. */
            if (use_hotspot && hotspot_valid) {
                int available = 0;
                for (int i = hotspot_start; i < hotspot_end; i++) {
                    if (!used_positions[i]) available++;
                }
                if (available > 0) {
                    int pick = rand() % available;
                    for (int i = hotspot_start; i < hotspot_end; i++) {
                        if (used_positions[i]) continue;
                        if (pick-- == 0) { position = i; break; }
                    }
                } else if (guidance_mode == 1) {
                    /* hotspot_only stays inside the declared region even if
                     * the mutation budget exceeds its width. */
                    position = hotspot_start + rand() % (hotspot_end - hotspot_start);
                }
            }
            if (position < 0) {
                int available = 0;
                for (int i = 0; i < SEQ_LEN; i++) {
                    if (!used_positions[i]) available++;
                }
                if (available <= 0) break;
                int pick = rand() % available;
                for (int i = 0; i < SEQ_LEN; i++) {
                    if (used_positions[i]) continue;
                    if (pick-- == 0) { position = i; break; }
                }
            }
            if (position < 0) break;
            used_positions[position] = 1;
            output[p][position] = random_amino_acid_different_from(seed[position], output[p][position]);
        }
        output[p][SEQ_LEN] = '\0';
    }
}
