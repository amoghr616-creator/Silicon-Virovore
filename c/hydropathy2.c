/*
 * hydropathy.c
 * Silicon Virovore Project — May 2026
 *
 * Implements the 19-residue Kyte-Doolittle hydropathy window scan 
 * and peak detection (> 1.6 threshold) for membrane-binding hook validation.
 * Includes a negative design heuristic penalty for soluble sequences.
 *
 * Compile: gcc -o hydro_demo hydropathy.c
 * Run:     ./hydro_demo
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

typedef struct {
    float max_peak;
    int peak_detected;         // 1 if max_peak > 1.6, otherwise 0
    int total_windows_scanned;
    float solvation_penalty;   // Distance penalty from target threshold (0.0 if passed)
} hydropathy_summary;

// Lookup array containing standard textbook values from the Kyte-Doolittle scale
float amino_acid_hydropathy[26] = {
    1.8,  0.0,  2.5, -3.5, -3.5,  2.8, -0.4, -3.2,  4.5,  0.0, 
   -3.9,  3.8,  1.9, -3.5,  0.0, -1.6, -3.5, -4.5, -0.8, -0.7, 
    0.0,  4.2, -0.9,  0.0, -1.3,  0.0
};

float get_amino_acid_hydropathy(char amino_acid)
{
    amino_acid = toupper(amino_acid);
    if (amino_acid < 'A' || amino_acid > 'Z') {
        return 0.0; // Fail-safe fallback for unexpected characters
    }
    return amino_acid_hydropathy[amino_acid - 'A'];
}

int scan_hydropathy_profile(const char *sequence, hydropathy_summary *summary_out)
{
    if (sequence == NULL || summary_out == NULL) {
        return -1;
    }

    int len = strlen(sequence);
    
    // Safety check: A sequence shorter than the window size cannot be processed
    if (len < 19) {
        fprintf(stderr, "WARNING: Sequence length (%d) too short for 19-residue window.\n", len);
        summary_out->max_peak = -99.0;
        summary_out->peak_detected = 0;
        summary_out->total_windows_scanned = 0;
        summary_out->solvation_penalty = 50.0; // Fixed high penalty for invalid inputs
        return 0;
    }

    float max_peak = -99.0;
    int windows_count = 0;

    /* ============================================================================
     * SLIDING WINDOW CORE ENGINE
     * ============================================================================
     * Outer loop 'i' represents the starting character of the current window.
     * It terminates at (len - 19) to prevent reading past the end of the string.
     */
    for (int i = 0; i <= len - 19; i++) 
    {
        float window_sum = 0.0;
        
        // Inner loop 'j' gathers the 19 characters inside the current window
        for (int j = 0; j < 19; j++) 
        {
            char aa = sequence[i + j];
            window_sum += get_amino_acid_hydropathy(aa);
        }
        
        float current_window_average = window_sum / 19.0;
        windows_count++;
        
        // Track the single highest hydrophobicity score encountered
        if (current_window_average > max_peak) 
        {
            max_peak = current_window_average;
        }
    }

    // Export properties to the structural pointer output
    summary_out->max_peak = max_peak;
    summary_out->total_windows_scanned = windows_count;
    summary_out->peak_detected = (max_peak > 1.6) ? 1 : 0;

    /* ============================================================================
     * NEGATIVE DESIGN: SOLVATION & AGGREGATION HEURISTIC PENALTY
     * ============================================================================
     */
    if (summary_out->peak_detected) 
    {
        summary_out->solvation_penalty = 0.0; // Meets criteria, zero penalty applied
    } 
    else 
    {
        // Computes how far below the 1.6 threshold the sequence is and scales by 5.0
        summary_out->solvation_penalty = (1.6f - max_peak) * 5.0f;
    }

    return 0;
}

int main(void)
{
    printf("=============================================================\n");
    printf("Silicon Virovore — Hydropathy Scanner Prototype\n");
    printf("May 2026\n");
    printf("=============================================================\n\n");

    // TEST 1: Highly hydrophobic simulated segment (Should Pass)
    printf("TEST 1: Transmembrane Hydrophobic Segment\n");
    printf("---\n");
    const char *membrane_test = "IVLIVLIVLIVLIVLIVLIVL"; 
    hydropathy_summary r1;

    if (scan_hydropathy_profile(membrane_test, &r1) == 0) {
        printf("Sequence: %s\n", membrane_test);
        printf("Windows Scanned: %d\n", r1.total_windows_scanned);
        printf("Highest Window Average: %.2f\n", r1.max_peak);
        printf("Membrane Hook Status: %s\n", r1.peak_detected ? "PASS (Valid Transmembrane Hook)" : "FAIL");
        printf("Solvation Penalty: %.2f\n\n", r1.solvation_penalty);
    }

    // TEST 2: Highly soluble cytosolic protein segment (Should Fail)
    printf("TEST 2: Highly Soluble Hydrophilic Segment\n");
    printf("---\n");
    const char *soluble_test = "KKRRKKRRKKRRKKRRKKRR"; 
    hydropathy_summary r2;

    if (scan_hydropathy_profile(soluble_test, &r2) == 0) {
        printf("Sequence: %s\n", soluble_test);
        printf("Windows Scanned: %d\n", r2.total_windows_scanned);
        printf("Highest Window Average: %.2f\n", r2.max_peak);
        printf("Membrane Hook Status: %s\n", r2.peak_detected ? "PASS" : "FAIL (Too soluble to bind membrane)");
        printf("Solvation Penalty: %.2f\n\n", r2.solvation_penalty);
    }

    printf("=============================================================\n");
    printf("All sliding window loops executed successfully!\n");
    printf("=============================================================\n");

    return 0;
}