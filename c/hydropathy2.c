/*
 * hydropathy.c
 * Silicon Virovore Project — June 2026 Focus
 *
 * Streamlined O(N) linear sliding window hydropathy matrix calculator.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "virovore.h"

// Clear, explicit designated initializers matching textbook KD values
const float KD_SCALE[256] = {
    ['A'] = 1.8,  ['R'] = -4.5, ['N'] = -3.5, ['D'] = -3.5,
    ['C'] = 2.5,  ['Q'] = -3.5, ['E'] = -3.5, ['G'] = -0.4,
    ['H'] = -3.2, ['I'] = 4.5,  ['L'] = 3.8,  ['K'] = -3.9,
    ['M'] = 1.9,  ['F'] = 2.8,  ['P'] = -1.6, ['S'] = -0.8,
    ['T'] = -0.7, ['W'] = -0.9, ['Y'] = -1.3, ['V'] = 4.2
};

int scan_hydropathy_profile(const char *sequence, int window_size, hydropathy_summary *summary_out)
{
    if (sequence == NULL || summary_out == NULL) {
        return -1;
    }

    int len = strlen(sequence);
    
    if (len < window_size) {
        summary_out->max_peak = -99.0;
        summary_out->peak_detected = 0;
        summary_out->total_windows_scanned = 0;
        summary_out->solvation_penalty = 50.0; 
        return 0;
    }

    double current_window_sum = 0.0;
    float max_peak = -99.0;
    int windows_count = 0;

    // 1. Establish initial window baseline sum
    for (int i = 0; i < window_size; i++) {
        char aa = toupper((unsigned char)sequence[i]);
        current_window_sum += KD_SCALE[(int)aa];
    }
    max_peak = (float)(current_window_sum / window_size);
    windows_count++;

    // 2. Linear O(N) Sliding Loop: Drop the trailing edge residue, add the leading edge residue
    for (int i = window_size; i < len; i++) {
        char leaving_aa  = toupper((unsigned char)sequence[i - window_size]);
        char entering_aa = toupper((unsigned char)sequence[i]);

        current_window_sum += KD_SCALE[(int)entering_aa] - KD_SCALE[(int)leaving_aa];
        float window_avg = (float)(current_window_sum / window_size);
        windows_count++;

        if (window_avg > max_peak) {
            max_peak = window_avg;
        }
    }

    // 3. Export data profiles back to pipeline state tracking mechanisms
    summary_out->max_peak = max_peak;
    summary_out->total_windows_scanned = windows_count;
    summary_out->peak_detected = (max_peak > 1.6f) ? 1 : 0;

    if (summary_out->peak_detected) {
        summary_out->solvation_penalty = 0.0f;
    } else {
        summary_out->solvation_penalty = (1.6f - max_peak) * 5.0f;
    }

    return 0;
}