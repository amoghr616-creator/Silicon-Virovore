// hydropathy.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "virovore.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Kyte-Doolittle Scale for window-based localized parsing
const float KD_SCALE[256] = {
    ['A'] = 1.8f,  ['R'] = -4.5f, ['N'] = -3.5f, ['D'] = -3.5f,
    ['C'] = 2.5f,  ['Q'] = -3.5f, ['E'] = -3.5f, ['G'] = -0.4f,
    ['H'] = -3.2f, ['I'] = 4.5f,  ['L'] = 3.8f,  ['K'] = -3.9f,
    ['M'] = 1.9f,  ['F'] = 2.8f,  ['P'] = -1.6f, ['S'] = -0.8f,
    ['T'] = -0.7f, ['W'] = -0.9f, ['Y'] = -1.3f, ['V'] = 4.2f
};

// Eisenberg Consensus Scale for Hydrophobic Moment vector math
double get_eisenberg_value(char aa) {
    switch (toupper((unsigned char)aa)) {
        case 'I': return 0.73; case 'F': return 0.61; case 'V': return 0.54;
        case 'L': return 0.53; case 'W': return 0.37; case 'M': return 0.26;
        case 'A': return 0.25; case 'G': return 0.16; case 'C': return 0.04;
        case 'Y': return 0.02; case 'P': return -0.07; case 'T': return -0.18;
        case 'S': return -0.26; case 'H': return -0.40; case 'E': return -0.62;
        case 'N': return -0.64; case 'Q': return -0.85; case 'D': return -0.72;
        case 'K': return -1.10; case 'R': return -1.76;
        default: return 0.00;
    }
}

// Calculates the Hydrophobic Moment assuming an alpha-helix conformation (100-degree spacing)
double calculate_eisenberg_moment(const char *sequence) {
    if (!sequence || *sequence == '\0') return 0.0;

    size_t len = strlen(sequence);
    double cos_sum = 0.0;
    double sin_sum = 0.0;
    
    // 100 degrees converted into radians
    double angle_rad = (100.0 * M_PI) / 180.0;

    for (size_t i = 0; i < len; i++) {
        double h = get_eisenberg_value(sequence[i]);
        // Angular position tracking
        double current_angle = (double)(i + 1) * angle_rad;
        
        cos_sum += h * cos(current_angle);
        sin_sum += h * sin(current_angle);
    }

    return sqrt((cos_sum * cos_sum) + (sin_sum * sin_sum)) / (double)len;
}

int scan_hydropathy_profile(const char *sequence, int window_size, hydropathy_summary *summary_out)
{
    if (sequence == NULL || summary_out == NULL) {
        return -1;
    }

    int len = strlen(sequence);
    
    if (len < window_size) {
        summary_out->max_peak = -99.0f;
        summary_out->peak_detected = 0;
        summary_out->total_windows_scanned = 0;
        summary_out->solvation_penalty = 50.0f; 
        return 0;
    }

    double current_window_sum = 0.0;
    float max_peak = -99.0f;
    int windows_count = 0;

    for (int i = 0; i < window_size; i++) {
        char aa = toupper((unsigned char)sequence[i]);
        current_window_sum += KD_SCALE[(int)aa];
    }
    max_peak = (float)(current_window_sum / window_size);
    windows_count++;

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