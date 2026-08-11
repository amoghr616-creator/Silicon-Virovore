#include <stdio.h>
#include <string.h>
#include <ctype.h>

// Replace these declarations with your actual function signatures from grep
extern void init_variant(void); 
extern double calculate_alignment_score(const char* seq);

static int is_valid_amino_acid(char c) {
    const char *valid_aa = "ACDEFGHIKLMNPQRSTVWY";
    return (strchr(valid_aa, toupper((unsigned char)c)) != NULL);
}

double c_check_sequence_fitness(const char* sequence) {
    if (sequence == NULL) return -999.0;
    
    int len = strlen(sequence);
    if (len < 9) return -1.0;
    
    for (int i = 0; i < len; i++) {
        if (!is_valid_amino_acid(sequence[i])) {
            return -1.0;
        }
    }
    
    // 1. Initialize any required C global tables/state
    // init_variant(); 
    
    // 2. Safely call internal routines
    // double align = calculate_alignment_score(sequence);
    // return align;

    return 1.0; // Baseline safe score for pipeline testing
}
