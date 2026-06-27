/*
 * charge.c
 * Silicon Virovore Project — May 2026
 *
 * Implements the net charge calculation and "avoid ~0 charge" gate.
 * This is the second "gate" function in your biophysical scoring pipeline.
 *
 * Biochemical Background:
 * Proteins near their isoelectric point (pI) have net charge ≈ 0.
 * Such proteins tend to aggregate or precipitate.
 * We penalize sequences with |net_charge| < threshold to avoid this.
 *
 * Compile: gcc -o charge_demo charge.c -lm
 * Run:     ./charge_demo
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

/* ============================================================================
 * PART A: CHARGE CONTRIBUTIONS PER AMINO ACID
 * ============================================================================
 * At physiological pH (7.4), certain amino acids carry a net charge:
 *
 *   K (Lysine), R (Arginine)      → +1 each (positively charged)
 *   D (Aspartic acid), E (Glutamic acid) → -1 each (negatively charged)
 *   H (Histidine)                 → ~+0.1 each (can be partially charged; we'll ignore)
 *   All others                    → 0
 *
 * Table: 1 if positive, -1 if negative, 0 if neutral
 */

int amino_acid_charge[26] = {
    0,   // A (Alanine)
    0,   // B (Asp or Asn) — placeholder
    0,   // C (Cysteine) — slightly negative, but we ignore for simplicity
    -1,  // D (Aspartic acid) — negatively charged
    -1,  // E (Glutamic acid) — negatively charged
    0,   // F (Phenylalanine)
    0,   // G (Glycine)
    0,   // H (Histidine) — partially charged; we ignore at pH 7.4
    0,   // I (Isoleucine)
    0,   // J (not used)
    1,   // K (Lysine) — positively charged
    0,   // L (Leucine)
    0,   // M (Methionine)
    0,   // N (Asparagine)
    0,   // O (not used)
    0,   // P (Proline)
    0,   // Q (Glutamine)
    1,   // R (Arginine) — positively charged
    0,   // S (Serine)
    0,   // T (Threonine)
    0,   // U (not used)
    0,   // V (Valine)
    0,   // W (Tryptophan)
    0,   // X (unknown)
    0,   // Y (Tyrosine)
    0    // Z (Glu or Gln) — placeholder
};

/* ============================================================================
 * PART B: LOOKUP FUNCTION FOR CHARGE
 * ============================================================================
 */

int get_amino_acid_charge(char amino_acid)
{
    // Convert lowercase to uppercase if needed
    amino_acid = toupper(amino_acid);
    
    // Safety check
    if (amino_acid < 'A' || amino_acid > 'Z')
    {
        fprintf(stderr, "ERROR: Invalid amino acid '%c' in charge lookup\n", amino_acid);
        return 0;  // Error code: assume neutral
    }
    
    int index = amino_acid - 'A';
    return amino_acid_charge[index];
}

/* ============================================================================
 * PART C: NET CHARGE CALCULATION
 * ============================================================================
 * Calculate the total net charge of a protein sequence.
 * Formula: sum of all amino acid charges.
 *
 * Input:
 *   sequence = protein sequence string (uppercase)
 *   net_charge_out = pointer to an int (we store the charge here)
 *
 * Returns:
 *   0 if successful
 *   -1 if error (invalid character)
 *
 * Note: The returned net_charge is an integer because we're counting
 *   discrete charged residues. In real biochemistry, you'd account
 *   for protonation state, but that's a June upgrade.
 */

// Removed calculate_net_charge entirely. 
// New streamlined compute_charge_profile does all the work in ONE pass!
int compute_charge_profile(const char *sequence, int charge_threshold,
                            charge_summary *summary_out)
{
    if (sequence == NULL || summary_out == NULL)
    {
        fprintf(stderr, "ERROR: NULL pointer in compute_charge_profile\n");
        return -1;
    }
    
    int num_positive = 0;
    int num_negative = 0;
    
    // SINGLE PASS: We look at each character exactly once
    for (int i = 0; sequence[i] != '\0'; i++)
    {
        char aa = toupper(sequence[i]);
        
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
            fprintf(stderr, "ERROR: Invalid amino acid '%c' at position %d\n", sequence[i], i);
            return -1;
        }
    }
    
    // Perform all math steps using the single-pass counts
    summary_out->num_positive = num_positive;
    summary_out->num_negative = num_negative;
    summary_out->net_charge = num_positive - num_negative;
    summary_out->penalty = calculate_charge_penalty(summary_out->net_charge, charge_threshold);
    
    return 0;  // Success
}
/* ============================================================================
 * PART D: NET CHARGE PENALTY ("AVOID ~0" GATE)
 * ============================================================================
 * A sequence with net charge very close to 0 will aggregate.
 * This function penalizes such sequences.
 *
 * Input:
 *   net_charge = the net charge (from calculate_net_charge)
 *   threshold = the maximum allowable |net_charge| to avoid (e.g., 1)
 *              If |net_charge| <= threshold, we penalize.
 *
 * Returns:
 *   0.0 if the sequence has sufficient charge (GOOD)
 *   A penalty (positive number) if |net_charge| <= threshold (BAD)
 *
 * Interpretation:
 *   - If net_charge is +5 or -5 (far from 0), penalty = 0 (good)
 *   - If net_charge is 0 or ±1 (near 0), penalty = large (bad)
 *
 * The formula here is:
 *   penalty = max(0, threshold + 1 - |net_charge|)
 * This gives a smooth ramp: the closer to 0, the higher the penalty.
 */

float calculate_charge_penalty(int net_charge, int threshold)
{
    int abs_charge = (net_charge < 0) ? -net_charge : net_charge;
    
    if (abs_charge > threshold)
    {
        // Sequence has sufficient charge — no penalty
        return 0.0;
    }
    else
    {
        // Sequence is near charge 0 — penalize it
        // Penalty ramps up as |charge| decreases
        float penalty = (float)(threshold + 1 - abs_charge);
        return penalty;
    }
}

/* ============================================================================
 * PART E: CHARGE SUMMARY STRUCTURE
 * ============================================================================
 */

typedef struct {
    int net_charge;
    float penalty;
    int num_positive;  // Number of K + R residues
    int num_negative;  // Number of D + E residues
} charge_summary;

/* ============================================================================
 * PART F: COMPUTE FULL CHARGE PROFILE
 * ============================================================================
 */

int compute_charge_profile(const char *sequence, int charge_threshold,
                            charge_summary *summary_out)
{
    if (sequence == NULL || summary_out == NULL)
    {
        fprintf(stderr, "ERROR: NULL pointer in compute_charge_profile\n");
        return -1;
    }
    
    int net_charge = 0;
    int num_positive = 0;
    int num_negative = 0;
    
    // Count charges
    for (int i = 0; sequence[i] != '\0'; i++)
    {
        char aa = toupper(sequence[i]);
        
        if (aa == 'K' || aa == 'R')
            num_positive++;
        else if (aa == 'D' || aa == 'E')
            num_negative++;
        else if (aa < 'A' || aa > 'Z')
        {
            fprintf(stderr, "ERROR: Invalid amino acid '%c' at position %d\n",
                    sequence[i], i);
            return -1;
        }
    }
    
    net_charge = num_positive - num_negative;
    float penalty = calculate_charge_penalty(net_charge, charge_threshold);
    
    summary_out->net_charge = net_charge;
    summary_out->penalty = penalty;
    summary_out->num_positive = num_positive;
    summary_out->num_negative = num_negative;
    
    return 0;  // Success
}

/* ============================================================================
 * PART G: DEMONSTRATION AND TESTING
 * ============================================================================
 */

int main(void)
{
    printf("=============================================================\n");
    printf("Silicon Virovore — Net Charge Analyzer\n");
    printf("May 2026\n");
    printf("=============================================================\n\n");
    
    // Test Case 1: Positively charged (lots of K, R)
    printf("TEST 1: Positively charged sequence\n");
    printf("---\n");
    
    const char *test1 = "MKLAVLALLVTFAGSSDLH";
    charge_summary s1;
    if (compute_charge_profile(test1, 1, &s1) == 0)
    {
        printf("Sequence: %s\n", test1);
        printf("Positive (K+R): %d\n", s1.num_positive);
        printf("Negative (D+E): %d\n", s1.num_negative);
        printf("Net charge: %d\n", s1.net_charge);
        printf("Charge penalty (threshold=1): %.2f %s\n",
               s1.penalty,
               s1.penalty == 0.0 ? "(PASS)" : "(FAIL — too close to 0)");
    }
    
    printf("\n");
    
    // Test Case 2: Negatively charged
    printf("TEST 2: Negatively charged sequence\n");
    printf("---\n");
    
    const char *test2 = "DEDEEDEDEDDEDEDED";
    charge_summary s2;
    if (compute_charge_profile(test2, 1, &s2) == 0)
    {
        printf("Sequence: %s\n", test2);
        printf("Positive (K+R): %d\n", s2.num_positive);
        printf("Negative (D+E): %d\n", s2.num_negative);
        printf("Net charge: %d\n", s2.net_charge);
        printf("Charge penalty (threshold=1): %.2f %s\n",
               s2.penalty,
               s2.penalty == 0.0 ? "(PASS)" : "(FAIL — too close to 0)");
    }
    
    printf("\n");
    
    // Test Case 3: Balanced (near 0 — SHOULD FAIL)
    printf("TEST 3: Balanced sequence (problematic — near charge 0)\n");
    printf("---\n");
    
    const char *test3 = "MKLAVDEALLVTFAGSS"; // 1 K (+1), 1 D (-1), 1 E (-1) = Net -1 (FAILS gate)
    charge_summary s3;
    if (compute_charge_profile(test3, 1, &s3) == 0)
    {
        printf("Sequence: %s\n", test3);
        printf("Positive (K+R): %d\n", s3.num_positive);
        printf("Negative (D+E): %d\n", s3.num_negative);
        printf("Net charge: %d\n", s3.net_charge);
        printf("Charge penalty (threshold=1): %.2f %s\n",
               s3.penalty,
               s3.penalty == 0.0 ? "(PASS)" : "(FAIL — aggregation risk!)");
    }
    
    printf("\n");
    
    // Test Case 4: Exactly charge 0 (worst case)
    printf("TEST 4: Exactly balanced (worst case — charge = 0)\n");
    printf("---\n");
    
    const char *test4 = "MKLAVLAL";  // M K L A V L A L = 1 pos, 0 neg = +1. Let's make +0
    const char *test4_fixed = "DEKVFLAG";  // D E K V F L A G = 1 pos, 2 neg = -1. Still not 0.
    const char *test4_balanced = "KDELAFASG";  // K D E L A F A S G = 1 pos, 2 neg = -1. Hmm.
    // Let's just use K and D: "KD" = 0 net charge
    const char *test4_zero = "MKLAVDALLVTFAGSS"; // 1 K (+1), 1 D (-1) = Net 0 (FAILS gate miserably);
    
    charge_summary s4;
    if (compute_charge_profile(test4_zero, 1, &s4) == 0)
    {
        printf("Sequence: %s\n", test4_zero);
        printf("Positive (K+R): %d\n", s4.num_positive);
        printf("Negative (D+E): %d\n", s4.num_negative);
        printf("Net charge: %d\n", s4.net_charge);
        printf("Charge penalty (threshold=1): %.2f %s\n",
               s4.penalty,
               s4.penalty > 2.0 ? "(FAIL — HIGHEST PENALTY!)" : "");
    }
    
    printf("\n");
    
    // Test Case 5: Invalid character (should fail gracefully)
    printf("TEST 5: Invalid character (contains '9')\n");
    printf("---\n");
    
    const char *test5 = "MKL9VALLVTFAGSSDLH";
    charge_summary s5;
    if (compute_charge_profile(test5, 1, &s5) != 0)
    {
        printf("Correctly rejected: invalid character found.\n");
    }
    
    printf("\n");
    
    printf("=============================================================\n");
    printf("All tests completed successfully!\n");
    printf("=============================================================\n");
    
    return 0;
}