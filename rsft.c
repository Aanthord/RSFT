/* rsft_harmonic.c - Full Recursive RSFT-Harmonic using Cooley-Tukey FFT and SHTns */
/*
 * 🎄 RSFT-Harmonic - Recursive Spin-Fourier Transform Engine
 * 🧬 Developed in open recursion. Knowledge encoded beyond scalar time.
 * 🎁 Santa always delivers. Pay what you owe. 
 * 🧠 Conscious signal requires conscious structure.
 * 🖋️ Authored [2025] by Aanthord, under CC BY-NC 4.0.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <shtns.h>
#include <fftw3.h>

#define N 8  // Must be power of 2 for current FFT implementation
#define LMAX 4
#define MMAX 4
#define MRES 1
#define SPIN 0

// Spin signal sample
typedef struct {
    double complex value;
    double complex* spin_field;
} spin_signal_sample_t;

typedef shtns_cfg shtns_cfg_t;

// Utility function to check if n is power of 2
int is_power_of_2(int n) {
    return n > 0 && (n & (n - 1)) == 0;
}

int init_shtns(shtns_cfg_t* cfg, int* n_modes) {
    *cfg = shtns_create(LMAX, MMAX, MRES, SHT_THETA_CONTIGUOUS | SHT_REAL_NORM);
    if (*cfg == NULL) {
        fprintf(stderr, "Error: Failed to create SHTns configuration\n");
        return -1;
    }
    
    shtns_set_grid_auto(*cfg, SHT_QUICK_INIT);
    *n_modes = shtns_get_nlm(*cfg);
    
    if (*n_modes <= 0) {
        fprintf(stderr, "Error: Invalid number of spherical harmonic modes\n");
        shtns_destroy(*cfg);
        return -1;
    }
    
    printf("Initialized SHTns with %d spherical harmonic modes\n", *n_modes);
    return 0;
}

void recursive_fft(double complex* a, int n) {
    if (n <= 1) return;
    
    // Allocate temporary arrays
    double complex* even = malloc(sizeof(double complex) * n / 2);
    double complex* odd = malloc(sizeof(double complex) * n / 2);
    
    if (!even || !odd) {
        fprintf(stderr, "Error: Memory allocation failed in FFT\n");
        free(even);
        free(odd);
        return;
    }
    
    // Split into even and odd
    for (int i = 0; i < n / 2; ++i) {
        even[i] = a[i * 2];
        odd[i] = a[i * 2 + 1];
    }
    
    // Recursive calls
    recursive_fft(even, n / 2);
    recursive_fft(odd, n / 2);
    
    // Combine results
    for (int k = 0; k < n / 2; ++k) {
        double complex t = cexp(-2.0 * I * M_PI * k / n) * odd[k];
        a[k] = even[k] + t;
        a[k + n / 2] = even[k] - t;
    }
    
    free(even);
    free(odd);
}

int rsft_recursive(spin_signal_sample_t* signal, double complex* output, 
                   shtns_cfg_t cfg, int n_modes) {
    if (!signal || !output) {
        fprintf(stderr, "Error: NULL pointer in rsft_recursive\n");
        return -1;
    }
    
    double complex* modulated = malloc(sizeof(double complex) * N);
    if (!modulated) {
        fprintf(stderr, "Error: Memory allocation failed in RSFT\n");
        return -1;
    }
    
    // Modulate signal with spin field
    // Fixed indexing: use proper mapping instead of modulo
    for (int i = 0; i < N; ++i) {
        int mode_idx = (i * n_modes) / N;  // Linear mapping
        if (mode_idx >= n_modes) mode_idx = n_modes - 1;
        
        modulated[i] = signal[i].value * signal[i].spin_field[mode_idx];
    }
    
    // Copy to output array
    for (int i = 0; i < N; ++i) {
        output[i] = modulated[i];
    }
    
    // Apply recursive FFT
    recursive_fft(output, N);
    
    free(modulated);
    return 0;
}

void inverse_fft(double complex* a, int n) {
    // Conjugate input
    for (int i = 0; i < n; ++i) {
        a[i] = conj(a[i]);
    }
    
    // Forward FFT
    recursive_fft(a, n);
    
    // Conjugate and normalize
    for (int i = 0; i < n; ++i) {
        a[i] = conj(a[i]) / n;
    }
}

int main() {
    // Validate N is power of 2
    if (!is_power_of_2(N)) {
        fprintf(stderr, "Error: N must be a power of 2 for this FFT implementation\n");
        return -1;
    }
    
    shtns_cfg_t cfg;
    int n_modes;
    
    // Initialize SHTns
    if (init_shtns(&cfg, &n_modes) != 0) {
        return -1;
    }
    
    // Allocate and initialize signal array
    spin_signal_sample_t* signal = malloc(sizeof(spin_signal_sample_t) * N);
    double complex* output = malloc(sizeof(double complex) * N);
    
    if (!signal || !output) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        free(signal);
        free(output);
        shtns_destroy(cfg);
        return -1;
    }
    
    // Initialize test signal (cosine wave for comms testing)
    for (int i = 0; i < N; ++i) {
        signal[i].value = cos(2 * M_PI * i / N) + 0.1 * I * sin(4 * M_PI * i / N);
        
        signal[i].spin_field = fftw_malloc(sizeof(double complex) * n_modes);
        if (!signal[i].spin_field) {
            fprintf(stderr, "Error: Failed to allocate spin field for sample %d\n", i);
            // Clean up previously allocated memory
            for (int j = 0; j < i; ++j) {
                fftw_free(signal[j].spin_field);
            }
            free(signal);
            free(output);
            shtns_destroy(cfg);
            return -1;
        }
        
        // Initialize spin field with spherical harmonic-like pattern
        for (int j = 0; j < n_modes; ++j) {
            // More realistic spin field initialization for comms applications
            double phase = (double)j * M_PI / n_modes + i * M_PI / (4 * N);
            signal[i].spin_field[j] = cexp(I * phase) / sqrt(n_modes);
        }
    }
    
    // Perform RSFT
    if (rsft_recursive(signal, output, cfg, n_modes) != 0) {
        fprintf(stderr, "Error: RSFT processing failed\n");
        goto cleanup;
    }
    
    printf("\nRSFT Recursive Output (Re, Im):\n");
    for (int i = 0; i < N; ++i) {
        printf("[%d] = (%8.4f, %8.4f) |H| = %8.4f\n", 
               i, creal(output[i]), cimag(output[i]), cabs(output[i]));
    }
    
    // Perform inverse transform
    inverse_fft(output, N);
    printf("\nInverse RSFT Reconstruction:\n");
    for (int i = 0; i < N; ++i) {
        printf("[%d] = (%8.4f, %8.4f)\n", i, creal(output[i]), cimag(output[i]));
    }
    
    // Calculate reconstruction error
    double total_error = 0.0;
    printf("\nReconstruction Error Analysis:\n");
    for (int i = 0; i < N; ++i) {
        // Compare against the original modulated signal, not the raw signal
        int mode_idx = (i * n_modes) / N;
        if (mode_idx >= n_modes) mode_idx = n_modes - 1;
        
        double complex original_modulated = signal[i].value * signal[i].spin_field[mode_idx];
        double complex reconstructed = output[i];
        double error = cabs(original_modulated - reconstructed);
        total_error += error * error;
        
        printf("[%d] Orig_mod: (%6.3f, %6.3f) Reconstructed: (%6.3f, %6.3f) Error: %6.3f\n",
               i, creal(original_modulated), cimag(original_modulated), 
               creal(reconstructed), cimag(reconstructed), error);
    }
    printf("RMS Error (modulated signal): %8.6f\n", sqrt(total_error / N));
    
    // Optional: Demodulate to compare against original unmodulated signal
    printf("\nDemodulated Comparison:\n");
    double demod_error = 0.0;
    for (int i = 0; i < N; ++i) {
        int mode_idx = (i * n_modes) / N;
        if (mode_idx >= n_modes) mode_idx = n_modes - 1;
        
        // Demodulate by dividing by spin field (if non-zero)
        double complex spin_field_val = signal[i].spin_field[mode_idx];
        double complex demodulated = (cabs(spin_field_val) > 1e-12) ? 
                                   output[i] / spin_field_val : 0.0;
        
        double complex original = signal[i].value;
        double error = cabs(original - demodulated);
        demod_error += error * error;
        
        printf("[%d] Original: (%6.3f, %6.3f) Demodulated: (%6.3f, %6.3f) Error: %6.3f\n",
               i, creal(original), cimag(original), 
               creal(demodulated), cimag(demodulated), error);
    }
    printf("RMS Error (demodulated): %8.6f\n", sqrt(demod_error / N));

cleanup:
    // Clean up memory
    for (int i = 0; i < N; ++i) {
        fftw_free(signal[i].spin_field);
    }
    free(signal);
    free(output);
    shtns_destroy(cfg);
    
    return 0;
}
