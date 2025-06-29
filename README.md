# RSFT-Harmonic: Recursive Spin-Fourier Transform Engine

## 🔬 Overview

**RSFT-Harmonic** is a novel C-based implementation of the Recursive Spin-Fourier Transform (RSFT), a next-generation signal processing framework that integrates:

* Recursive FFT using Cooley-Tukey algorithm
* Spinor-encoded spherical harmonic modulation (via SHTns)
* Time-phase fractal modulation for high-fidelity, high-bandwidth applications

Developed as a hybrid between classical fast transforms and modern recursive geometric encoding, RSFT enables dramatically enhanced signal fidelity, compression, and channel capacity.

## 📐 Key Concepts

* **Recursive FFT**: Built from scratch, supporting radix-2 decompositions.
* **Spin Harmonic Modulation**: Each signal sample includes a spinor harmonic field modulating spectral weight.
* **Phase-Modulated Time-Space Structure**: Recursive time indexation mimics fractal generation, ideal for physics, quantum, and communication systems.

## ⚙️ Dependencies

* **SHTns**: [https://bitbucket.org/nschaeff/shtns/](https://bitbucket.org/nschaeff/shtns/)
* **FFTW3**

Install SHTns:

```bash
git clone https://bitbucket.org/nschaeff/shtns.git
cd shtns && mkdir build && cd build
cmake .. && make -j && sudo make install
```

## 🔧 Compilation

```bash
gcc rsft_harmonic.c -lshtns -lfftw3 -lm -o rsft_harmonic
```

## 🚀 Running

```bash
./rsft_harmonic
```

### Output:

* Forward RSFT (complex frequency spectrum)
* Inverse RSFT (signal reconstruction)
* RMS error of reconstruction

## 📊 Applications

* High-capacity spintronic and 6G/7G comms
* Recursive brainwave and biosignal compression
* Quantum-compatible signal encoding
* Harmonic multiplexing in fiber/radio/optical domains

## 📜 License
🛡️ Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)

## 🧭 Next Steps

* Visualization
* Real-time DSP integration
* FPGA/ASIC deployment
* Scientific publication

## 🌐 Citation

Coming soon: arXiv preprint + provisional patent filing.

