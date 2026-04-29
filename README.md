# PicoDst_Ee_spectra_analyzer
C++/ROOT analyzer for STAR PicoDst data: event/track QA, electron PID, e+e- invariant mass spectra, conversion rejection, and mixed-event background subtraction.



The main project file is:

picoAnalyzer.cpp


## Features

- reads STAR PicoDst input with `StPicoDstReader`;
- applies basic event selection using the reconstructed primary vertex;
- applies primary-track quality cuts;
- performs electron PID using TPC `nSigmaElectron` and TOF `1 / beta` information;
- fills event, track, PID, and detector QA histograms;
- suppresses conversion-electron candidates with a low-mass opposite-sign pair cut;
- builds the same-event invariant-mass spectrum of `e+e-` pairs;
- estimates the combinatorial background using mixed events;
- normalizes and subtracts the mixed-event background;
- stores all output histograms in a ROOT file.

## Requirements

To build and run the analyzer, you need:

- a C++17-compatible compiler;
- ROOT 6;
- the STAR PicoDst environment;
- `libStPicoDst`;
- PicoDst headers, including:
 - `StPicoDstReader.h`;
 - `StPicoDst.h`;
 - `StPicoEvent.h`;
 - `StPicoTrack.h`;
 - `StPicoBTofPidTraits.h`;
 - other PicoDst headers included in `picoAnalyzer.cpp`.

Before compilation, load the ROOT and STAR environment required by your local setup.

## Build

A basic compilation command is:

```bash
g++ -std=c++17 picoAnalyzer.cpp -o picoAnalyzerStandalone \
   $(root-config --cflags --libs) \
   -lStPicoDst
```

If the compiler cannot find PicoDst headers or libraries automatically, provide the paths explicitly:

```bash
g++ -std=c++17 picoAnalyzer.cpp -o picoAnalyzerStandalone \
   $(root-config --cflags --libs) \
   -I/path/to/StPicoDst/include \
   -L/path/to/StPicoDst/lib \
   -lStPicoDst
```

Replace `/path/to/StPicoDst/include` and `/path/to/StPicoDst/lib` with the actual paths on your system.

## Run

The executable expects two command-line arguments:

```bash
./picoAnalyzerStandalone inputFileName outputFileName.root
```

where:

- `inputFileName` is a PicoDst ROOT file or an input file list supported by `StPicoDstReader`;
- `outputFileName.root` is the output ROOT file containing histograms.

Example:

```bash
./picoAnalyzerStandalone input.list output.root
```

After the analysis finishes, open the output file in ROOT:

```bash
root output.root
```

Example ROOT session:

```cpp
TFile *f = new TFile("output.root");
hInvMassEEFinal->Draw();
```

## Analysis workflow

The analysis is performed in two passes over the input data.

### First pass

The first pass is used to prepare the event centrality binning. It:

1. reads PicoDst events;
2. fills the `refMult` distribution;
3. applies basic event cuts;
4. fills the `refMult` distribution after event selection;
5. calculates centrality-bin boundaries from the integral of the `refMult` distribution.

### Second pass

The second pass performs the main analysis. It:

1. applies event cuts;
2. loops over tracks in each event;
3. applies primary-track quality cuts;
4. applies electron PID cuts;
5. builds a list of electron and positron candidates;
6. rejects conversion-electron candidates;
7. fills the same-event `e+e-` invariant-mass spectrum;
8. fills the mixed-event background spectrum;
9. normalizes the mixed-event background;
10. writes all histograms to the output ROOT file.

## Event selection

An event is accepted if the reconstructed primary vertex satisfies:

```
|Vz| < 40 cm
sqrt(Vx^2 + Vy^2) < 2 cm
```

## Track selection

A track is accepted if it satisfies:

```
isPrimary == true
|eta| < 1.0
0.2 < pT < 50 GeV/c
DCA < 1.0 cm
nHitsFit / nHitsMax > 0.52
nHitsFit >= 20
nHitsDedx >= 11
```

## Electron PID selection

An electron candidate is accepted if it satisfies:

```
-1.9 < nSigmaElectron < 3
track has TOF PID information
|1 / beta - 1| < 0.03
|btofYLocal| < 10
```

The code also fills BEMC-related `P/E` histograms. The corresponding BEMC cut is present in the source code but is currently commented out.

## Conversion-electron rejection

To suppress electrons from photon conversions, electron candidates are checked pairwise. If a candidate forms an opposite-sign pair with invariant mass

```
M(e+e-) < 0.05 GeV/c^2
```

then this candidate is rejected from the final same-event and mixed-event spectra.

## Event mixing

The combinatorial background is estimated with event mixing. Events are grouped by:

- `Vz` bin;
- centrality bin determined from `refMult`.

For each `(Vz, centrality)` bin, the analyzer stores a pool of electron candidates from previous events.



## Main parameters

The main binning and event-pool parameters are defined at the beginning of the source file:

```cpp
const int nVzBins = 10;
const int nCentBins = 9;
const int poolDepth = 10;
```

The centrality percentiles are defined as:

```cpp
double centEdges[] = {0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80};
```

## Output histograms

The analyzer writes event, track, PID, detector, and physics histograms to the output ROOT file.

### Event QA

- `hRefMult` reference multiplicity after event cuts;
- `hVxVy_before` primary-vertex `Vx` versus `Vy` after event cuts;
- `hVz_before` primary-vertex `Vz` after event cuts;
- `hVtxXvsY` primary-vertex `Z` distribution.

### Track QA

- `hGlobalPtot` global-track momentum distribution after selection;
- `hPrimaryPtot` primary-track momentum distribution after selection;
- `hTransvMomentum` `phi` versus `pT` for different charge signs;
- `hNfit_before`, `hNfit_after` `nHitsDedx` before and after track cuts;
- `hDCA_before`, `hDCA_after` `pT` before and after track cuts;
- `hEta_before`, `hEta_after` `pT` versus `eta` before and after track cuts.

### PID QA

- `hNSigmaElectron` `nSigmaPion` distribution;
- `hNSigmaKaon` `nSigmaProton` distribution;
- `hDedxVsPq_before` TPC `dE/dx` versus `p/q` after PID cuts;
- `hNsigmaE_vs_Pq_before` `nSigmaElectron` versus `p/q` after PID cuts;
- `hInvBeta_vs_Pq_before` TOF `1 / beta` versus `p/q` after PID cuts;
- `hEoverP_vs_Pq_before` BEMC `P/E` versus `p/q` after PID cuts;
- `hTofBeta` BTOF tray hit distribution;
- `hBTowAdc` FMS ADC distribution;
- `hETofToT` EPD ADC distribution.

