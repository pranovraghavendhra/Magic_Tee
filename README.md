# Magic Tee (Hybrid Tee) Formation and Analysis

A MATLAB model of a magic tee, a four-port waveguide junction that combines an E-plane tee and an H-plane tee. The project builds the ideal and non-ideal S-matrix, checks its properties, and plots the results along with a 3D model of the junction.

## Ports

| Port | Name | Role |
|---|---|---|
| 1 | Collinear arm | Input / output |
| 2 | Collinear arm | Input / output |
| 3 | H-arm | Sum port |
| 4 | E-arm | Difference port |

## What the code does

- Builds the ideal 4x4 S-matrix of a magic tee.
- Builds a non-ideal version with insertion loss, finite isolation and port mismatch.
- Checks symmetry (reciprocity), unitarity (losslessness), passivity and isolation.
- Computes output power for five excitation cases.
- Sweeps the input phase difference to show sum and difference operation.
- Computes WR-90 (TE10) cutoff frequency, guide wavelength, phase and group velocity.
- Draws a rotatable 3D model of the tee with labelled ports.

## Requirements

- MATLAB R2018b or newer (`xline` needs R2018b).
- No extra toolboxes.

## How to run

1. Clone or download this repository.
2. Open MATLAB and set the repository folder as the current folder.
3. Run:

```matlab
MagicTee_Analysis
```

The Command Window prints the S-matrices, property checks, power table and waveguide data, and six figures open.

## Outputs

| Figure | Content |
|---|---|
| 1 | 3D geometry of the magic tee (WR-90 cross-section) |
| 2 | S-matrix magnitudes, ideal vs non-ideal |
| 3 | Output power distribution for five excitations |
| 4 | Sum and difference output vs input phase |
| 5 | WR-90 guide wavelength, phase and group velocity |
| 6 | Key S-parameters in dB |

## Ideal S-matrix

```
S = (1/sqrt(2)) * [ 0   0   1   1
                    0   0   1  -1
                    1   1   0   0
                    1  -1   0   0 ]
```

- Input at Port 3: equal split to Ports 1 and 2, in phase. Port 4 isolated.
- Input at Port 4: equal split to Ports 1 and 2, 180 degrees out of phase. Port 3 isolated.
- Ports 1 and 2 together: sum at Port 3, difference at Port 4.

## Limitations

- The non-ideal values (0.15 dB insertion loss, 28 dB isolation, 22 dB return loss) are assumed, not measured.
- The non-ideal model rescales the matrix to stay passive, so the delivered power is lower than the nominal insertion loss implies.
- The S-matrix is a circuit-level model. It is not derived from the 3D geometry, and the geometry has no internal matching elements.
- The model is not frequency-dependent.
- Results have not been validated against a field solver such as HFSS or CST.

## Applications

Power divider and combiner, balanced mixer, duplexer, impedance bridge, monopulse radar comparator.

## Repository structure

```
MagicTee_Analysis.m           Main MATLAB file
MagicTee_Presentation.pptx    Project presentation
README.md                     This file
```
