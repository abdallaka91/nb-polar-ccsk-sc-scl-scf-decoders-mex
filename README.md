# NB Polar CCSK Frequency-Domain Simulator

MATLAB Monte Carlo simulator for non-binary polar codes with CCSK modulation and frequency-domain decoding through a Windows MEX/DLL decoder.

## Main Files

- `NB_polar_CCSK.m`: main simulator script.
- `support_files/nbdecode_mex*.mexw64`: MATLAB MEX decoder entry point.
- `support_files/libNbScFFTdec.dll`: decoder library loaded by the MEX file.
- `support_files/matrices/`: reliability matrices indexed by GF size and code length.
- `support_files/ccsk_sequences.txt`: CCSK reference sequences.
- `support_files/`: MATLAB helper functions.

## Basic Use

Edit the configuration section at the top of `NB_polar_CCSK.m`:

```matlab
code_length = 128;         % N
transmitted_length = 128;  % Ns
information_length = 90;   % K
gf_size = 64;              % q
snr_db = -7.5;
decoder_type = 'dec4';
max_frames = 4e4;
frames_per_call = 2000;
num_threads = 4;
```

Then run:

```matlab
NB_polar_CCSK
```

The MEX expects decoder input as `q x Ns x frames_per_call` single-precision probabilities. The simulator forms this internally.

## Decoder Files Needed On Windows

The MEX file and DLL are stored together in `support_files/`, which the simulator adds to the MATLAB path:

```text
nbdecode_mex1.mexw64
libNbScFFTdec.dll
```

Depending on the target machine, MinGW runtime DLLs may also be needed beside the MEX file in `support_files/`:

```text
libgcc_s_seh-1.dll
libstdc++-6.dll
libwinpthread-1.dll
```

## Supported Decoders

See the comment block at the top of `NB_polar_CCSK.m` for the supported SC, SCL, and SCF decoder strings and their parameters.

## Authors

Authors and contributors include Abdallah Abdallah, Bertrand Le Gal, Camille Moniere, and Emmanuel Boutillon. The non-binary polar decoder software used by this simulator is mainly from Bertrand Le Gal and collaborators.

## Licence

This project is distributed under the CeCILL-B free software licence. See `support_files/LICENSE` and the official licence text:

https://cecill.info/licences/Licence_CeCILL-B_V1-en.html

## Acknowledgement

This work has been funded by the French ANR project MIOT under grant Projet-ANR-24-CE93-0017.

Project website: https://project.inria.fr/miot/