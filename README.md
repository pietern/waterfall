# waterfall

Interactive frequency spectrum analysis.

## Usage

Dependencies:

* GTK+ 3
* FFTW 3

To build, run `make` in the project root.

To execute, run `./waterfall` and specify a signal file to analyze.

### Supported formats

* Wave (WAV)
* Raw 8-bit signed integer
* Raw 16-bit signed integer
* Raw 32-bit signed integer
* Raw 32-bit float

Wave files are expected to have 2 channels (in-phase and quadrature).

Raw files are expected to have interleaved in-phase and quadrature samples.

## License

Simplified BSD. See `LICENSE` file.
