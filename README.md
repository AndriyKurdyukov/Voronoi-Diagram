<img width="320" height="320" alt="voronoi" src="https://github.com/user-attachments/assets/0e533b03-2e9d-4f8e-8991-666b01852149" />

# A Fortune's algorithm implementation

## Description

An implementation of Fortune's algorithm for Voronoi diagrams in a 2D plane. The implementation is simplified compared to original algorithm, since a linear list is used
instead of a binary tree of arcs. Because the search in a linear list of arcs takes a runtime of O(n), this results in a worse runtime of O(n^2) instead of possible O(n log(n)) if using a binary tree of arcs instead. \
Also some special edge cases are not yet handled, e.g the case of two generator points having maximum and equal  y-coordinate, although I might add a handling later(Edit: starting with v0.1.0, handling for the case of equal y-coordinates is implemented).

## Dependencies

- SFML 2.6 graphics library
- GNU Make

## Usage

- "./all [-a] filename" - build voronoi diagram, with flag "-a" for optional computation of voronoi cell areas(in pixels^2).
The infinite areas are shown as "-1".
"filename" is expected to be a file with point coordinates in the following format: "x,y" per every new line.
- "./all -v" - output program version to the console.

