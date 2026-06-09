<img width="320" height="320" alt="voronoi" src="https://github.com/user-attachments/assets/0e533b03-2e9d-4f8e-8991-666b01852149" />

# A Fortune's algorithm implementation

## Description

An implementation of Fortune's algorithm for Voronoi diagrams in a 2D plane. The implementation is simplified compared to original algorithm, since a linear list is used
instead of a binary tree of arcs. This results in a worse runtime of O(n^2 log(n)) instead of possible O(n log(n)). \
Also some special edge cases are not yet handled, e.g the case of two generator points having maximum and equal  y-coordinate, although I might add a handling later.

## Usage

TBD
