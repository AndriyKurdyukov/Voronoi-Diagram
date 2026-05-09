
<img width="320" height="320" alt="voronoi" src="https://github.com/user-attachments/assets/0e533b03-2e9d-4f8e-8991-666b01852149" />



An implementation of Fortune's algorithm for Voronoi diagrams in a 2D plane. The implementation is simplified comapred to original algorithm, since a linear list is used
instead of a binary tree of arcs. This results in a worse runtime. 
Also some special edge cases are not handled, e.g the case of two generator points having maximum and equal  y-coordinate,although I might add
a handling later.
