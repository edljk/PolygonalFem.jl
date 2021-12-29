# PolygonalFem

[![Build Status](https://travis-ci.com/edljk/PolygonalFem.jl.svg?branch=main)](https://travis-ci.com/edljk/PolygonalFem.jl)
[![Coverage](https://codecov.io/gh/edljk/PolygonalFem.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/edljk/PolygonalFem.jl)

Adapted code from the article / matlab's code `The virtual element method in 50 lines of matlab`. See https://arxiv.org/pdf/1604.06021.pdf.

## Install PolygonalFem

```julia
>] add https://github.com/edljk/PolygonalFem.jl
```

## Benchmark

```julia
> using PolygonalFem
> PolygonalFem.vem("squarepolmesh_lesscoarse", 1000)
```

### Results with random Voronoi cells

<table>
<tr>
    <td>
        <img src="https://github.com/edljk/PolygonalFem.jl/blob/main/test/figures/squarepolmesh_coarse_100.png" style=width:200px>
     </td>
     <td>
         <img src="https://github.com/edljk/PolygonalFem.jl/blob/main/test/figures/squarepolmesh_coarse_1000.png"  style=width:200px>
    </td>
     <td>
         <img src="https://github.com/edljk/PolygonalFem.jl/blob/main/test/figures/squarepolmesh_coarse_10000.png"  style=width:200px>
    </td>
</tr>
</table>

### Results with regular Voronoi cells (mesh obtained by 1000 iterations of Lloyd's  algorithm)
<table>
<tr>
    <td>
         <img src="https://github.com/edljk/PolygonalFem.jl/blob/main/test/figures/squarepolmesh_100.png"  style=width:200px>
     </td>
     <td>
         <img src="https://github.com/edljk/PolygonalFem.jl/blob/main/test/figures/squarepolmesh_1000.png"  style=width:200px>
    </td>
     <td>
         <img src="https://github.com/edljk/PolygonalFem.jl/blob/main/test/figures/squarepolmesh_10000.png"  style=width:200px>
    </td>
</tr>
</table>


### Results with regular Voronoi cells (mesh obtained by 5 iterations of Lloyd's algorithm)
<table>
<tr>
    <td>
        <img src="https://github.com/edljk/PolygonalFem.jl/blob/main/test/figures/squarepolmesh_lesscoarse_100.png" style=width:200px>
     </td>
     <td>
         <img src="https://github.com/edljk/PolygonalFem.jl/blob/main/test/figures/squarepolmesh_lesscoarse_1000.png"  style=width:200px>
    </td>
     <td>
         <img src="https://github.com/edljk/PolygonalFem.jl/blob/main/test/figures/squarepolmesh_lesscoarse_10000.png"  style=width:200px>
    </td>
</tr>
</table>