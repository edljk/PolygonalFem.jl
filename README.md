# PolygonalFem

[![Build Status](https://travis-ci.com/edljk/PolygonalFem.jl.svg?branch=main)](https://travis-ci.com/edljk/PolygonalFem.jl)
[![Coverage](https://codecov.io/gh/edljk/PolygonalFem.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/edljk/PolygonalFem.jl)

Adapted code from the article / matlab's code `The virtual element method in 50 lines of matlab`. See https://arxiv.org/pdf/1604.06021.pdf.

### Results with regular Voronoi cells (obtained with ~ 1000 Lloyd's iterations)
<table>
<tr>
    <td>
        <img src="https://user-images.githubusercontent.com/14992507/147603571-bb214ddc-2c5f-406c-aff6-e7e4810f6269.png" style=width:200px>
     </td>
     <td>
         <img src="https://user-images.githubusercontent.com/14992507/147603595-a5d61609-8b81-4a62-b9ff-86fccec29967.png"  style=width:200px>
    </td>
</tr>
</table>

### Results with coarse Voronoi cells (obtained with ~ 100 Lloyd's iterations)

<table>
<tr>
    <td>
        <img src="https://github.com/edljk/PolygonalFem.jl/blob/main/test/figures/squarepolmesh_coarse_100.png" style=width:200px>
     </td>
     <td>
         <img src="https://github.com/edljk/PolygonalFem.jl/blob/main/test/figures/squarepolmesh_coarse_1000.png"  style=width:200px>
    </td>
</tr>
</table>