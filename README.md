# 3D → 2D Converter

A MATLAB script that loads a triangle‐mesh in RAW format, applies world, view, and projection transforms, performs backface & frustum culling, depth‐sorts triangles, and renders a Lambert‐shaded 2D projection.

---

## Table of Contents

- [Features](#features)  
- [Prerequisites](#prerequisites)  
- [Usage](#usage)  
- [Configuration](#configuration)  
- [Example Scenes](#example-scenes)  
- [Directory Structure](#directory-structure)  
- [License](#license)  
- [Author](#author)  

---

## Features

- **World Transform**: rotation, scaling, translation  
- **Camera**: look-at view matrix from configurable eye position  
- **Projection**: perspective frustum (FOV, aspect, near/far)  
- **Culling**: backface and view-frustum culling  
- **Depth Sorting**: Painter’s algorithm (sort by average depth)  
- **Shading**: Lambertian diffuse lighting  
- **Input**: any `.raw` mesh file with one triangle per line (`x1 y1 z1 x2 y2 z2 x3 y3 z3`)

---

## Prerequisites

- MATLAB R2018a or later  
- A RAW mesh file (`.raw`), e.g.:
