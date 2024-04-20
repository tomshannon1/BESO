# BESO
## Topology Optimization for Python using Bidirectional Evolutionary Structure Optimization (BESO)

This is a Python package that allows you to run 3D structural topology optimization problems using the bi-directional evolutionary structure optimization (BESO) method. This code was written based on work done for the 2D case by Xie and Huang in the book, Evolutionary Topology Optimization of Continuum Structures: Methods and Applications. You can learn more about this method and topology optimization in my undergraduate thesis: A comparison of topology optimization methods for the design of a cantilever beam.

Below is an example of a cantilever beam generated with this optimization process. 

![image](https://github.com/tomshannon1/BESO/assets/18470042/515f0146-b37a-42f6-8fa9-dc86835b7850)

This is a final three-dimensional design of a cantilever designed from this optimization method. The software exports VTK files for every design iteration for easier computations and faster export times. However, the user can then modify these files with software, such as Paraview, to make the design look like the model above and suitable for 3D printing. The above image has a poisson surface recontruction filter to smooth jagged edges of the VTK file. A simple marching cubes algorithm could also be used for a similar effect.

![](https://lh6.googleusercontent.com/w3bY1uCfacg6dtadv0kjLqBv6_srCsDfL5-wGSmNVUGUlAAkzM3ktf9j7yQ_e43cHBnUfMLz3u4Hw357oZ4bJGKPXOeHWQXK7Y54rwI5Ipp8QuDFziJoqi8WCO8vMp45qnS7SBksYwQ)

This is an animation of the design process. This design process starts with a 3D grid of voxel elements, composed of eight nodal points. The user specifies the nodal points on which a given load on the structure is placed. Based on this load, the optimization process subtracts inefficient material and adds material in place where the material is needed. You can give a particular volume fraction, which is the percentage of material to keep from the original design. This design above is a cantilever beam set to be designed at 15% volume fraction.

## Future work

This software is quite resource intensive, depending on the number of elements used in the design. Therefore, this software would highly benefit from running on a GPU instance. Additionally, this could benefit by using sparse matrices. In the past, I have tried to refactor to make it work with sparse matrices, but the solvers did not produce the correct output and made physically incorrect designs.  

## Contributers
<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tr>
    <td align="center"><a href="https://github.com/tomshannon1"><img src="https://avatars3.githubusercontent.com/u/18470042?s=460&v=4" width="100px;" alt=""/><br /><sub><b>tomshannon1</b></sub></a><br />
  </td>
  </tr>
</table>
<!-- markdownlint-enable -->
<!-- prettier-ignore-end -->
<!-- ALL-CONTRIBUTORS-LIST:END -->

If you want to contribute, message me (@tomshannon1) or send a pull request!
