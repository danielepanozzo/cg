![Build](https://github.com/danielepanozzo/cg/workflows/Build/badge.svg)

# Graduate Computer Graphics, CSCI-GA 2270-001 Spring 2023

### Course Instructor
*Daniele Panozzo*

60 5th Ave, 5th Floor

[panozzo@nyu.edu](mailto:panozzo@nyu.edu)

URL: [http://cs.nyu.edu/~panozzo/](http://cs.nyu.edu/~panozzo/)

### Assistants
*Arvi Gjoka*
[arvi.gjoka@nyu.edu](mailto:arvi.gjoka@nyu.edu)

*Jiacheng Dai*
[jd4705@nyu.edu ](mailto:jd4705@nyu.edu )

Office Hours (starting from February 2nd): 
* **Mondays**  Arvi 9-10AM, Discord
* **Tuesdays** Arvi 3-4PM, 60 5th Ave, 502
* **Wednesdays** Jiacheng 2-3PM, 60 5th Ave, 502 
* **Thursdays** Daniele 3-4PM, 60 5th Ave, 504  
* **Fridays** Jiacheng 2-3PM, 60 5th Ave, 502

### Lectures:
Thursdays, 4:55 PM - 6:55 PM, GCASL, Room 269

## Course Description

This course provides an introduction to the field of Computer Graphics. We will cover the basic mathematical concepts, such as 2D and 3D transformations, study the interaction of light with geometry to derive  shading models, and implement rendering algorithms such as ray tracing and rasterization. We will investigate how these fundamental components are integrated in current graphics processors and study the corresponding programming APIs. This course will also include a brief introduction to C++.

By the end of the course, the student must be able to:

* Explain and apply the fundamental mathematical concepts used in  image synthesis algorithms
* Implement a basic rendering system based on ray tracing
* Implement a basic rendering pipeline based on rasterization
* Develop graphics programs in C++ using the Simple Directmedia Library (SDL)

*Textbook*:
Fundamentals of Computer Graphics, 4th Edition
December 18, 2015 by A K Peters/CRC Press
Textbook - 734 Pages - 541 Color
ISBN 9781482229394

## Assignments

In this course there will be 4 mandatory coding assignments that will account for 70% of the grade. There will also be 2 optional assignments, for a total of 15% extra points.

Please use the script `prepare_assignment.sh` to zip the submission folder, specifying the folder to zip as the first argument and the writeup document as the second document. For example, `sh prepare_assignment.sh Assignment_1 Assignment_1/writeup.pdf`. This will delete an existing zip file with the same name (in the case above, any existing `Assignment_1.zip` will be deleted).
This may not work on Windows, so in that case make sure that you zip the parent folder and either indicate or provide the writeup separately in the submission. Please archive in to a `.zip` file in that case.

## Course Notes:

The course notes will be adjusted along the way. The order is indicative and might change.

[01 - Introduction](https://www.icloud.com/keynote/0Bi3HXvG70bpshIbt1t9PnGmw#01_-_Introduction_to_Computer_Graphics)

[02 - Basic Linear Algebra](https://www.icloud.com/keynote/0bR6rH_qhMGyack3AvOLN9KpA#02_-_Basic_Linear_Algebra)

[03 - Introduction to C++](https://www.icloud.com/keynote/0g2wBvEMQe7c4KRNidmCT44rQ#03_-_C++)

[04 - Images](https://www.icloud.com/keynote/078471RTY56oFkHbVjhquf4Lg#04_-_Images)

[05 - Ray Tracing](https://www.icloud.com/keynote/0Xt7leP_xqOA9pEE24U9-q5vg#05_-_Ray_Tracing)

[06 - Spatial Data Structures](https://www.icloud.com/keynote/0WGDZa8VZoXxqlLSq2gp_G_Rw#06_-_Spatial_Data_Structures)

[07 - Procedural Synthesis](https://www.icloud.com/keynote/0RV7ZnHhuQCWAHlj29VpmmLKQ#07_-_Procedural_Synthesis)

[08 - 2D Transformations](https://www.icloud.com/keynote/0hdbFFSx6TrJSzmICBf4Yjo2g#08_-_2D_Transformations)

[09 - Viewing Transformations](https://www.icloud.com/keynote/0DlviF0tU_vb8pn-w6qoN3OaA#09_-_Viewing_Transformations)

[10 - Rasterization - Theory](https://www.icloud.com/keynote/0gwK2pQbGYorL7xXpoAlb8xog#10_-_Rasterization_-_Theory)

[11 - Rasterization - Implementation](https://www.icloud.com/keynote/0WuGGx7-YzpkpxN5lyMvQyHew#11_-_Rasterization_-_Implementation)

[12 - Perspective Projection](https://www.icloud.com/keynote/0qLSBn6y3y4Fn-tvY_ZOtKlFQ#12_-_Perspective_Projection)

[13 - Texture Mapping](https://www.icloud.com/keynote/0DYWwAzDUEQC5AMmDMIrFyI0w#13_-_Texture_Mapping)

[14 - Parametrization](https://www.icloud.com/keynote/06w_Cnj7E81JLB76Pl1tSjiqg#14_-_Parametrization)

[15 - Designing Interpolating Curves](https://www.icloud.com/keynote/0ztivLSI82_YZgJLkU3MsemQQ#15_-_Designing_Interpolating_Curves)

[16 - Designing Approximating Curves](https://www.icloud.com/keynote/0KOu9icbYhEryoNeTVsjftmRQ#16_-_Designing_Approximating_Curves)

[17 - Designing Surfaces](https://www.icloud.com/keynote/0fjx3PAnYzwWgrgiIXEJI7K6g#17_-_Designing_Surfaces)

[18 - Mesh Deformation](https://www.icloud.com/keynote/0i9XXnime7phhNqYRM1m-a7FA#18_-_Mesh_Deformation)

[19 - Mesh Data Structures](https://www.icloud.com/keynote/07fbGW6rsRjxWrmTu3M8kDgvA#19_-_Mesh_Data_Structures)


## Assignments

[General Instructions](https://github.com/danielepanozzo/cg/tree/master/RULES.md)

[Assignment 1: Introduction to C++, Geometry (Optional)](https://github.com/danielepanozzo/cg/tree/master/Assignment_1)

[Assignment 2: Introduction to Raytracing and Shading](https://github.com/danielepanozzo/cg/tree/master/Assignment_2)

[Assignment 3: Basic Raytracing Effects (Optional)](https://github.com/danielepanozzo/cg/tree/master/Assignment_3)

[Assignment 4: Ray Tracing Triangle Meshes and AABB Trees](https://github.com/danielepanozzo/cg/tree/master/Assignment_4)

[Assignment 5: Rasterization](https://github.com/danielepanozzo/cg/tree/master/Assignment_5)

[Assignment 6: 2D Shape Editor](https://github.com/danielepanozzo/cg/tree/master/Assignment_6)

Assignment 6 can be replaced by a final project.
