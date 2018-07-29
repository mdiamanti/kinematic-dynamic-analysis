# kinematic-dynamic-analysis
Current repository includes MATLAB files related to the kinematic and dynamic analysis of Kinova Jaco<sup>2</sup> 6 DOF robotic arm. Α PD controller scheme is, also, included so as to test the behavior of the computed dynamic model.

The appropriate kinematic parameters and other technical specifications of Jaco<sup>2</sup> 6 DOF can be found within the datasheet `JACO²-6DOF-Advanced-Specification-Guide.pdf`. Other useful information relative to j2n6s300* model was extracted from [Github kinova-ros package](https://github.com/Kinovarobotics/kinova-ros).

## Details about MATLAB files
-`kinematic_analysis_j2n6s300.m` file contains the kinematic analysis, while `dynamic_analysis_j2n6s300.m` file contains the dynamic analysis of Jaco<sup>2</sup> 6 DOF robotic arm. 

-All the other files constitute functions, either autogenerated from MATLAB or not. 

-The autogenerated ones are used in the context of the Simulink PD controller scheme, `pd_controller_j2n6s300.slx`. In this way, there is no need to execute any other file before using this Simulink model.

</br>
</br>
* j2n6s300 consists an abbreviation that refers to the model Jaco<sup>2</sup> 6 DOF with 3 robot fingers.
