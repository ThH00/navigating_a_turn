# navigating_a_turn

This repository contains three codes:
1. `angled_chaplygin_sleigh`: simulating a modified Chaplygin sleigh whose knife-edge makes a constant nonzero angle with the axis of the sleigh
2. `three_wheeled_vehicle`: simulating a three wheeled vehicle modeled as a rigid body with an angled knife edge constraint at the front wheel
3. `constrained_rolling_disk`: simulating a rolling disk that is constrained to make a constant angle with the vertical and to trace a circle

## angled_chaplygin_sleigh
`get_eom.m` derives the Gibbs-Appell equations of motion of the angled Chaplygin sleigh symbolically

`solve_eom.m` solves the obtained equations of motion. the dimensions of the sleigh and other constants can be changed in this file. running this file also produces an animation of the motion.

## three_wheeled_vehicle
`get_eom.m` derives the Gibbs-Appell equations of motion of the three wheel vehicle rigid body numerically

`solve_eom.m` solves the obtained equations of motion. the dimensions of the three-wheeled vehicle can be changes in this file.

`get_animation.m` obtaines an animation of the results

## constrained_rolling_disk
`generalized_alpha_animation.py` solves for the constrained motion of the rolling disk numerically. the dimensions of the disk and other constants can be edited in this script. this simulation results are saved in the `outputs` folder.

`get_animation.m` animates the results. the disk dimensions should also be imputted in this file.