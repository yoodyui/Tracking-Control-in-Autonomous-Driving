# Companion MATLAB Code

This folder contains the MATLAB scripts used with the textbook
*Tracking Control in Autonomous Driving*. The examples are organized by
chapter so that readers can find the code that matches each worked example.

## How to Run

Open any example script in MATLAB and run it directly. Each script calls
`setup_paths.m` automatically, which adds the shared `helpers` folder to the
MATLAB path.

If MATLAB cannot find a helper function, run this once from the `codes` folder:

```matlab
setup_paths
```

## Folder Layout

- `chapter03_kinematic_modeling/`: kinematic bicycle modeling and linearized perturbation examples.
- `chapter04_vehicle_dynamics/`: dynamic bicycle model, nonlinear/linearized comparison, and discrete-time perturbation examples.
- `chapter05_longitudinal_pid/`: longitudinal PID speed-tracking example.
- `chapter06_non_model_based/`: Pure Pursuit, Stanley, PID, ILC, and spatial ILC tracking examples.
- `chapter07_model_based/`: TV-LQR and MPC examples.
- `chapter08_cbf/`: CBF and ECBF safety-filter examples.
- `helpers/`: shared model, waypoint, distance, CTE, and geometry utilities.

## Main Example Scripts

### Chapter 3

- `Ex3_1_LinearizedModelSim.m`
- `Ex3_2_LinearizedModelSim2Inputs.m`
- `Ex3_3_kinematic_bicycle_discrete.m`

### Chapter 4

- `Ex4_1_Dynamic_vs_Kinematic.m`
- `Ex4_2_FullDynamics.m`
- `Ex4_3_dynamics_bicycle_discrete.m`

### Chapter 5

- `Ex5_1_PIDLongitudinal.m`

### Chapter 6

- `Ex6_1_Pure_Pursuit.m`
- `Ex6_1_Pure_Pursuit_Long.m`
- `Ex6_2_Stanley.m`
- `Ex6_2_Stanley_Long.m`
- `Ex6_3_PID.m`
- `Ex6_3_PID_Long.m`
- `Ex6_4_PID_dynamics.m`
- `Ex6_5_PID_ILC.m`
- `Ex6_5_PID_ILC_Long.m`
- `Ex6_5_PID_ILC_Long_ZPF.m`
- `Ex6_6_PID_ILC_Long_Spatial.m`

### Chapter 7

- `Ex7_1_LQR_kinematic.m`
- `Ex7_2_LQR_dynamics.m`
- `Ex7_3_MPC_dynamicsStateErr.m`
- `Ex7_4_MPC_dynamicsCTE.m`

### Chapter 8

- `Ex8_1_PID_CBFKinematic.m`
- `Ex8_2_PID_ECBFDynamics.m`

## Notes

The public companion-code repository is intentionally limited to the chapter
folders, `helpers`, and `setup_paths.m`. Older exploratory scripts, temporary
backups, and unused variants are excluded so that the repository follows the
textbook example sequence clearly.

When adding a new example, place the script in the chapter folder that matches
the book and put reusable model or geometry routines in `helpers`.
