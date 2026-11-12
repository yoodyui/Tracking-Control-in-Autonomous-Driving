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
- `Ex7_1_LQR_kinematic_weights.m` (the three-weight comparison in Example 7.1)
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

## Alignment with the revised manuscript (12 November 2026)

- Example 3.3 now uses the manuscript's nominal speed, initial conditions, inputs,
  and independently propagated continuous-time/Euler state histories.
- Example 5.1 runs both the unlimited-acceleration and limited-acceleration cases.
- The unfiltered time-indexed ILC scripts reset and truncate lap logs and learn
  only from samples actually executed. Their scalar, unshifted learning law is
  a baseline demonstration, not a guaranteed-convergent design.
- In `Ex7_1_LQR_kinematic_weights.m`, set `remove_actuator_limits = true` for
  Exercise 7.3. The reference and nominal input trajectory remain fixed so the
  effect of applied-command saturation can be isolated. Unconstrained plots are
  written in a separate subfolder and do not overwrite the worked-example plots.

The printed listings sometimes omit setup/plotting and use shorter variable
names. Repository scripts are intended as complete runnable versions.

### MATLAB dependencies

Control System Toolbox is needed for `ctrb`/`dlqr`; Optimization Toolbox for
`quadprog`; Signal Processing Toolbox for `butter`/`filtfilt` (and `rms` where
required by the installed MATLAB version). Some scripts use `wrapToPi`, commonly
provided by Mapping Toolbox. MATLAB version and toolbox availability must be
checked on the machine used for numerical reproduction.

Validation of this revision included source comparison and independent numerical
checks of model Jacobians, PID saturation/anti-windup, and the scalar CBF example.
A full MATLAB run was not completed in the editing environment because MATLAB
terminated at startup. Reproduction of all textbook plots is still required.
