# Model Predictive Control (MPC)

## About MPC

Model Predictive Control (MPC) uses the dynamic process of the dynamical model and the state during a finite time-horizon. This way, by knowing the state for the next n steps, we can optimize the actuators and anticipate changes in the trajectory by taking into consideration the dynamical model of the car.

## Compiling and executing the project

- Clone the repo and cd to it on a Terminal.
- Create the build directory: `mkdir build`.
- `cd build`
- `cmake ..`
- `make`: this will create the executable `mpc`, that contains a program that reads the input from the simulators and computes the mpc.

`mpc` communicates with a Udacity simulator that reads the input data and produces estimations for the object and trajectories. See [Udacity's seed project](https://github.com/udacity/CarND-MPC-Project) for more information on how to install the simulator.

Open the simulator and run `./mpc` to enable the connection. Select Project 5 from the simulator.

## Implementation

The project computes a few values; some of them are actuators to modify the trajectory of the car, and the others are trajectory information:

- Steering angle: what is the angle we need to move the car to follow the desired trajectory.
- Throttle: how much we need to accelerate or break the car to follow the desired trajectory.
- Trajectory to follow: to be represented in the simulator to get an idea of the movement the car should follow.
- Trajectory the car would follow with the current actuator (mpc predicted trajectory)

The code tries to find the actuators that, given certain dynamics of the model, would minimize the trajectory error, by using an approximation of the trajectory using polynomial curve fitting. This part is implemented in MPC::Solve.

I tried to organized the code to make it easier to follow, by dividing a bit the code into logical methods, to make the whole structure more readable.


## Problems during the implementation and Results

During the implementation, there were a few things that made my work a bit harder; the current implementation does not work when I increment the number of points in the future to 20, and if I decrease the dt. I didn't play a lot with that but maybe using timestamps of the input values would allow more precise dt? 

Another bug I had for a while, was the use of numerical limits for lower and upper bounds for the vars. I wasted a lot of time assuming that numeric_limits<double>::min() and max() were what I needed, but the code did not work (weird actuator signals).

Once these two bugs were addressed, the car could drive around the track flawlessly. I tried increasing the speed up to 90mph with no problems. I ended up implementing the MPC with dt = 0.1 and N = 10, which gave good results.

![Simulator result turning left](images/turn_left.png)
![Simulator result top speed](images/top_speed.png)

