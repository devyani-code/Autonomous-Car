# Autonomous-Car
Worked upon the problem of robot sensation and movement and probabilistic techniques to account for inaccuracies in both.

Localization.py: Contains an python code for Grid based approach of localization.

State_estimation_Kalman_filter.py: Implementation of powerful method of localization called the Kalman filter. It uses sequential matrix operations to predict the position and velocity of a robot whose motion is linear but has unknown position and velocity.

Particle_filter_localization.py:. In this technique, different possible positions of the robot (the "particles") are randomly generated and then given a score based on how compatible these positions are with the robot's observed measurements. The particles then are reselected with a probability corresponding to their score, and this process (after many iterations) causes particles to cluster around a position that is highly likely to be the robot's actual position
