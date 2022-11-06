#For the purpose of this assume that the robot can move only left, right, up, or down. It cannot move diagonally. Also, for this assignment, the robot will never overshoot its destination square; it will either make the movement or it will remain stationary.
# The function localize takes the following arguments:
#
# colors:
#        2D list, each entry either 'R' (for red cell) or 'G' (for green cell)
#
# measurements:
#        list of measurements taken by the robot, each entry either 'R' or 'G'
#
# motions:
#        list of actions taken by the robot, each entry of the form [dy,dx],
#        where dx refers to the change in the x-direction (positive meaning
#        movement to the right) and dy refers to the change in the y-direction
#        (positive meaning movement downward)
#        NOTE: the *first* coordinate is change in y; the *second* coordinate is
#              change in x
#
# sensor_right:
#        float between 0 and 1, giving the probability that any given
#        measurement is correct; the probability that the measurement is
#        incorrect is 1-sensor_right
#
# p_move:
#        float between 0 and 1, giving the probability that any given movement
#        command takes place; the probability that the movement command fails
#        (and the robot remains still) is 1-p_move; the robot will NOT overshoot
#        its destination in this exercise
#
# The function should RETURN (not just show or print) a 2D list (of the same
# dimensions as colors) that gives the probabilities that the robot occupies
# each cell in the world.
#
# Compute the probabilities by assuming the robot initially has a uniform
# probability of being in any cell.
#
# Also assume that at each step, the robot:
# 1) first makes a movement,
# 2) then takes a measurement.
#
# Motion:
#  [0,0] - stay
#  [0,1] - right
#  [0,-1] - left
#  [1,0] - down
#  [-1,0] - up
def move_2D(motion,p,p_move):
    row=len(p)
    col=len(p[0])
    U_y=motion[0]
    U_x=motion[1]
    aux=[[0.0 for i in range(col)]for j in range(row)] 
    p_stay=1-p_move
    for i in range(row):
        for j in range(col):
            aux[i][j]=p_move*p[(i-U_y)%len(p)][(j-U_x)%len(p[0])]+(p_stay)*p[i][j]
    
    return aux
    
def sense_2D(colors,p,measurements,sensor_right):
    row=len(p)
    col=len(p[0])
    s=0
    aux=[[0.0 for i in range(col)]for j in range(row)] 
    for i in range(row):
        for j in range(col):
            hit=(measurements==colors[i][j])
            aux[i][j]=p[i][j]*(hit*sensor_right+(1-hit)*(1-sensor_right))
            s+=aux[i][j]
    for i in range(row):
        for j in range(col):
            aux[i][j]=aux[i][j]/s
    
    return aux

def localize(colors,measurements,motions,sensor_right,p_move):
    # initializes p to a uniform distribution over a grid of the same dimensions as colors
    pinit = 1.0 / float(len(colors)) / float(len(colors[0]))
    p = [[pinit for row in range(len(colors[0]))] for col in range(len(colors))]
    
    # >>> Insert your code here <<<
    for k in range(len(measurements)):
        p=move_2D(motions[k],p,p_move)
        p=sense_2D(colors,p,measurements[k],sensor_right)
    return p
    
    
def show(p):
    rows = ['[' + ','.join(map(lambda x: '{0:.5f}'.format(x),r)) + ']' for r in p]
    print('[' + ',\n '.join(rows) + ']')
    

colors = [['R','G','G','R','R'],
          ['R','R','G','R','R'],
          ['R','R','G','G','R'],
          ['R','R','R','R','R']]
measurements = ['G','G','G','G','G']
motions = [[0,0],[0,1],[1,0],[1,0],[0,1]]
p = localize(colors,measurements,motions,sensor_right = 0.7, p_move = 0.8)
show(p) 