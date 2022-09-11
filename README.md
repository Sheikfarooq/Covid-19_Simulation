
# Covid-19 spread simulation

Simulating the people (points), considering they are  moving in a space (box).
The people are walking at a random speed (𝑣_𝑥, 𝑣_𝑦). When the persons bump into the wall of the box, they bounce off it back again.
At the beginning of the simulation, only a few individuals are sick (represented in red colored dots).
Once these sick persons (red dots) come into the proximity of healthy persons (green dots), these also become infected (becomes red).
After a certain time, the sick persons (red dots) recover (becomes blue) and can no longer infect or become infected.


## Initialization

 Each individual (number of individuals given in each task) has been initialized with the following properties:
 
Random 𝑥 and 𝑦 position in a given area

Random 𝑣_x and 𝑣_𝑦 speeds in the interval [−1,1]

10% of the individuals are initially infected, the others are healthy
## Tools

MATLAB 2018b
## Color Reference

| Color             | Description                                                                |
| ----------------- | ------------------------------------------------------------------ |
| Red | Infected |
| Green | Healthy |
| Blue | Recovered |

## Result
