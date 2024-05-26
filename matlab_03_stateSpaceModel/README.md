To run this version, run the `main.m` file. This program takes a different approach that the genetic algoritm. The genetic algoritm searched for constant parameters of the model transfer function, while this program searches for the state-space model matrices.
The `initial_guess` variable is already set to a good starting point. You may change it if you want to.
In the `learn.m` file, you may give the funtion a constant mass, spring constant, or other parameters. The program will then focus on the other parameters better.
The `G`, `G_2`, `p` and `P` parameters control the learning depth and speed. Changing them effects the speed and accuracy of the program.
This program was also optimized to rung well on GPU. Use the `main_gpu.m` file to run the program on an nvidia gpu with CUDA support.