# Project Conclusions

## The graph

![Time required per input given a set of processes](img/graph.png)

## The "Why"
As we can see, as the number of processes increases, we can see an increase in the time taken to
complete the program, as process creation, processing, commnication, and termination results in
an overhead for the program that thus causes it to slow as `p -> inf`, where `p = number_of_processes`.

Another trend to see is that as the datasets increase with a constant number of processes, we see that
for the datasets given, the time required for them to complete is around the same, with a few instances
of being faster than the smaller inputs. A possibility for why this is so is because the "huge" dataset
given to us is actually quite small when my algorithm is used, since I've priorised the minimisation of
communication between processes to the point that there is only a single `MPI` call that sends and receives,
namely `MPI_Allgather()`. This means that it will scale quite well with datasets as they grow.
