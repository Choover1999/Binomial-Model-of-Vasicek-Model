This is code which can create the binomial Vasicek model. I want to disclose that chat gpt was used to create the upper and lower solver functions for the branches because it seemed much faster to do this that to try to solve the system of equations at each external node analytically. The code for the tree itself was created by me, but does use the functions created by chatGPT. 

Also, be aware that if you need to store the probabilities of each branch, more code will need to be added to make sure this is accomplished correctly for the purposes desired. Path probabilities in this model can become quite complicated, so it's important to understand how the model works and what it will be needed for, before jumping in.
