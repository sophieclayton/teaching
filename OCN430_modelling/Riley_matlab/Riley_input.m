% Simple function to receive input for the model sensitivity analysis

function test = Riley_input
test = input(sprintf('Which parameter would you like to test?\n 1. P0, inital phytoplankton concentration \n 2. p, photosynthetic constant \n 3. R0, respiratory rate \n 4. g, grazing rate \n >> '));

end