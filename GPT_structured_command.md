## Model Instructions

You are a research assistant writing R code for an agent-based model that simulates school
and household transmission. Here is what the code should do. Model description: We
assume that there are exactly 500 students. You should organize them into exactly 222
households. Each household contains exactly 2 adults, and therefore 944 unique
individuals. The model is seeded assuming 5 exposures, and you should start your
simulation at Day 1. On each school day (M-F), students go to class and mix with all
members of their class, with an average transmission probability of 0.01 per day. On all
days, household members mix. Once infected, an individual cannot be re-infected.

You will write R code for this agent-based model in the following steps:

1. initiate a data frame where each row corresponds to one individual. Generate a type
   column indicating student or household member. Initialize everyone to be
   susceptible. The students are assigned randomly into the 222 households based on
   siblings, following a Poisson distribution with a mean of 2 students per household. This
   means some households will have more than 2 students and some fewer, but make
   sure every household has at least one student. Create a household ID, called "HH_ID",
   which denotes the household belonging of students and other household members.
   Also, assign the 500 students randomly into 25 classes, each with exactly 20 students.

2. Write a subfunction that simulates and tracks the transmission within classroom. The
   function takes in the ID of a current infection and output the IDs of secondary
   infections. Note that students only mix in classrooms on Monday through Friday. The
   infected student only mixes with the other individuals within the same classroom. The
   transmission rate per class contact follows a Uniform(0, 0.02) per day. Then, use the
   drawn transmission rate to determine if a classroom contact has become infected by
   drawing from a binomial, where and the mean equals to the drawn transmission rate
   multiplied by the number of susceptible contacts. The result should be a vector of
   length equal to the number of contacts and each element has value either 1 denoting
   an infection or 0 otherwise. The duration of incubation periods follows Poisson(2 days)
   and the duration of infectiousness follows Poisson(5 days). Make sure you correctly
   simulate the number of days following a Poisson distribution which should have equal
   mean and variance. Keep track of the number of secondary infections in classrooms
   over weekdays vs. weekends. Keep track of the actual transmission rate in classrooms.

3. Similarly, write a subfunction that simulates and tracks the transmission within
   household. The function takes in the ID of a current infection and output the IDs of
   secondary infections. Note that all students and household members mix in household
   on all 7 days a week. The infected individual only mixes with the other individuals within
   the same household. The transmission rate in the household should be 2 times the
   drawn transmission rate in the classroom. Then, similarly, use the household
   transmission rate to determine if a household contact has become infected by drawing
   from a binomial, where the mean equals to 2 times the drawn classroom transmission
   rate multiplied by the number of susceptible contacts. The result should be a vector of
   length equal to the number of contacts and each element has value either 1 denoting
   an infection or 0 otherwise. The duration of incubation periods still follows Pois(2 days)
   and the duration of infectiousness follows Pois(5 days). Keep track of the number of
   secondary infections in households over weekdays vs. weekends. Keep track of the
   actual transmission rate in households.

4. Write a simulation function that randomly initialize 5 seed infections from all
   individuals. Call the function "run_model", which takes 4 parameters: data frame of
   agents (called "df"), number of time steps ("t"), multiplier on household transmission
   compared to classroom ("rel_HH"), and initial seed infections ("seeds"). Simulate for
   1000 runs and in each simulation, the model runs for 30 days. Make sure to output the
   data frame of all individuals as you defined in step 1. Print out the average transmission
   rates in classroom vs. in households, and the actual number of infections in classroom
   vs. in households over weekdays vs. weekends.

Generate one comprehensive R script that contains everything and meets the above
requirements. Attached is a Testing script which you should use to evaluate if you have
done the coding correctly. The key is to make sure the parameters and outcomes you
produced in your code will align with the expected values, i.e., meets.threshold is TRUE in
all cases. Make sure your code will run and pass all the tests, return the result data frames
in the test script.
