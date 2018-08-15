################################################################################################
#
#           7. Using ms and msFST to generate a neutral FST distribution                 
#
#
################################################################################################

# Perform neutral simulations with ms, and calculate FST between populations using msFST
# Calculate FST for 100,000 SNPs
# The sampling scheme (5 from Puyacatengo sulfidic, 6 from all others) is consistant with the sample size in our empirical dataset
# Population sizes (-n), migration rates (-ma), and timing of lineage joining (-ej) are based off of the best fit demographic model from dadi
# These models are included in the supplemental figures of the paper
# Use a loop to perform 1000 simulations each

for i in {1..1000}

do

ms 12 50394 -s 1 -I 2 6 6 -ma x 1.756 0.202 x -n 1 3.29 -n 2 1.49 -ej 1.076 2 1 -en 1.076 1 1 | msstats -I 2 6 6 > taco_simulations"$i".txt
ms 11 42607 -s 1 -I 2 5 6 -ma x 13.73 5.90 x -n 1 0.426 -n 2 0.514 -ej 0.121 2 1 -en 0.121 1 1 | msstats -I 2 5 6 > puya_simulations"$i".txt

done

