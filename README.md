# cplexNMF
Implementation of Alternating Least Square (ALS) Non-negative Matrix Factorization (NMF) algorithm in CPLEX.
The code implements the ALS NMF with a spatial neighborhood regularization according to the following equation:

![equation](http://www.sciweavers.org/tex2img.php?eq=%5Cbegin%7Baligned%7D%0A%26%20%5Cunderset%7B%5Cpmb%7BW%7D%2C%20%5Cpmb%7BH%7D%7D%7B%5Ctext%7Bminimize%7D%7D%0A%26%20%26%20f%28%5Cpmb%7BW%7D%2C%5Cpmb%7BH%7D%29%20%3D%20%5Cfrac%7B1%7D%7B2%7D%20%5C%7CV%20-%20W%5Ctimes%20H%5C%7C_F%5E2%20%2B%20%5Cfrac%7B%5Clambda%7D%7B2%7D%20%5Csum_%7Bi%2Cj%5Cin%5COmega%7D%5Csum_%7Bl%3D1%7D%5E%7Bn%7D%28h_%7Bil%7D-h_%7Bjl%7D%29%5E2%5C%5C%0A%26%20%5Ctext%7Bsubject%20to%7D%0A%26%20%26%20%5Cforall%7Ei%2Cj%3A%20W_%7Bi%2Cj%7D%2C%20H_%7Bi%2Cj%7D%20%5Cge%200%0A%5Cend%7Baligned%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

Where ![equation](http://www.sciweavers.org/tex2img.php?eq=%5COmega&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0) represent the neighborhood elements in the data Matrix V.
