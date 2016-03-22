# c_domains

c_domains is a console application for micromagnetic simulations. It counts a magnetic configuration of thin ferromagnetic film deposited onto ferroelastic crystal. Especially, it takes into consideration magnetoelastic effects resulting from spontaneous strains in ferroelastic domains. 

In c_domains there are three predefined ferroelastic substrates (with surface (001)) with ID: 1 - gadolinium molybdate (GMO), 2 - lithium caesium sulphate (LCS) and 3 - potassium dihydrogen phosphate (KDP).

Application counts iteratively magnetic configuration of a thin film at equilibrium for a given temperature, substrate and ferroelastic domain pattern. It uses effective fields and it determines total free energy for this purpose. For physical details see Graczyk P. et al., Journal of Alloys and Compounds 656, 825-829 (2016).
