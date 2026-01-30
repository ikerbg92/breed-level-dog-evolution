README

This file contains all the necessary information to use the data and code to reproduce the analyses for Bargas-Galarraga, I., Vilà, C. and Gonzalez-Voyer, A. In press. High investment into reproduction is associated with reduced lifespan in dogs. American Naturalist.

The complete life history dataset can be identified by the name "data", and is provided in two versions: 

1. data (excel file) contains headers with self-explanatory column names.

2. data (.csv file) contains data without headers to easily export the dataset to R.

Both "data" files are identical except by the headers and include the following information by column (ordered left to right):

A. breed code name from Parker et al. 2017
B. breed code name from Garamszegi et al. 2020
C and D. common breed name and synonym (if any)

E to AG. Adult weight data by source listed in the supplementary material (Table S2). Adult weight data can be identified by the letter "w" at the beginning of the cell name. 
Some sources report adult weight as a range, thus numbers "1" and "2" following the letter "w" refer to the lower and upper values of the range. Column Q contains  data reported by the AKC,
which is the data that we used as "adult weight data" in the paper (as stated in methods). This is the column you would use to reproduce the analyses.

AH to AS. Lifespan data by source listed in the supplementary material (Table S1). lifespan data can be identified by the letter "l" at the beginning of the cell name. 
Some sources report lifespan as range, thus numbers "1" and "2" following the letter "l" refer to the lower and upper values of the range. Column AS contains the lifespan data used in the analyses of our article (as stated in the Methods). This is the column you would use to reproduce the analyses.

AT to BP. Litter size data by source listed in the supplementary material (Table S4). litter size data can be identified by the letters "ls" at the beginning of the cell name. 
Column BN contains litter size data used in the analyses of our article (as stated in the Methods).

BQ to CB. Neonate weight data by source listed in the supplementary material (Table S3). Neonate weight data can be identified by the letters "bw" (for "birth weight") 
at the beginning of the cell name. Column CB contains neonate weight data used in the analyses of our article (as stated in the Methods).

CC. Reproductive investment data calculated as = birth weight in kg (bwkg) multiplied by litter size (ls) = "bwkgxls".

CD. to CI. growth data from Hawthorne et al. (2004) and Posada et al. (2014). These data were used to compare with our growth estimate (CJ) as stated in methods.

CJ. Growth data calculated as stated in methods.


Source details are listed in supplementary materials*



3. code (R file): R code used for analyses. Work was done with R version 4.1.2 (2021-11-01).

As stated in the Methods section of the article, to simultaneously control for shared ancestry and hybridization of dog breeds we used SNP data from Parker et al 2017. Because we did not have data for all traits for all breeds. For linear models data is available for 92 breeds, snps.csv and hap.csv are the matrices to include as random factors in these analyses. 

4. snps (.csv file): shared ancestry data based on genomic data from Parker et al. (2017), following Garamszegi et al (2020), as stated in methods.
5. hap (.csv file): haplotype-sharing data based on genomic data from Parker et al. (2017), following Garamszegi et al (2020), as stated in methods.

For the life history traits (weight, life span, litter size and neonate weight) we had data for a variable number of breeds thus the matrices below must be used for each of the analyses to estimate the influence of shared ancestry and geneflow.

6. hapdfw (.csv file): specific haplotype sharing matrix for "dfw" object in code
7. snpsdfw (.csv file): specific snps matrix for "dfw" object in code

8. hapdfl (.csv file): specific haplotype sharing matrix for "dfl" object in code
9. snpsdfl (.csv file): specific snps matrix for "dfl" object in code

10. hapdfls (.csv file): specific haplotype sharing matrix for "dfls" object in code
11. snpsdfls (.csv file): specific snps matrix for "dfls" object in code

12. hapdfbw (.csv file): specific haplotype sharing matrix for "dfbw" object in code
13. snpsdfbw (.csv file): specific snps matrix for "dfbw" object in code


Literature cited (sources for the life history data)

American Kennel Club (AKC) (n.d.). Retrieved February 3, 2021, from https://www.akc.org/
Bell, J., K. Cavanagh, L. Tilley, and F. W. Smith. 2012. Veterinary medical guide to dog and cat breeds. CRC press 656 pp. 
Borge, K. S., R. Tønnessen, A. Nødtvedt, and A. Indrebø. 2011. Litter size at birth in purebred dogs-A retrospective study of 224 breeds. Theriogenology 75:911-919.
Clark, R. D. 2017a. Medical, Genetic & Behavioral Risk Factors of Herding Breeds. Xlibris Corporation. 743868 Library of Congress Control Number: 2017902418
Clark, R. D. 2017b Medical, Genetic & Behavioral Risk Factors of the Terrier Breeds. Xlibris Corporation. 743865 
Clark, R. D. 2017c. Medical, Genetic & Behavioral Risk Factors of the Toy Breeds. Xlibris Corporation. 743866 Library of Congress Control Number: 2017902569
Dog Breed Info (n.d.). Retrieved February 3, 2021, from https://www.dogbreedinfo.com
Easy Pet MD. (n.d.). Retrieved February 3, 2021, from http://www.easypetmd.com
Fan, R., G. Olbricht, X. Baker, and C. Hou. 2016. Birth mass is the key to understanding the negative correlation between lifespan and body size in dogs. Aging (Albany NY), 8:3209.
Gerstner, G. E., M. Cooper, and & P. Helvie. 2010. Chewing rates among domestic dog breeds. Journal of Experimental Biology 213:2266-2272.
Goleman, M., M. Karpiński, P. Czyżowski, and L. Drozd. 2015. Litter size variation in Polish selected small dog breeds. Italian Journal of Animal Science 14:3953-.
Groppetti, D., A. Pecile, C. Palestrini, S. P. Marelli, and P. Boracchi. 2017. A national census of birth weight in purebred dogs in Italy. Animals 7:43.
Groppetti, D., G. Ravasio, V. Bronzo, and A. Pecile. 2015. The role of birth weight on litter size and mortality within 24 h of life in purebred dogs: What aspects are involved? Animal reproduction science 163:112-119.
Gubbels, E. J., J. Scholten, L. Janss, and J. Rothuizen. 2009. Relationship of cryptorchidism with sex ratios and litter sizes in 12 dog breeds. Animal reproduction science 113:187-195.
Johnston, S. D., Root Kustritz, M. V., & Olson, P. S. (2001). Canine and feline theriogenology. Saunders.
Hawthorne, A. J., D. Booles, P. A.Nugent, G. Gettinby, and J. Wilkinson. 2004. Body-weight changes during growth in puppies of different breeds. The Journal of nutrition 134:2027S-2030S.
Kraus, C., S. Pavard, and D. E. Promislow. 2013. The size–life span trade-off decomposed: why large dogs die young. The American Naturalist 181:492-505.
Leroy, G., F. Phocas, B. Hedan, E. Verrier, and X. Rognon. 2015. Inbreeding impact on litter size and survival in selected canine breeds. The Veterinary Journal 203:74-78.
Michell, A. R. 1999. Longevity of British breeds of dog and its relationships with‐sex, size, cardiovascular variables and disease. Veterinary Record 145:625-629.
Mila, H., A. Grellet, A. Feugier, and S. Chastant-Maillard. 2015. Differential impact of birth weight and early growth on neonatal mortality in puppies. Journal of animal science 93:4436-4442.
Nielen, A. L. J., S. Van Der Beek, G. J.  Ubbink, and B.W. Knol. 2001. Epidemiology: Population parameters to compare dog breeds: Differences between five Dutch purebred populations. Veterinary Quarterly 23:43-49.
O’Neill, D. G., D. B. Church, P. D. McGreevy, P. C. Thomson, and D. C. Brodbelt. 2013. Longevity and mortality of owned dogs in England. The Veterinary Journal 198:638-643.
Ograk, Y. Z. 2009. Researches on litter size in Kangal breed of Turkish shepherd dogs. J. Anim Vet Adv 8:674676.
Okkens, A. C., T. W. M. Hekerman, J. W. A. De Vogel, and B. Van Haaften. 1993. Influence of litter size and breed on variation in length of gestation in the dog. Veterinary Quarterly 15:160-161.
Schrack, J., G. Dolf, I. M. Reichler, and C. Schelling. 2017. Factors influencing litter size and puppy losses in the Entlebucher Mountain dog. Theriogenology 95:163-170.
Speakman, J. R., A. Van Acker, and E.J. Harper. 2003. Age‐related changes in the metabolism and body composition of three dog breeds and their relationship to life expectancy. Aging cell 2:265-275.
Thomassen, R., G. Sanson, A. Krogenaes, J. A. Fougner, K. A. Berg, & W. Farstad. 2006. Artificial insemination with frozen semen in dogs: a retrospective study of 10 years using a non-surgical approach. Theriogenology, 66:1645-1650.
Wildt, D. E., E. J. Baas, P. K. Chakraborty, T. L. Wolfle, and A. P. Stewart. 1982. Influence of inbreeding on reproductive performance, ejaculate quality and testicular volume in the dog. Theriogenology 17:445-452.
Yilmaz, O. 2007. Turkish kangal (Karabash) shepherd dog. Impress Printing Comp. Ankara (Turkey).





