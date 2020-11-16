# LAMBDA_232 = 4.934E-11#Amelin and Zaitsev, 2002.
LAMBDA_232 = 4.9475E-11 #old value, used in isoplot
ERR_LAMBDA_232 = 0.015E-11#Amelin and Zaitsev, 2002.

LAMBDA_235 = 9.8485E-10
#ERR_LAMBDA_235 =

LAMBDA_238 = 1.55125E-10
#ERR_LAMBDA_238 =

# U238_U235 = 137.817 #https://doi.org/10.1016/j.gca.2018.06.014
U238_U235 = 137.88 #old value, used in isoplot

ERR_U238_U235 = 0.031 ##https://doi.org/10.1016/j.gca.2018.06.014

lambdas = [LAMBDA_238, LAMBDA_235, LAMBDA_232]


isotope_ratios = ["U238_Pb206", "U235_Pb207", "Th232_Pb208", "Pb206_Pb207"]

minerals = ["zircon", "baddeleyite", "perovskite", "monazite", "apatite"]

EarthAge = 4600

sqrt2pi = 2.506628274631

listColumnNames = ['Th232/Pb208', 'err.Th232/Pb208',
                       '206Pb/207Pb', 'err.206Pb/207Pb',
                       '235U/207Pb', 'err.235U/207Pb',
                       '238U/206Pb', 'err.238U/206Pb',
                       'corr.',
                       'Age_232Th/208Pb', 'err.Age_232Th/208Pb',
                       'Age_206Pb/207Pb', 'err.Age_206Pb/207Pb',
                       'Age_238U/206Pb', 'err.Age_238U/206Pb',
                       'Age_235U/207Pb', 'err.Age_235U/207Pb',
                       '%disc.67-86', '%disc.86-57', 'is grain good?']

