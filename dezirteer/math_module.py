# -*- coding: utf-8 -*-
from bisect import *
from const import *
import operator
import functools
from math import sqrt, log, fabs
import scipy.stats
#from scipy import std, exp
from numpy import std, exp
import random

pbpb_table =  []
concordia_table = []


def calc_rho(rat68, rat68err, rat75, rat75err, rat76, rat76err):
    any_error = False
    corr_coef_75_68 = (rat68err / rat68) / (rat75err / rat75)
    corr_coef_86_76 = (rat68err / rat68) * (rat76 / rat76err)
    if corr_coef_75_68 > 1:
        corr_coef_75_68 = 1 / corr_coef_75_68
        any_error = True
    if corr_coef_86_76 > 1:
        corr_coef_86_76 = 1 / corr_coef_86_76
        any_error = True
    return corr_coef_75_68, corr_coef_86_76, any_error


def t_student(alpha, gl):
    return scipy.stats.t.ppf(1 - (alpha / 2), gl)


def calc_ratio(age):
    # print(age)
    pb207_u235 = exp(LAMBDA_235 * age * 1000000) - 1
    pb206_u238 = exp(LAMBDA_238 * age * 1000000) - 1
    if age != 0:
        pb207_pb206 = (pb207_u235 / pb206_u238) * (1 / U238_U235)
        u238_pb206 = 1 / pb206_u238
    else:
        pb207_pb206 = 0
        u238_pb206 = 1000
    pb208_th232 = exp(LAMBDA_232 * age * 1000000) - 1
    return [pb206_u238, pb207_u235, pb207_pb206, u238_pb206, pb208_th232]


def compb(age, n):  # Stacey & Cramers 2 stage pb evolution model
    r4570=calc_ratio(4570)
    r3700=calc_ratio(3700)
    rage=calc_ratio(age)
##    if n == 3:
##        return compb(age, 1) / compb(age, 0)  # 7/6c
##    elif n == 4:
##        return compb(age, 2) / compb(age, 0)  # 8/6c
##    else: 
    if age <= 3700:
        if n == 0:
            return 11.152 + 9.735 * (r3700[0] - rage[0])  # 6/4c
        elif n == 1:
            return 12.998 + 9.735 / U238_U235 * (r3700[1] - rage[1])  # 7/4c
        elif n == 2:
            return 31.23 + 36.837 * (r3700[4] - rage[4])  # 8/4c
    else:
        if n == 0:
            return 9.307 + 7.1925 * (r4570[0] - rage[0])  # 6/4c
        elif n == 1:
            return 10.294 + 7.192 / U238_U235 * (r4570[1] - rage[1])  # 7/4c
        elif n == 2:
            return 29.476 + 32.208 * (r4570[4] - rage[4])  # 8/4c


def pb4cor(fc, pb_pb4, pb_uth, lam):  # universal 204pb corr function for all measured ratios
    if pb_pb4[0] > 0 and pb_uth[0] > 0:
        # fc = pb_pb4_c / pb_pb4[0]
        ratio = pb_uth[0]*(1-fc)
        tmp1_int = (pb_pb4[1] / pb_pb4[0]) ** 2
        tmp2_int = (pb_uth[1] / pb_uth[0]) ** 2
        ratio_err_int = sqrt(tmp1_int / (pb_pb4[0] * (1 - fc)) ** 2 + tmp2_int) * ratio
        if pb_pb4[2] != -1 and pb_uth[2] != -1:
            tmp1_prop = (pb_pb4[2] / pb_pb4[0]) ** 2
            tmp2_prop = (pb_uth[2] / pb_uth[0]) ** 2
            ratio_err_prop = sqrt((sqrt(tmp1_prop) / pb_pb4[0] / (1 - fc)) ** 2 + tmp2_prop) * ratio
            age_err_prop = ratio_err_prop / (1 + ratio) / lam / 1000000
        if ratio > 0:
            age = log(1 + ratio) / lam / 1000000
            age_err_int = ratio_err_int / (1 + ratio) / lam / 1000000
        else:
            age = -1
            age_err_int = -1
    else:
        age = -1
        age_err_int = -1
        age_err_prop = -1
        ratio = -1
        ratio_err_int = -1
        ratio_err_prop = -1
    return age, age_err_int, age_err_prop, ratio, ratio_err_int, ratio_err_prop

def andersen(xt2, yt2, zt2, c7, c8, x, y, z, u):
    k = U238_U235
    dt = 1
    a = 0
    b = 4500
    c = 1000

    while dt > 0.001: # solution of equation 7 (Andersen, 2002)
        if eq7(c, xt2, yt2, zt2, c7, c8, x, y, z, u, k) < 0:
            b = c
        else:
            a = c
        dt = b - a
        c = dt / 2 + a
    t1 = a
    t1_ratios = calc_ratio(t1)
    xt1 = t1_ratios[1]
    yt1 = t1_ratios[0]
    # zt1 = eq7(t1, xt2, yt2, zt2, c7, c8, x, y, z, u, k)[3]
    fc = (-y * xt1 + y * xt2 + yt2 * xt1 + x * yt1 - x * yt2 - xt2 * yt1) / (
                -y * xt1 + y * xt2 + y * c7 * k * yt1 - y * c7 * k * yt2) * 100
    yr = y * (1 - fc)
    xr = x - y * c7 * k * fc
    zr = z - y * c8 * u * fc
    fl = (yt1 - yr) / (yt1 - yt2)
    return t1, fc, xr, yr, zr, fl  # corrected age, fract. of common lead, radiogenic ratios and fratc of pb loss


def eq7(t1, xt2, yt2, zt2, c7, c8, x, y, z, u, k):
    xt1 = calc_ratio(t1)[1]
    yt1 = calc_ratio(t1)[0]
    zt1 = calc_ratio(t1)[4]
    zero = (y * (xt1 - xt2) - yt2 * xt1 + x * (yt2 - yt1) + xt2 * yt1) / (
            xt1 - xt2 - c7 * k * yt1 + c7 * k * yt2) - (
                   z * (yt2 - yt1) + zt2 * yt1 + y * (zt1 - zt2) - yt2 * zt1) / (
                   zt1 - zt2 - c8 * u * yt1 + c8 * u * yt2)
    # zero = math.fabs(zero)
    return zero

def pbc_corr(zir, corr_type, *args):  # returns Pbc-corrected ages
    d = 1
    corr_age = [-1, -1, -1]
    mr68 = zir.pb206_u238
    mr75 = zir.pb207_u235
    mr82 = zir.pb208_th232
    mr76 = zir.pb207_pb206
    mr64 = zir.pb206_pb204
    mr74 = zir.pb207_pb204
    mr84 = zir.pb208_pb204
    mr28 = zir.th232_u238 #[0.5, 0.05, 0.1]
    mr86 = mr82[0] * mr28[0] / mr68[0]
    #mru84 = zir.u238_pb204
    #mrth24 = zir.th232_pb204
    age = zir.calc_age(0)[0]
    com64 = compb(age, 0)
    com74 = compb(age, 1)
    com84 = compb(age, 2)
    com76 = com74 / com64
    com86 = com84 / com64
    if corr_type == 1:  # 204
        if mr64[0] != -1 and mr68[0] != -1:
            a68 = pb4cor(com64/mr64[0], mr64, mr68, LAMBDA_238)
        else:
            a68 = [-1, -1, -1, -1, -1, -1]
        if mr74[0] != -1 and mr75[0] != -1:
            a75 = pb4cor(com76/mr76[0]*com64/mr64[0], mr74, mr75, LAMBDA_235)
        else:
            a75 = [-1, -1, -1, -1, -1, -1]

        if mr84[0] != -1 and mr82[0] != -1:
            a82 = pb4cor(com86/mr86*com64/mr64[0], mr84, mr82, LAMBDA_232)
        else:
            a82 = [-1, -1, -1, -1, -1, -1]
        a76 = [-1, -1, -1, -1, -1, -1]
        a76[3] = a75[3] / a68[3] / U238_U235
        a76[0] = find_age(a76[3])
        if a68[0] > 0 and a75[0] > 0 and a76[0] > 0:
            t = a76[0] * 1000000
            e5 = exp(LAMBDA_235 * t)
            e8 = exp(LAMBDA_238 * t)
            a76[4] = sqrt((a75[4] * a68[2]) ** 2 + (a68[4] * a75[2]) ** 2) / (U238_U235 * a68[3] ** 2)
            a76[1] = U238_U235 * a76[4] / (LAMBDA_235 * e5/(e8 - 1) - LAMBDA_238 * e8 * (e5 - 1) / (e8 - 1)**2)
            a76[1] = a76[1] / 1000000

            a76[5] = sqrt((a75[5] * a68[2]) ** 2 + (a68[5] * a75[2]) ** 2) / (U238_U235 * a68[3] ** 2)
            a76[2] = U238_U235 * a76[5] / (
                        LAMBDA_235 * e5 / (e8 - 1) - LAMBDA_238 * e8 * (e5 - 1) / (e8 - 1) ** 2)
            a76[2] = a76[2] / 1000000
        else:
            a76 = [-1, -1, -1, -1, -1, -1]
        
        # a4c = [a68[0], a75[0], a76[0], a82[0]]
        # a4c_err_int = [a68[1], a75[1], a76[1], a82[1]]
        # a4c_err_prop = [a68[2], a75[2], a76[2], a82[2]]
        # corr_age = [a4c, a4c_err_int, a4c_err_prop]
        corr_age = [a68[0], a68[1], a68[2]]
        # print(corr_age)
    elif corr_type == 2 and mr76[0] > 0 and mr68[0] > 0:  # 207
        a = 0
        b = 4500
        c = 1000
        # age
        while d > 0.001:
            r = calc_ratio(c)
            x1 = r[3]
            x2 = (r[2] - com76) / (mr68[0] * (mr76[0] - com76))
            if x2 > x1:
                b = c
            else:
                a = c
            c = (b - a) / 2 + a
            d = b - a
        corr_age[0] = a

        # error
        e1 = exp(LAMBDA_238 * a) - 1
        e2 = exp(LAMBDA_235 * a) - 1
        rv = U238_U235 ** 2 * ((mr68[0] * mr76[1]) ** 2 + (mr76[0] * mr68[1]) ** 2)
        d = (U238_U235 * com76 * LAMBDA_238 * (e1 + 1) - LAMBDA_235 * (e2 + 1)) ** 2
        n1 = 0  # n1=(U238_U235*(zir.pb206_u238[0]-calc_ratio(t)[0])*) #commonly it's assumed r76c_err=0 => n1=0
        n2 = U238_U235 ** 2 * com76 * (com76 - 2 * mr76[0]) * mr68[1] ** 2
        n = n1 + n2 + rv
        corr_age[1] = sqrt(n / d) / 1000000

        if mr68[2] != -1 and mr76[2] != -1:
            rv = U238_U235 ** 2 * ((mr68[0] * mr76[2]) ** 2 + (mr76[0] * mr68[2]) ** 2)
            n2 = U238_U235 ** 2 * com76 * (com76 - 2 * mr76[0]) * mr68[2] ** 2
            n = n1 + n2 + rv
            corr_age[2] = sqrt(n / d) / 1000000
        else:
            corr_age[2] = -1
        
    elif corr_type == 3 and mr82[0] > 0 and mr68[0] > 0 and mr28[0] > 0:  # 208
        t = 1000
        # age
        a = 0
        b = 4500
        c = 1000
        while d > 0.001:
            r = calc_ratio(c)
            x1 = r[3]
            x2 = (r[3] * r[4] * mr28[0] - com86) / (mr68[0] * (mr86 - com86))
            if x2 > x1:
                b = c
            else:
                a = c
            c = (b - a) / 2 + a
            d = b - a

        corr_age[0] = a
        
        # error
        e1 = exp(LAMBDA_238 * a * 10**6)
        e2 = exp(LAMBDA_232 * a * 10**6)
        c1 = mr82[0] + 1 - e2
        c2 = mr28[0] * com64 / com84
        d = (c2 * LAMBDA_232 * e2 - LAMBDA_238 * e1) ** 2
        n1 = c1 ** 2 * (mr28[1]/com86) ** 2
        n2 = c2 * mr82[1] + mr68[1] ** 2
        n = n1 + n2
        corr_age[1] = sqrt(n / d) / 1000000
        
        if mr68[2] != -1 and mr76[2] != -1 and mr82[2] != -1:
            n1 = c1 ** 2 * (com64 / com84 * mr28[2]) ** 2
            n2 = c2 * mr82[2] + mr68[2] ** 2
            n = n1 + n2
            corr_age[2] = sqrt(n / d) / 1000000
        else:
            corr_age[2] = -1

    elif corr_type == 4:  # and
        # age
        t2 = 0  # NEED CORRECTION!!!  age of pb lost, must entered by user
        t2_ratios = calc_ratio(t2)
        xt2 = t2_ratios[1]
        yt2 = t2_ratios[0]
        zt2 = t2_ratios[4]
        c7 = 15.628 / 18.7
        c8 = 38.63 / 18.7
        rho = zir.corr_coef_75_68
        x = mr75[0]
        y = mr68[0]
        z = mr82[0]
        u = 1/mr28[0]

        # print(x,y,z,eq7(500, xt2, yt2, zt2, c7, c8, x, y, z, u, U238_U235))
        corr_data = andersen(xt2, yt2, zt2, c7, c8, x, y, z, u)
        corr_age[0] = corr_data[0]

        # mc_fc = []
        mc_ages_int = []
        mc_ages_prop = []
        # print(x,y,z,corr_age[0])
        # sigma errors
        for i in range(1000): # Monte-Carlo statistics collection
            mcx = random.normalvariate(0, 1)
            mcy = mcx * rho + sqrt(1 - rho**2) * random.normalvariate(0, 1)
            mcz = random.normalvariate(0, 1)
            mcx_int = mcx * mr75[1] + x
            mcy_int = mcy * mr68[1] + y
            mcz_int = mcz * mr82[1] + z
            mc_ages_int.append(andersen(xt2, yt2, zt2, c7, c8, mcx_int, mcy_int, mcz_int, u)[0])

            if mr75[2] != -1 and mr68[2] != -1 and mr82[2] != -1:
                mcx_prop = mcx * mr75[2] + x
                mcy_prop = mcy * mr68[2] + y
                mcz_prop = mcz * mr82[2] + z
                mc_ages_prop.append(andersen(xt2, yt2, zt2, c7, c8, mcx_prop, mcy_prop, mcz_prop, u)[0])
            else:
                mc_ages_prop = -1
            # mc_fc.append(andersen(xt2, yt2, zt2, c7, c8, mkx, mky, mkz, u)[1])

        corr_age[1] = std(mc_ages_int)
        corr_age[2] = std(mc_ages_prop)

    else:
        corr_age = [-1, -1, -1]
    return corr_age

def sumproduct(*lists):
    return sum(functools.reduce(operator.mul, data) for data in zip(*lists))


def fill_pbpb_table():  # fills the table of 207pb-206pb ages
    if len(pbpb_table) == 0:
        i = 5
        while i <= 4600 * 5:
            age = (1 / U238_U235) * (exp(LAMBDA_235 * i * 1000000) - 1) / (exp(LAMBDA_238 * i * 1000000) - 1)
            pbpb_table.append(age)
            i += 5


def fill_concordia_table():
    if len(concordia_table) == 0:
        t = 1
        ratios = [-1, -1, -1]
        while t <= EarthAge:
            mln_age = t * 1000000
            e8 = exp(LAMBDA_238 * mln_age) - 1
            e5 = exp(LAMBDA_235 * mln_age) - 1
            ratios.append(e8)
            ratios.append(e5)
            pbpb_age = (1 / U238_U235) * e5 / e8
            ratios.append(pbpb_age)
            concordia_table.append(ratios)
            ratios = []
            t += 1


# calculates pb-pb age from the table of pb-pb ratios,filled in fill_Pb206207_table()
def find_age(pLeadRatio):
    return bisect(pbpb_table, pLeadRatio) * 5


class Filters(object):  # describes filters that should be applied to data in Analysis_set object
    def __init__(self, filter_by_uconc=[False, 1000], which_age=[0, 1000], use_pbc=0,
                 filter_by_err=[False, 0.1], include207235Err=False,
                 pos_disc_filter=0.2, neg_disc_filter=-0.1, disc_type=[4, 1000],
                 sample_name_filter=[], unc_type='1', filter_by_commPb=[False, 0.1], minAgeCrop=0, maxAgeCrop=EarthAge):
        self.__filter_by_uconc = filter_by_uconc
        self.__which_age = which_age
        self.__use_pbc = use_pbc
        self.__filter_by_err = filter_by_err
        self.__include207235Err = include207235Err
        self.__pos_disc_filter = pos_disc_filter
        self.__neg_disc_filter = neg_disc_filter
        self.__disc_type = disc_type
        self.__sample_name_filter = sample_name_filter
        self.__unc_type = unc_type #1 for internal, 2 for propagated
        self.__filter_by_commPb = filter_by_commPb
        self.__minAgeCrop = minAgeCrop
        self.__maxAgeCrop = maxAgeCrop


    @property
    def filter_by_uconc(self):
        return self.__filter_by_uconc

    @filter_by_uconc.setter
    def filter_by_uconc(self, value):
        self.__filter_by_uconc = value

    @property
    def which_age(self):
        return self.__which_age

    @which_age.setter
    def which_age(self, value):
        self.__which_age = value

    @property
    def use_pbc(self):
        return self.__use_pbc

    @use_pbc.setter
    def use_pbc(self, value):
        self.__use_pbc = value

    @property
    def filter_by_err(self):
        return self.__filter_by_err

    @filter_by_err.setter
    def filter_by_err(self, value):
        self.__filter_by_err = value

    @property
    def include207235Err(self):
        return self.__include207235Err

    @include207235Err.setter
    def include207235Err(self, value):
        self.__include207235Err = value

    @property
    def pos_disc_filter(self):
        return self.__pos_disc_filter

    @pos_disc_filter.setter
    def pos_disc_filter(self, value):
        self.__pos_disc_filter = value

    @property
    def neg_disc_filter(self):
        return self.__neg_disc_filter

    @neg_disc_filter.setter
    def neg_disc_filter(self, value):
        self.__neg_disc_filter = value

    @property
    def disc_type(self):
        return self.__disc_type

    @disc_type.setter
    def disc_type(self, value):
        self.__disc_type = value

    @property
    def sample_name_filter(self):
        return self.__sample_name_filter

    @sample_name_filter.setter
    def sample_name_filter(self, value):
        self.__sample_name_filter = value

    @property
    def unc_type(self):
        return self.__unc_type

    @unc_type.setter
    def unc_type(self, value):
        self.__unc_type = value

    @property
    def filter_by_commPb(self):
        return self.__filter_by_commPb

    @filter_by_commPb.setter
    def filter_by_commPb(self, value):
        self.__filter_by_commPb = value

    @property
    def minAgeCrop(self):
        return self.__minAgeCrop

    @minAgeCrop.setter
    def minAgeCrop(self, value):
        self.__minAgeCrop = value

    @property
    def maxAgeCrop(self):
        return self.__maxAgeCrop

    @maxAgeCrop.setter
    def maxAgeCrop(self, value):
        self.__maxAgeCrop = value


class Analysis(object):
    def __init__(self, analysis_name="", exposure_time="",
                 pb206_u238=(0, 0, 0), pb207_u235=(0, 0, 0), corr_coef_75_68=0, corr_coef_86_76=0,
                 pb208_th232=(0, 0, 0), pb207_pb206=(0, 0, 0), u_conc=(0, 0, 0), pbc=(0, 0, 0), pb206_pb204=(0, 0, 0),
                 pb207_pb204=(0, 0, 0), pb208_pb204=(0, 0, 0), th232_pb204=(0, 0, 0), u238_pb204=(0, 0, 0),
                 sigma_level=0, final_U_Th_Ratio = (0, 0, 0), th232_u238=(0, 0, 0)):
        self.__analysis_name = analysis_name
        self.__exposure_time = exposure_time
        self.__pb206_u238 = pb206_u238
        self.__pb207_u235 = pb207_u235
        self.__corr_coef_75_68 = corr_coef_75_68
        self.__corr_coef_86_76 = corr_coef_86_76
        self.__pb208_th232 = pb208_th232
        self.__pb207_pb206 = pb207_pb206
        self.__u_conc = u_conc
        self.__pbc = pbc
        self.__pb206_pb204 = pb206_pb204
        self.__pb207_pb204 = pb207_pb204
        self.__pb208_pb204 = pb208_pb204
        self.__th232_pb204 = th232_pb204
        self.__u238_pb204 = u238_pb204
        self.__sigma_level = sigma_level
        self.__final_U_Th_Ratio = final_U_Th_Ratio
        self.__th232_u238 = th232_u238

    def __repr__(self):
        return self.analysis_name

    @property
    def analysis_name(self):
        return self.__analysis_name

    @analysis_name.setter
    def analysis_name(self, value):
        self.__analysis_name = value

    @property
    def exposure_time(self):
        return self.__exposure_time

    @exposure_time.setter
    def exposure_time(self, value):
        self.__exposure_time = value

    @property
    def pb206_u238(self):
        return self.__pb206_u238

    @pb206_u238.setter
    def pb206_u238(self, value):
        self.__pb206_u238 = value

    @property
    def pb207_u235(self):
        return self.__pb207_u235

    @pb207_u235.setter
    def pb207_u235(self, value):
        self.__pb207_u235 = value

    @property
    def corr_coef_75_68(self):
        return self.__corr_coef_75_68

    @corr_coef_75_68.setter
    def corr_coef_75_68(self, value):
        self.__corr_coef_75_68 = value

    @property
    def corr_coef_86_76(self):
        return self.__corr_coef_86_76

    @corr_coef_86_76.setter
    def corr_coef_86_76(self, value):
        self.__corr_coef_86_76 = value

    @property
    def pb208_th232(self):
        return self.__pb208_th232

    @pb208_th232.setter
    def pb208_th232(self, value):
        self.__pb208_th232 = value

    @property
    def pb207_pb206(self):
        return self.__pb207_pb206

    @pb207_pb206.setter
    def pb207_pb206(self, value):
        self.__pb207_pb206 = value

    @property
    def u_conc(self):
        return self.__u_conc

    @u_conc.setter
    def u_conc(self, value):
        self.__u_conc = value

    @property
    def pbc(self):
        return self.__pbc

    @pbc.setter
    def pbc(self, value):
        self.__pbc = value

    @property
    def pb206_pb204(self):
        return self.__pb206_pb204

    @pb206_pb204.setter
    def pb206_pb204(self, value):
        self.__pb206_pb204 = value

    @property
    def pb207_pb204(self):
        return self.__pb207_pb204

    @pb207_pb204.setter
    def pb207_pb204(self, value):
        self.__pb207_pb204 = value

    @property
    def pb208_pb204(self):
        return self.__pb208_pb204

    @pb208_pb204.setter
    def pb208_pb204(self, value):
        self.__pb208_pb204 = value

    @property
    def th232_pb204(self):
        return self.__th232_pb204

    @th232_pb204.setter
    def th232_pb204(self, value):
        self.__th232_pb204 = value

    @property
    def u238_pb204(self):
        return self.__u238_pb204

    @u238_pb204.setter
    def u238_pb204(self, value):
        self.__u238_pb204 = value

    @property
    def sigma_level(self):
        return self.__sigma_level

    @sigma_level.setter
    def sigma_level(self, value):
        self.__sigma_level = value

    @property
    def final_U_Th_Ratio(self):
        return self.__final_U_Th_Ratio

    @final_U_Th_Ratio.setter
    def final_U_Th_Ratio(self, value):
        self.__final_U_Th_Ratio = value

    @property
    def th232_u238(self):
        return self.__th232_u238

    @th232_u238.setter
    def th232_u238(self, value):
        self.__th232_u238 = value

    def u238_pb206(self):
        rat238206 = 1 / self.pb206_u238[0]
        return [rat238206, rat238206 * (self.pb206_u238[1] / self.pb206_u238[0]),
                rat238206 * (self.pb206_u238[2] / self.pb206_u238[0])]

    # calculates age ± error from isotopic value and uncertainty
    def calc_age(self, isotopic_system):
        age_err_int = -1
        age_err_prop = -1
        try:
            if isotopic_system == 0 and self.pb206_u238[0] > 0:
                age = (1 / lambdas[isotopic_system]) * log(self.pb206_u238[0] + 1) / 1000000
                age_err_int = (1 / lambdas[isotopic_system]) * self.pb206_u238[1] / 1000000
                age_err_prop = (1 / lambdas[isotopic_system]) * self.pb206_u238[2] / 1000000

            elif isotopic_system == 1 and self.pb207_u235[0] > 0:
                age = (1 / lambdas[isotopic_system]) * log(self.pb207_u235[0] + 1) / 1000000
                age_err_int = (1 / lambdas[isotopic_system]) * self.pb207_u235[1] / 1000000
                age_err_prop = (1 / lambdas[isotopic_system]) * self.pb207_u235[2] / 1000000

            elif isotopic_system == 2 and self.pb208_th232[0] > 0:
                age = (1 / lambdas[isotopic_system]) * log(self.pb208_th232[0] + 1) / 1000000
                age_err_int = (1 / lambdas[isotopic_system]) * self.pb208_th232[1] / 1000000
                age_err_prop = (1 / lambdas[isotopic_system]) * self.pb208_th232[2] / 1000000

            elif isotopic_system == 3 and self.pb207_pb206[0] > .04605:
                # .04605 corresponds to age67 ~ 0
                age = find_age(self.pb207_pb206[0]) * 1000000
                C1 = 1 / U238_U235
                C2 = LAMBDA_235
                C3 = LAMBDA_238
                df_int = self.pb207_pb206[1]
                df_prop = self.pb207_pb206[2]
                dfdt = C1 * (C3 * exp(C3 * age) * (exp(C2 * age) - 1) - C2 * exp(C2 * age) *
                             (exp(C3 * age) - 1)) / ((exp(C3 * age) - 1) ** 2)
                age_err_int = abs(df_int / dfdt / 1000000)
                age_err_prop = abs(df_prop / dfdt / 1000000)
                age = age / 1000000
            else:
                age = -1
            return [age, age_err_int, age_err_prop]
        except ValueError:
            pass

    # calculates two type of discordance: (1) between 206_238 and 207_235 and (2) between 206_238 and 206_207
    def calc_discordance(self, disc_type, age_cutoff): #0: u-conc, #1 - fixed limit, #2 - 67-86, #3- 57-86 #4 - the one with the lesser value
        age_206_238 = self.calc_age(0)[0]
        age_207_235 = self.calc_age(1)[0]
        age_207_206 = self.calc_age(3)[0]
        disc_68_57 = age_207_235 / age_206_238 - 1
        disc_68_76 = 1 - age_206_238 / age_207_206

        if disc_type == 0:
            return -1

        if disc_type == 1:
            if age_206_238 > age_cutoff:
                return disc_68_57
            else:
                return disc_68_76

        if disc_type == 2:
            return disc_68_76

        elif disc_type == 3:
            return disc_68_57

        else: #if disc_type[0] == 4
            if abs(disc_68_57) < abs(disc_68_76):
                return disc_68_57
            else:
                return disc_68_76

    # true if should use U238/Pb206, false if should use Pb206/Pb207
    def use_206_238(self, age_cutoff):
        return self.calc_age(0)[0] < age_cutoff

    # checks if a grain passes user-defined Filters
    def is_grain_good(self, pFilter):
        do_uconc = pFilter.filter_by_uconc[0]
        uconc_ppm_cutoff = pFilter.filter_by_uconc[1]
        which_age = pFilter.which_age[0]
        age_fixed_limit = pFilter.which_age[1]
        do_err = pFilter.filter_by_err[0]
        do_207235_err = pFilter.include207235Err
        err_cutoff = pFilter.filter_by_err[1]
        type_disc = pFilter.disc_type[0]
        pos_disc_cutoff = pFilter.pos_disc_filter
        neg_disc_cutoff = pFilter.neg_disc_filter
        sample_name_filter = pFilter.sample_name_filter
        min_age_crop = pFilter.minAgeCrop
        max_age_crop = pFilter.maxAgeCrop
        is_grain_in_chosen_sample = True
        is_age_good = True
        age_206_238 = self.calc_age(0)
        age_207_206 = self.calc_age(3)

        # cut out negative ratios
        if self.pb206_u238[0] < 0 or self.pb207_u235[0] < 0 or self.pb207_pb206[0] < 0:
            are_ratios_positive = False
        else:
            are_ratios_positive = True

        # filter by Uconc
        if do_uconc and (self.u_conc[0] > uconc_ppm_cutoff):
            is_uconc_good = False
        else:
            is_uconc_good = True

        # filter by sample name
        if sample_name_filter != []:
            if parse_sample_analysis(self.analysis_name)[0] in sample_name_filter:#str(self.analysis_name).rpartition(p_divider)[0] in sample_name_filter:
                is_grain_in_chosen_sample = True
            else:
                is_grain_in_chosen_sample = False

        # decide on the default age system
        if which_age == 0:  # from the lesser error
            if age_206_238[int(pFilter.unc_type)] > age_207_206[int(pFilter.unc_type)]:
                age_68_67 = 1
                this_age = 3
            else:
                age_68_67 = 0
                this_age = 0
        elif (which_age == 1 and age_206_238[0] > age_fixed_limit) or which_age == 2:  # fixed limit, age>limit
            age_68_67 = 1
            this_age = 3
        elif (which_age == 1 and age_206_238[0] < age_fixed_limit) or which_age == 3:  # fixed limit, age<limit
            age_68_67 = 0
            this_age = 0
        else:
            age_68_67 = -1

        if age_68_67 == 0:
            age = age_206_238[0]
            err = age_206_238[int(pFilter.unc_type)]
        else:
            age = age_207_206[0]
            err = age_207_206[int(pFilter.unc_type)]

        # filter by minimum and maximum age crops
        if min_age_crop > 0: # TEMP TO BE DELETED
            pass

        if (age > min_age_crop) and (age < max_age_crop):
            is_age_good = True
        else:
            is_age_good = False

        # filter by measurement error
        if do_err:
            if err / age < err_cutoff:  # checks whether error is within limit
                is_err_good = True
            else:
                is_err_good = False
            if do_207235_err == 1:
                age207235 = self.calc_age(1)[0]
                age207235err = self.calc_age(1)[int(pFilter.unc_type)]
                if age207235err / age207235 < err_cutoff:
                    is_207235err_good = True
                else:
                    is_207235err_good = False
            else:
                is_207235err_good = True  # if no need to filter by 207-235 age error
        else:
            is_err_good = True  # if no need to filter by age error
            is_207235err_good = True

        # filter by discordance #0: u-conc, #1 - fixed limit, #2 - 57-86, #3- 67-86 #4 - the one with the lesser value
        disc = self.calc_discordance(type_disc, pFilter.disc_type[1]) #22222
        '''if type_disc == 0:
            disc = -1  # this is for future implementation of 'based on Uconc' algorithm
        elif type_disc == 1:
            disc = self.calc_discordance(1, )
        elif type_disc == 2:
            disc = self.calc_discordance(2)
        elif type_disc == 3:
            disc = self.calc_discordance(3)  #11111
        else: #if 4
            disc = self.calc_discordance(4)  # 11111'''

        if (disc < pos_disc_cutoff) & (disc > neg_disc_cutoff):
            is_disc_good = True
        else:
            is_disc_good = False

        return are_ratios_positive & is_uconc_good & is_err_good & is_207235err_good & is_disc_good & \
               is_grain_in_chosen_sample & is_age_good, this_age


class AnalysesSet(object):
    def __init__(self, analyses_list, name):
        self.__analyses_list = analyses_list
        self.__name = name
        self.__good_set = {}
        self.__bad_set = []
        self.__min_206_238 = 0
        self.__max_206_238 = 0
        self.__min_207_235 = 0
        self.__max_207_235 = 0
        self.__max_238_206 = 0
        self.__min_238_206 = 0
        self.__min_207_206 = 0
        self.__max_207_206 = 0
        self.__min_age = 0
        self.__max_age = 4500

    def __repr__(self):
        return str(self.analyses_list)

    def __len__(self):
        return len(self.analyses_list)

    @property
    def min_age(self):
        return self.__min_age

    @min_age.setter
    def min_age(self, value):
        self.__min_age = value

    @property
    def max_age(self):
        return self.__max_age

    @max_age.setter
    def max_age(self, value):
        self.__max_age = value

    @property
    def min_206_238(self):
        return self.__min_206_238

    @min_206_238.setter
    def min_206_238(self, value):
        self.__min_206_238 = value

    @property
    def max_206_238(self):
        return self.__max_206_238

    @max_206_238.setter
    def max_206_238(self, value):
        self.__max_206_238 = value

    @property
    def min_207_235(self):
        return self.__min_207_235

    @min_207_235.setter
    def min_207_235(self, value):
        self.__min_207_235 = value

    @property
    def max_207_235(self):
        return self.__max_207_235

    @max_207_235.setter
    def max_207_235(self, value):
        self.__max_207_235 = value

    @property
    def min_238_206(self):
        return self.__min_238_206

    @min_238_206.setter
    def min_238_206(self, value):
        self.__min_238_206 = value

    @property
    def max_238_206(self):
        return self.__max_238_206

    @max_238_206.setter
    def max_238_206(self, value):
        self.__max_238_206 = value

    @property
    def min_207_206(self):
        return self.__min_207_206

    @min_207_206.setter
    def min_207_206(self, value):
        self.__min_207_206 = value

    @property
    def max_207_206(self):
        return self.__max_207_206

    @max_207_206.setter
    def max_207_206(self, value):
        self.__max_207_206 = value

    @property
    def name(self):
        return self.__name

    @name.setter
    def name(self, value):
        self.__name = value

    @property
    def analyses_list(self):
        return self.__analyses_list

    @analyses_list.setter
    def analyses_list(self, value):
        self.__analyses_list = value

    @property
    def good_set(self):
        return self.__good_set

    @property
    def bad_set(self):
        return self.__bad_set

    @property
    def filters(self):
        return self.__filters



    # sorts data into good and bad sets depending on Filters settings. Returns several parameters of the good set:
    # number of grains, weighted average age ± uncertainty (±1s and 95%), MSWD, max and min ages, max and min conc values
    def good_bad_sets(self, p_filter):
        index = 0
        max_age = 0
        min_age = 5000

        max_206_238 = 0
        min_206_238 = 2

        max_207_235 = 0
        min_207_235 = 100

        max_238_206 = 0
        min_238_206 = 2000

        max_207_206 = 0
        min_207_206 = 1


        ages = []
        errs_inv_sq = []
        wa_age = 0
        aux_list = []
        wa_age_err = 0
        wa_age_err_scatter = 0
        mswd = 0
        sumproduct_ages_err = 0
        number_of_good_grains = 0
        sum_of_inv_sq_err = 0

        k = 3

        self.good_set.clear()
        self.bad_set.clear()
        while index <= len(self.analyses_list) - 1:
            zircon = self.analyses_list[index]

            #parsed_analysis = parse_sample_analysis(str(zircon))
            l_is_grain_good = Analysis.is_grain_good(zircon, p_filter)
            if l_is_grain_good[0]:
                z_age = zircon.calc_age(l_is_grain_good[1])
                z_206_238 = zircon.pb206_u238
                z_207_235 = zircon.pb207_u235
                z_238_206 = zircon.u238_pb206()
                z_207_206 = zircon.pb207_pb206

                if z_age[0] > max_age:
                    max_age = int(z_age[0])
                if z_age[0] < min_age:
                    min_age = int(z_age[0])

                if z_206_238[0] > max_206_238:
                    max_206_238 = z_206_238[0]+k*z_206_238[1]
                if z_206_238[0] < min_206_238:
                    min_206_238 = z_206_238[0]-k*z_206_238[1]

                if z_207_235[0] > max_207_235:
                    max_207_235 = z_207_235[0]+k*z_207_235[1]
                if z_207_235[0] < min_207_235:
                    min_207_235 = z_207_235[0]-k*z_207_235[1]

                if z_238_206[0] > max_238_206:
                    max_238_206 = z_238_206[0]+k*z_238_206[1]
                if z_238_206[0] < min_238_206:
                    min_238_206 = z_238_206[0]-k*z_238_206[1]

                if z_207_206[0] > max_207_206:
                    max_207_206 = z_207_206[0]+k*z_207_206[1]
                if z_207_206[0] < min_207_206:
                    min_207_206 = z_207_206[0]-k*z_207_206[1]

                self.__good_set.update({zircon: z_age})
                number_of_good_grains += 1
                ages.append(z_age[0])
                errs_inv_sq.append(1 / (z_age[1] ** 2))
            else:
                self.__bad_set.append(zircon)
            index += 1

        if number_of_good_grains > 1:
            wa_age = sumproduct(ages, errs_inv_sq) / sum(errs_inv_sq)
            wa_age_err = sqrt(1 / sum(errs_inv_sq))
            aux_list = [(x - wa_age) ** 2 for x in ages]
            mswd = sumproduct(errs_inv_sq, aux_list) / (number_of_good_grains - 1)
        elif number_of_good_grains == 1:
            wa_age = ages[0]
            wa_age_err = 0
            mswd = 0
        elif number_of_good_grains == 0:
            wa_age = 0
            wa_age_err = 0
            mswd = 0
        wa_age_err_scatter = wa_age_err * sqrt(mswd) * t_student(0.05, number_of_good_grains - 1)

        self.__min_age = min_age
        self.__max_age = max_age


        self.__min_206_238 = min_206_238
        self.__max_206_238 = max_206_238

        self.__min_207_235 = min_207_235
        self.__max_207_235 = max_207_235

        self.__min_238_206 = min_238_206
        self.__max_238_206 = max_238_206

        self.__min_207_206 = min_207_206
        self.__max_207_206 = max_207_206


        return [number_of_good_grains, wa_age, wa_age_err, wa_age_err_scatter, mswd, max_age, min_age]


    # calculates probability density function for a given age
    def pdp_calc(self, p_age_needed, unc_type):
        if bool(self.good_set):
            sum_gauss = 0
            for key, value in self.good_set.items():
                error = 2 * value[unc_type]
                age = value[0]
                sum_gauss = sum_gauss + (1 / (error * sqrt2pi)) * exp((-(p_age_needed - age) ** 2) / (2 * error ** 2))
            return 1 / len(self.good_set) * sum_gauss

    # calculates kernel density estimate for a given age
    def kde_calc(self, p_age_needed, p_bandwidth):
        if bool(self.good_set):
            sum_gauss = 0
            for key, value in self.good_set.items():
                age = value[0]
                sum_gauss = sum_gauss + (1 / (p_bandwidth * sqrt2pi)) * exp(
                    (-(p_age_needed - age) ** 2) / (2 * (p_bandwidth) ** 2))
            return 1 / (p_bandwidth * len(self.good_set)) * sum_gauss

    # fills and returns a list of pdp's for ages from 0 to EarthAge
    def pdp(self, unc_type):
        index = 0
        list_pdp = []
        list_peaks = []
        if bool(self.good_set):
            while index < EarthAge:
                list_pdp.append(self.pdp_calc(index, unc_type))
                if index > 1 and list_pdp[index - 2] < list_pdp[index - 1] and list_pdp[index] < list_pdp[index - 1]:  # peak recognizing
                    list_peaks.append(index - 1)
                index += 1
            return [list_pdp, list_peaks]

    # fills and returns a list of kde's for ages from 0 to EarthAge
    def kde(self, p_bandwidth):
        index = 0
        # stores non-normalized kde values:
        temp_list = []
        list_peaks = []
        # stores cumulate kde
        ckde = 0
        if bool(self.good_set):
            while index < EarthAge:
                curr_kde = self.kde_calc(index, p_bandwidth)
                temp_list.append(curr_kde)
                ckde += curr_kde
                if index > 1 and temp_list[index - 2] < temp_list[index - 1] and temp_list[index] < temp_list[
                    index - 1]:  # peak recognizing
                    list_peaks.append(index - 1)
                index += 1
            list_kde = [i * (1 / ckde) for i in temp_list]
            return [list_kde, list_peaks]

    # fills and returns a list of cumulative pdp's for ages from 0 to EarthAge
    def cpdp(self, unc_type):
        if bool(self.good_set):
            pdp = self.pdp(unc_type)[0]
            list_pdp = []
            list_pdp.append(pdp[0])
            for index in range(1, len(pdp)):
                list_pdp.append(list_pdp[index - 1] + pdp[index])
            return list_pdp

    # fills and returns a list of cumulative kde's for ages from 0 to EarthAge
    def ckde(self, p_bandwidth):
        if bool(self.good_set):
            kde = self.kde(p_bandwidth)[0]
            list_ckde = []
            list_ckde.append(kde[0])
            for index in range(1, len(kde)):
                list_ckde.append(list_ckde[index - 1] + kde[index])
            return list_ckde


# calculates KS d-value
def d_value(list1, list2):
    d = 0
    if list1 is not None and list2 is not None:
        for age in range(len(list1)):
            i = abs(list1[age] - list2[age])
            if i > d:
                d = i
    else:
        d = -1
    return d


# calculates KS p-value
def p_value(d_val, n1, n2):
    if d_val != -1:
        ne = n1 * n2 / (n1 + n2)
        j = 1
        probks = 0
        sum_j = 0
        lam = sqrt(ne) * d_val
        while j < 100:
            sum_j += ((-1) ** (j - 1)) * exp(-2 * (j ** 2) * (lam ** 2))
            j += 1
        probks = sum_j * 2
        if probks > 1:
            probks = 1
    else:
        probks = -1
    return probks


# goes through the analyses names in AnalysesSet, returns list of samples
def same_sample_set(p_set):
    prev_str = ""
    lset = []
    list_of_analyses_set = []
    i = 0
    p_set.analyses_list.sort(key=lambda x: x.analysis_name, reverse=False)
    while i < len(p_set):
        an = p_set.analyses_list[i]
        temp_str = parse_sample_analysis(str(an))[0]#str(an).rpartition(p_str)[0]
        if (temp_str != prev_str) and (prev_str == ""):
            lset.append(an)
            prev_str = temp_str
        elif (temp_str != prev_str) and (prev_str == ""):
            lset.append(an)
            list_of_analyses_set.append(l_an_set)
        elif (temp_str != prev_str) and (prev_str != ""):
            l_an_set = AnalysesSet(list(reversed(lset)), prev_str)
            list_of_analyses_set.append(l_an_set)
            del l_an_set
            lset = []
            lset.append(an)
            prev_str = temp_str
            if i + 1 == len(p_set): #for a case with the last sample with a single grain
                list_of_analyses_set.append(AnalysesSet(lset, prev_str))
        elif (temp_str == prev_str) and (i + 1 < len(p_set)):
            lset.append(an)
        elif (temp_str == prev_str) and (i + 1 == len(p_set)):
            lset.append(an)
            l_an_set = AnalysesSet(list(reversed(lset)), prev_str)
            list_of_analyses_set.append(l_an_set)
        i += 1
    return list_of_analyses_set


def conf_lim(sigma_level):
    if sigma_level == 1:
        return 0.6826
    elif sigma_level == 2:
        return 0.95
    else:
        return -1

def parse_sample_analysis(full_name):
    last_underscore = full_name.rfind('_')
    last_dash = full_name.rfind('-')
    last_comma = full_name.rfind(',')
    last_dot = full_name.rfind('.')
    pos = max(last_comma, last_dash, last_dot, last_underscore)
    sample_number = full_name[:pos]
    analysis_number = full_name[pos+1:]
    return sample_number, analysis_number

def calc_peaks_weight(peaks, an_set):
    peak_weight = dict.fromkeys(peaks, 0)
    for zircon, zircon_age in an_set.good_set.items():
        for peak in peak_weight:
            if abs(zircon_age[0] - peak) < (zircon_age[1] * 2):  # if difference between a peak and a given age is <2σ
                peak_weight[peak] = peak_weight[peak] + 1
            else:
                pass
    weight_sum = sum(peak_weight.values())
    for peak in peak_weight:
        peak_weight[peak] = round(peak_weight[peak] / weight_sum, 2)
    return peak_weight
