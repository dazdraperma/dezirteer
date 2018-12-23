from bisect import *
from const import *
import operator
import functools
from math import sqrt, exp, log
import scipy.stats

pbpb_table = []
concordia_table = []


def t_student(alpha, gl):
    return scipy.stats.t.ppf(1 - (alpha / 2), gl)


def calc_ratio(age):
    pb207_u235 = exp(LAMBDA_235 * age * 1000000) - 1
    pb206_u238 = exp(LAMBDA_238 * age * 1000000) - 1
    pb207_pb206 = (pb207_u235 / pb206_u238) * (1 / U238_U235)
    u238_pb206 = 1 / pb206_u238
    pb208_th232 = exp(LAMBDA_232 * age * 1000000) - 1
    return [pb206_u238, pb207_u235, pb207_pb206, u238_pb206, pb208_th232]


def compb(age, n):  # Stacey & Cramers 2 stage pb evolution model
    if n == 3:
        return compb(age, 1) / compb(age, 0)  # 7/6c
    elif n != 3 and age <= 3700:
        if n == 0:
            return 11.152 + 9.735 * (calc_ratio(4570)[0] - calc_ratio(age)[0])  # 6/4c
        elif n == 1:
            return 12.998 + 9.735 / U238_U235 * (calc_ratio(4570)[1] - calc_ratio(age)[1])  # 7/4c
        elif n == 2:
            return 31.23 + 36.837 * (calc_ratio(4570)[4] - calc_ratio(age)[4])  # 8/4c
    else:
        if n == 0:
            return 9.307 + 7.1925 * (calc_ratio(4570)[0] - calc_ratio(age)[0])  # 6/4c
        elif n == 1:
            return 10.294 + 7.192 / U238_U235 * (calc_ratio(4570)[1] - calc_ratio(age)[1])  # 7/4c
        elif n == 2:
            return 29.476 + 32.208 * (calc_ratio(4570)[4] - calc_ratio(age)[4])  # 8/4c


def pbc_corr(zir, corr_type, *args):  # returns Pbc-corrected ages
    delta = 1
    corr_age = [-1, -1]
    if corr_type == 0:  # 204
        # parts of common lead in measured
        f6 = compb(zir.calc_age(0), 0) / zir.pb206_pb204[0]
        f7 = compb(zir.calc_age(0), 1) / zir.pb207_pb204[0]
        f8 = compb(zir.calc_age(0), 2) / zir.pb208_pb204[0]

        # 204pb corrected ratios
        r68 = zir.pb206_u238[0] * (1 - f6)
        r75 = zir.pb207_u235[0] * (1 - f7)
        r82 = zir.pb208_th232[0] * (1 - f8)
        r76 = r75 / r68 * 1 / U238_U235

        # 204pb corrected ages
        a68 = log(r68 + 1) / LAMBDA_238
        a75 = log(r75 + 1) / LAMBDA_235
        a82 = log(r82 + 1) / LAMBDA_232
        a76 = find_age(r76)

        # age errors
        tmp64 = (zir.pb206_pb204[1] / zir.pb206_pb204[0]) ** 2
        tmp74 = (zir.pb207_pb204[1] / zir.pb207_pb204[0]) ** 2
        tmp84 = (zir.pb208_pb204[1] / zir.pb208_pb204[0]) ** 2
        r68er = sqrt(
            (sqrt(tmp64) / zir.pb206_pb204[0] / (1 - f6)) ** 2 + (zir.pb206_u238[1] / zir.pb206_u238[0]) ** 2) * r68
        r75er = sqrt(
            (sqrt(tmp74) / zir.pb207_pb204[0] / (1 - f7)) ** 2 + (zir.pb207_u235[1] / zir.pb207_u235[0]) ** 2) * r75
        r82er = sqrt(
            (sqrt(tmp84) / zir.pb208_pb204[0] / (1 - f8)) ** 2 + (zir.pb208_th232[1] / zir.pb208_th232[0]) ** 2) * r82
        a68er = r68er / (1 + r68) / LAMBDA_238 / 1000000
        a75er = r75er / (1 + r75) / LAMBDA_235 / 1000000
        a82er = r82er / (1 + r82) / LAMBDA_232 / 1000000
        a76er = zir.pb207_pb206[1]

    elif corr_type == 1:  # 207
        t = 1000
        # age
        while delta > 0.001:
            f = U238_U235 * compb(zir.calc_age(0), 3) * (zir.pb206_u238[0] - calc_ratio(t)[0]) - zir.pb207_u235[0] + \
                calc_ratio(t)[1]
            deriv = LAMBDA_235 * (calc_ratio(t)[1] + 1) - U238_U235 * compb(zir.calc_age(0), 3) * LAMBDA_238 * (
                        calc_ratio(t)[0] + 1)
            delta = -f / deriv
            t = t + delta
        corr_age[0] = t
        # error
        r75var = U238_U235 ** 2 * (
                    (zir.pb206_u238[0] * zir.pb207_pb206[1]) ** 2 + (zir.pb207_pb206[0] * zir.pb206_u238[1]) ** 2)
        denom = (U238_U235 * compb(zir.calc_age(0), 3) * LAMBDA_238 * (calc_ratio(t)[0] + 1) - LAMBDA_235 * (
                    calc_ratio(t)[1] + 1)) ** 2
        n1 = 0  # n1=(U238_U235*(zir.pb206_u238[0]-calc_ratio(t)[0])*) #commonly it's assumed r76c_err=0 => n1=0
        n2 = U238_U235 ** 2 * compb(zir.calc_age(0), 3) * (compb(zir.calc_age(0), 3) - 2 * zir.pb207_pb206[0]) * \
             zir.pb206_u238[1] ** 2
        n = n1 + n2 + r75var
        corr_age[1] = sqrt(n / denom)

    elif corr_type == 2:  # 208
        t = 1000
        # age
        while delta > 0.001:
            f = zir.pb206_u238[0] - (calc_ratio(t)[0] + 1) + 1 - zir.th232_pb204[0] / zir.u238_pb204[0] * compb(
                zir.calc_age(0), 0) / \
                compb(zir.calc_age(0), 2) * (zir.pb208_th232[0] - (calc_ratio(t)[4] + 1) + 1)
            deriv = zir.th232_pb204[0] / zir.u238_pb204[0] * compb(zir.calc_age(0), 0) / compb(zir.calc_age(0),
                                                                                               2) * LAMBDA_232 * \
                    (calc_ratio(t)[4] + 1) - LAMBDA_238 * (calc_ratio(t)[0] + 1)
            delta = -f / deriv
            t = t + delta  #
        corr_age[0] = t
        # error
        c1 = 0
        c2 = zir.th232_pb204[0] / zir.u238_pb204[0] * compb(zir.calc_age(0), 0) / compb(zir.calc_age(0), 2)
        n1 = 0
        n2 = c2 * zir.pb208_th232[1] * zir.pb206_u238[1] ** 2
        n = n1 + n2
        denom = (c2 * LAMBDA_232 * (calc_ratio(t)[4] + 1) - LAMBDA_238 * (calc_ratio(t)[0] + 1)) ** 2
        corr_age[1] = sqrt(n / denom)

    elif corr_type == 3:  # and
        # age
        t1 = zir.pb207_u235
        xt2 = calc_ratio(t2)[1]
        yt2 = calc_ratio(t2)[0]
        zt2 = calc_ratio(t2)[4]
        c7=15.628/18.7
        c8=38.63/18.7
        x=zir.pb207_u235[0]
        y=zir.pb206_u238[0]
        z=zir.pb208_th232[0]
        u=zir.u238_pb204/zir.th232_pb204
        k=U238_U235
        corr_age[0] = andersen(t1, xt2, yt2, zt2, c7, c8, x, y, z, u, k)[0]
        mkages=[]
        mkfc=[]
        #sigma errors
        for i in range(100):
            mkx=random.normalvariate(0,1) 
            mky=mkx*zir.corr_coef_75_68+sqrt(1-zir.corr_coef_75_68**2)*random.normalvariate(0,1)
            mkz=random.normalvariate(0,1)
            mkx=mkx*zir.pb207_u235[1]+x
            mky=mky*zir.pb206_u238[1]+y
            mkz=mkz*zir.pb208_th232[1]+z
            t1=calc_age(mkx)[0]
            mkages.append(andersen(t1,xt2,yt2,zt2,c7,c8,mkx,mky,mkz,u)[0])
            mkfc.append(andersen(t1,xt2,yt2,zt2,c7,c8,mkx,mky,mkz,u)[1])
        ageer=np.std(mkages)
        fcer=np.std(mkfc)
        corr_age[1]=ageer

    else:
        corr_age = [-1, -1]
    return corr_age


def zero(t1, xt2, yt2, zt2, c7, c8, x, y, z, u, k):
    xt1 = calc_ratio(t1)[1]
    yt1 = calc_ratio(t1)[0]
    zt1 = calc_ratio(t1)[4]
    eq7 = (y * (xt1 - xt2) - yt2 * xt1 + x * (yt2 - yt1) + xt2 * yt1) / (
            xt1 - xt2 - c7 * k * yt1 + c7 * k * yt2) - (
                  z * (yt2 - yt1) + zt2 * yt1 + y * (zt1 - zt2) - yt2 * zt1) / (
                  zt1 - zt2 - c8 * u * yt1 + c8 * u * yt2)
    return eq7


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
        t = 100
        ratios = []
        while t <= EarthAge:
            ratios.append(exp(LAMBDA_238 * t * 1000000) - 1)
            ratios.append(exp(LAMBDA_235 * t * 1000000) - 1)
            pbpb_age = (1 / U238_U235) * (exp(LAMBDA_235 * t * 1000000) - 1) / (exp(LAMBDA_238 * t * 1000000) - 1)
            ratios.append(pbpb_age)
            concordia_table.append(ratios)
            ratios = []
            t += 1


# calculates pb-pb age from the table of pb-pb ratios,filled in fill_Pb206207_table()
def find_age(pLeadRatio):
    return bisect(pbpb_table, pLeadRatio) * 5


class Filters(object):  # describes filters that should be applied to data in Analysis_set object
    def __init__(self, filter_by_uconc=[False, 1000], which_age=[1, 1000], use_pbc=False,
                 filter_by_err=[False, 0.1], include207235Err=False,
                 pos_disc_filter=0.2, neg_disc_filter=-0.1, disc_type=1,
                 sample_name_filter=[], unc_type='0'):
        self.__filter_by_uconc = filter_by_uconc
        self.__which_age = which_age
        self.__use_pbc = use_pbc
        self.__filter_by_err = filter_by_err
        self.__include207235Err = include207235Err
        self.__pos_disc_filter = pos_disc_filter
        self.__neg_disc_filter = neg_disc_filter
        self.__disc_type = disc_type
        self.__sample_name_filter = sample_name_filter
        self.__unc_type = unc_type

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


def imported_file(p_file_name):
    file_type = ""
    length = 0
    with open(p_file_name) as f_in:
        lines = list(line for line in (l.strip() for l in f_in) if line)  # deleting empty lines

        if any("GLITTER" in s for s in lines):
            file_type = "glitter"
            lines.pop(2)
            lines.pop(1)
            lines.pop(0)
            pos_isotopic_ratios_line = find_in_glitter(lines, lambda x: '_Isotopic_ratios.' in x)
            pos_errors_line = find_in_glitter(lines, lambda x: '_Isotopic_ratios:_1_sigma_uncertainty' in x)
            temp_list = lines[0: pos_errors_line * 2]
            lines = temp_list
            length = int(len(lines) / 2) - 1
        else:
            for str in range(len(lines)):
                lines[str] = lines[str].replace(" ", "")
                lines[str] = lines[str].replace("\t\t", "\tNAN\t")
            file_type = "iolite"
            length = len(lines)
    return [lines, file_type, length]


def header_pos(imported_list):
    if imported_list[1] == "iolite":
        l_list = []
        file_header = imported_list[0][0].split()
        file_header.pop(2)
        l_list.append(file_header.index('Sourcefile'))
        l_list.append(file_header.index('Duration(s)'))

        if 'Final206_238' in file_header:
            l_list.append([file_header.index('Final206_238'), file_header.index('Final206_238_Int2SE')])
        else:
            l_list.append([-1, -1])

        if 'Final207_235' in file_header:
            l_list.append([file_header.index('Final207_235'), file_header.index('Final207_235_Int2SE')])
        else:
            l_list.append([-1, -1])

        if 'ErrorCorrelation_6_38vs7_35' in file_header:
            l_list.append(file_header.index('ErrorCorrelation_6_38vs7_35'))
        else:
            l_list.append(0.99)

        if 'Final208_232' in file_header:
            l_list.append([file_header.index('Final208_232'), file_header.index('Final208_232_Int2SE')])
        else:
            l_list.append([-1, -1])

        if 'Final207_206' in file_header:
            l_list.append([file_header.index('Final207_206'), file_header.index('Final207_206_Int2SE')])
        else:
            l_list.append([-1, -1])

        if 'Approx_U_PPM' in file_header:
            l_list.append([file_header.index('Approx_U_PPM'), file_header.index('Approx_U_PPM_Int2SE')])
        else:
            l_list.append([-1, -1])

        if 'Pb204' in file_header:
            l_list.append([file_header.index('Pb204'), file_header.index('Pb204_Int2SE')])
        else:
            l_list.append([-1, -1])

        if 'ErrorCorrelation_38_6vs7_6' in file_header:
            l_list.append(file_header.index('ErrorCorrelation_38_6vs7_6'))
        else:
            l_list.append(0.99)

        if 'Final206_204' in file_header:
            l_list.append([file_header.index('Final206_204'), file_header.index('Final206_204_Int2SE')])
        else:
            l_list.append([-1, -1])

        if 'Final207_204' in file_header:
            l_list.append([file_header.index('Final207_204'), file_header.index('Final207_204_Int2SE')])
        else:
            l_list.append([-1, -1])

        if 'Final208_204' in file_header:
            l_list.append([file_header.index('Final208_204'), file_header.index('Final208_204_Int2SE')])
        else:
            l_list.append([-1, -1])

        if 'Final232_204' in file_header:
            l_list.append([file_header.index('Final232_204'), file_header.index('Final232_204_Int2SE')])
        else:
            l_list.append([-1, -1])

        if 'Final238_204' in file_header:
            l_list.append([file_header.index('Final238_204'), file_header.index('Final238_204_Int2SE')])
        else:
            l_list.append([-1, -1])

    elif imported_list[1] == "glitter":
        pass
    return l_list


def find_in_glitter(lst, predicate):
    return next((i for i, j in enumerate(lst) if predicate(j)), -1)


def file_to_analysis(imp_file, index):
    full_data = imp_file[0]

    if imp_file[1] == 'iolite':  # iolite routine
        str_temp = full_data[index].replace("no value", "-1")
        str_temp = str_temp.replace("novalue", "-1")
        str_temp = str_temp.replace("NAN", "-1")
        an = str_temp.split()  # necessary analysis, split
        an.pop(4)
        an.pop(3)
        header = header_pos(imp_file)
        analysis_name = an[header[0]]
        exposure_time = float(an[header[1]])
        sigma_level = 2

        pb206_u238 = []
        pb206_u238.append(float(an[header[2][0]]))
        pb206_u238.append(float(an[header[2][1]]) / sigma_level)

        pb207_u235 = []
        pb207_u235.append(float(an[header[3][0]]))
        pb207_u235.append(float(an[header[3][1]]) / sigma_level)

        if header[4] != -1:
            corr_coef_75_68 = float(an[header[4]])
        else:
            corr_coef_75_68 = 0.99

        if header[9] != -1:
            corr_coef_86_76 = float(an[header[9]])
        else:
            corr_coef_86_76 = 0.99

        pb208_th232 = []
        pb208_th232.append(float(an[header[5][0]]))
        pb208_th232.append(float(an[header[5][1]]) / sigma_level)

        pb207_pb206 = []
        pb207_pb206.append(float(an[header[6][0]]))
        pb207_pb206.append(float(an[header[6][1]]) / sigma_level)

        u_conc = []
        if header[7][0] != -1:
            u_conc.append(float(an[header[7][0]]))
            u_conc.append(float(an[header[7][1]]) / sigma_level)
        else:
            u_conc.append(-1)
            u_conc.append(-1)

        pbc = []
        if header[8][0] != -1:
            pbc.append(float(an[header[8][0]]))
            pbc.append(float(an[header[8][1]]) / sigma_level)
        else:
            pbc.append(-1)
            pbc.append(-1)

        pb206_pb204 = []
        if header[10][0] != -1:
            pb206_pb204.append(float(an[header[10][0]]))
            pb206_pb204.append(float(an[header[10][1]]) / sigma_level)
        else:
            pb206_pb204.append(-1)
            pb206_pb204.append(-1)

        pb207_pb204 = []
        if header[11][0] != -1:
            pb207_pb204.append(float(an[header[11][0]]))
            pb207_pb204.append(float(an[header[11][1]]) / sigma_level)
        else:
            pb207_pb204.append(-1)
            pb207_pb204.append(-1)

        pb208_pb204 = []
        if header[12][0] != -1:
            pb208_pb204.append(float(an[header[12][0]]))
            pb208_pb204.append(float(an[header[12][1]]) / sigma_level)
        else:
            pb208_pb204.append(-1)
            pb208_pb204.append(-1)

        th232_pb204 = []
        if header[13][0] != -1:
            th232_pb204.append(float(an[header[13][0]]))
            th232_pb204.append(float(an[header[13][1]]) / sigma_level)
        else:
            th232_pb204.append(-1)
            th232_pb204.append(-1)

        u238_pb204 = []
        if header[14][0] != -1:
            u238_pb204.append(float(an[header[14][0]]))
            u238_pb204.append(float(an[header[14][1]]) / sigma_level)
        else:
            u238_pb204.append(-1)
            u238_pb204.append(-1)

    else:  # glitter routine
        file_len = imp_file[2]
        sigma_level = 1.0  # check if all glitter files have this by default
        pos_isotopic_ratios_line = find_in_glitter(full_data, lambda x: '_Isotopic_ratios.' in x)
        pos_errors_line = find_in_glitter(full_data, lambda x: '_Isotopic_ratios:_1_sigma_uncertainty' in x)
        an = full_data[index + 1].split()
        an_err = full_data[index + 2 + file_len].split()

        analysis_name = an[0]
        exposure_time = 0

        u_conc = [-1, -1]
        pbc = [-1, -1]

        pb207_pb206 = []
        pb207_pb206.append(float(an[1]))
        pb207_pb206.append(float(an_err[1]) / sigma_level)

        pb206_u238 = []
        pb206_u238.append(float(an[2]))
        pb206_u238.append(float(an_err[2]) / sigma_level)

        pb207_u235 = []
        pb207_u235.append(float(an[3]))
        pb207_u235.append(float(an_err[3]) / sigma_level)

        pb208_th232 = []
        pb208_th232.append(float(an[4]))
        pb208_th232.append(float(an_err[4]) / sigma_level)

        pb206_pb204 = [-1, -1]
        pb207_pb204 = [-1, -1]
        pb208_pb204 = [-1, -1]
        th232_pb204 = [-1, -1]
        u238_pb204 = [-1, -1]

        corr_coef_75_68 = (pb206_u238[1] / pb206_u238[0]) / (pb207_u235[1] / pb207_u235[0])
        corr_coef_86_76 = ((1 / pb206_u238[1]) / (1 / pb206_u238[0])) / (pb207_pb206[1] / pb207_pb206[0])

    l_analysis = Analysis(analysis_name, exposure_time, pb206_u238, pb207_u235, corr_coef_75_68, corr_coef_86_76,
                          pb208_th232, pb207_pb206, u_conc, pbc, pb206_pb204, pb207_pb204, pb208_pb204,
                          th232_pb204, u238_pb204, sigma_level)
    return l_analysis


class Analysis(object):
    def __init__(self, analysis_name="", exposure_time="",
                 pb206_u238=(0, 0), pb207_u235=(0, 0), corr_coef_75_68=0, corr_coef_86_76=0, pb208_th232=(0, 0),
                 pb207_pb206=(0, 0), u_conc=(0, 0), pbc=(0, 0), pb206_pb204=(0, 0), pb207_pb204=(0, 0),
                 pb208_pb204=(0, 0), th232_pb204=(0, 0), u238_pb204=(0, 0), sigma_level=0):
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

    def u238_pb206(self):
        rat238206 = 1 / self.pb206_u238[0]
        return [rat238206, rat238206 * (self.pb206_u238[1] / self.pb206_u238[0])]

    # calculates age ± error from isotopic value and uncertainty
    def calc_age(self, isotopic_system):
        age_err = -1
        try:
            if isotopic_system == 0 and self.pb206_u238[0] > 0:
                age = (1 / lambdas[isotopic_system]) * log(self.pb206_u238[0] + 1) / 1000000
                age_err = (1 / lambdas[isotopic_system]) * self.pb206_u238[1] / 1000000

            elif isotopic_system == 1 and self.pb207_u235[0] > 0:
                age = (1 / lambdas[isotopic_system]) * log(self.pb207_u235[0] + 1) / 1000000
                age_err = (1 / lambdas[isotopic_system]) * self.pb207_u235[1] / 1000000

            elif isotopic_system == 2 and self.pb208_th232[0] > 0:
                age = (1 / lambdas[isotopic_system]) * log(self.pb208_th232[0] + 1) / 1000000
                age_err = (1 / lambdas[isotopic_system]) * self.pb208_th232[1] / 1000000

            elif isotopic_system == 3 and self.pb207_pb206[0] > .04605:
                # .04605 corresponds to age67 ~ 0
                age = find_age(self.pb207_pb206[0]) * 1000000
                C1 = 1 / U238_U235
                C2 = LAMBDA_235
                C3 = LAMBDA_238
                df = self.pb207_pb206[1]
                dfdt = C1 * (C3 * exp(C3 * age) * (exp(C2 * age) - 1) - C2 * exp(C2 * age) *
                             (exp(C3 * age) - 1)) / ((exp(C3 * age) - 1) ** 2)
                age_err = abs(df / dfdt / 1000000)
                age = age / 1000000
            else:
                age = -1
            return [age, age_err]
        except ValueError:
            pass

    # calculates two type of discordance: (1) between 206_238 and 207_235 and (2) between 206_238 and 206_207
    def calc_discordance(self, between_68_57):
        if between_68_57:
            return self.calc_age(1)[0] / self.calc_age(0)[0] - 1

        else:
            return 1 - self.calc_age(0)[0] / self.calc_age(3)[0]

    # true if should use U238/Pb206, false if should use Pb206/Pb207
    def use_206_238(self, age_cutoff):
        return self.calc_age(0)[0] < age_cutoff

    # checks if a grain passes user-defined Filters
    def is_grain_good(self, pFilter: Filters):
        do_uconc = pFilter.filter_by_uconc[0]
        uconc_ppm_cutoff = pFilter.filter_by_uconc[1]
        which_age = pFilter.which_age[0]
        age_fixed_limit = pFilter.which_age[1]
        do_err = pFilter.filter_by_err[0]
        do_207235_err = pFilter.include207235Err
        err_cutoff = pFilter.filter_by_err[1]
        type_disc = pFilter.disc_type
        pos_disc_cutoff = pFilter.pos_disc_filter
        neg_disc_cutoff = pFilter.neg_disc_filter
        sample_name_filter = pFilter.sample_name_filter
        is_grain_in_chosen_sample = True

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
            if str(self.analysis_name).rpartition('_')[0] in sample_name_filter:
                is_grain_in_chosen_sample = True
            else:
                is_grain_in_chosen_sample = False

        # decide on the default age system
        if (which_age == 1 and self.calc_age(0)[0] > age_fixed_limit) or which_age == 2:
            age_68_67 = 1
            this_age = 3
        elif (which_age == 1 and self.calc_age(0)[0] < age_fixed_limit) or which_age == 3:
            age_68_67 = 0
            this_age = 0
        else:
            age_68_67 = -1  # this is for future implementation of 'based on Uconc' algorithm

        # filter by measurement error
        if do_err:
            age = self.calc_age(age_68_67 * 3)[0]  # calculates either 68 or 67 age
            err = self.calc_age(age_68_67 * 3)[1]  # calculates correspondent error
            if err / age < err_cutoff:  # checks whether error is within limit
                is_err_good = True
            else:
                is_err_good = False
            if do_207235_err == 1:
                age207235 = self.calc_age(1)[0]
                age207235err = self.calc_age(1)[1]
                if age207235err / age207235 < err_cutoff:
                    is_207235err_good = True
                else:
                    is_207235err_good = False
            else:
                is_207235err_good = True  # if no need to filter by 207-235 age error
        else:
            is_err_good = True  # if no need to filter by age error
            is_207235err_good = True

        # filter by discordance #0: u-conc, #1 - fixed limit, #2 - 57-86, #3- 67-86
        if type_disc == 0:
            disc = -1  # this is for future implementation of 'based on Uconc' algorithm
        elif type_disc == 1:
            disc = self.calc_discordance(not age_68_67)
        elif type_disc == 2:
            disc = self.calc_discordance(True)
        elif type_disc == 3:
            disc = self.calc_discordance(False)  # 12345

        if (disc < pos_disc_cutoff) & (disc > neg_disc_cutoff):
            is_disc_good = True
        else:
            is_disc_good = False

        return are_ratios_positive & is_uconc_good & is_err_good & is_207235err_good & is_disc_good & \
               is_grain_in_chosen_sample, this_age


class AnalysesSet(object):
    def __init__(self, analyses_list, name):
        self.__analyses_list = analyses_list
        self.__name = name
        self.__good_set = {}
        self.__bad_set = []

    def __repr__(self):
        return str(self.analyses_list)

    def __len__(self):
        return len(self.analyses_list)

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
    # number of grains, weighted average age ± uncertainty (±1s and 95%), MSWD, max and min ages
    def good_bad_sets(self, p_filter: Filters):
        index = 0
        max_age = 0
        min_age = 5000
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
        self.good_set.clear()
        self.bad_set.clear()
        while index <= len(self.analyses_list) - 1:
            zircon = self.analyses_list[index]
            l_is_grain_good = Analysis.is_grain_good(zircon, p_filter)
            if l_is_grain_good[0]:
                z_age = zircon.calc_age(l_is_grain_good[1])
                if z_age[0] > max_age:
                    max_age = int(z_age[0])
                if z_age[0] < min_age:
                    min_age = int(z_age[0])
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
        return [number_of_good_grains, wa_age, wa_age_err, wa_age_err_scatter, mswd, max_age, min_age]

    # calculates probability density function for a given age
    def pdp_calc(self, p_age_needed):
        if bool(self.good_set):
            sum_gauss = 0
            for key, value in self.good_set.items():
                error = 2 * value[1]
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
    def pdp(self):
        index = 0
        list_pdp = []
        list_peaks = []
        if bool(self.good_set):
            while index < EarthAge:
                list_pdp.append(self.pdp_calc(index))
                if index > 1 and list_pdp[index - 2] < list_pdp[index - 1] and list_pdp[index] < list_pdp[
                    index - 1]:  # peak recognizing
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
    def cpdp(self):
        if bool(self.good_set):
            pdp = self.pdp()[0]
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
def same_sample_set(p_set: AnalysesSet, p_str):
    prev_str = ""
    lset = []
    list_of_analyses_set = []
    i = 0
    p_set.analyses_list.sort(key=lambda x: x.analysis_name, reverse=False)
    while i < len(p_set):
        an = p_set.analyses_list[i]
        temp_str = str(an).rpartition(p_str)[0]
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


def calc_peaks_weight(peaks: [], an_set: AnalysesSet):
    peak_weight = dict.fromkeys(peaks, 0)
    for zircon, zircon_age in an_set.good_set.items():
        for peak in peak_weight:
            if abs(zircon_age[0] - peak) < (zircon_age[1] * 2):  # if difference between a peak and a given age is <2σ
                peak_weight[peak] = peak_weight[peak] + 1
            else:
                pass
    return peak_weight