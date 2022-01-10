from math_module import *
from math import sqrt
try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle



def save_object(obj, filename):
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

def load_object(filename):
    file = open(filename, 'rb')
    object_file = pickle.load(file)
    file.close()
    return object_file

#this routine imports a file, checks whether it was originated in Iolite or Glitter and returns that value;
#deletes empty lines if present, returns non-empty lines and their number.
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

        elif any("Dezirteer_template" in s for s in lines):
            file_type = "template"
            length = len(lines)

        else:
            for str in range(len(lines)):
                lines[str] = lines[str].replace(" ", "")
                lines[str] = lines[str].replace("\t\t", "\tNAN\t")
            file_type = "iolite"
            length = len(lines)
    return [lines, file_type, length]


# this routine
def header_pos(imported_list):
    if imported_list[1] == "iolite":
        l_list = []
        file_header = imported_list[0][0].split()
        file_header.pop(2)
        l_list.append(file_header.index('Sourcefile'))
        l_list.append(file_header.index('Duration(s)'))

        if 'Final206_238' in file_header:
            if 'Final206_238_Prop2SE' in file_header:
                prop = file_header.index('Final206_238_Prop2SE')
            else:
                prop = -1
            l_list.append([file_header.index('Final206_238'), file_header.index('Final206_238_Int2SE'), prop])
        else:
            l_list.append([-1, -1, -1])

        if 'Final207_235' in file_header:
            if 'Final207_235_Prop2SE' in file_header:
                prop = file_header.index('Final207_235_Prop2SE')
            else:
                prop = -1
            l_list.append([file_header.index('Final207_235'), file_header.index('Final207_235_Int2SE'), prop])
        else:
            l_list.append([-1, -1, -1])

        if 'ErrorCorrelation_6_38vs7_35' in file_header:
            l_list.append(file_header.index('ErrorCorrelation_6_38vs7_35'))
        else:
            l_list.append(-1)

        if 'Final208_232' in file_header:
            if 'Final208_232_Prop2SE' in file_header:
                prop = file_header.index('Final208_232_Prop2SE')
            else:
                prop = -1
            l_list.append([file_header.index('Final208_232'), file_header.index('Final208_232_Int2SE'), prop])
        else:
            l_list.append([-1, -1, -1])

        if 'Final207_206' in file_header:
            if 'Final207_206_Prop2SE' in file_header:
                prop = file_header.index('Final207_206_Prop2SE')
            else:
                prop = -1
            l_list.append([file_header.index('Final207_206'), file_header.index('Final207_206_Int2SE'), prop])
        else:
            l_list.append([-1, -1, -1])

        if 'Approx_U_PPM' in file_header:
            if 'Approx_U_PPM_Prop2SE' in file_header:
                prop = file_header.index('Approx_U_PPM_Prop2SE')
            else:
                prop = -1
            l_list.append([file_header.index('Approx_U_PPM'), file_header.index('Approx_U_PPM_Int2SE'), prop])
        else:
            l_list.append([-1, -1, -1])

        if 'Pb204' in file_header:
            if 'Pb204_Prop2SE' in file_header:
                prop = file_header.index('Pb204_Prop2SE')
            else:
                prop = -1
            l_list.append([file_header.index('Pb204'), file_header.index('Pb204_Int2SE'), prop])
        else:
            l_list.append([-1, -1, -1])

        if 'ErrorCorrelation_38_6vs7_6' in file_header:
            l_list.append(file_header.index('ErrorCorrelation_38_6vs7_6'))
        else:
            l_list.append(-1)

        if 'Final206_204' in file_header:
            if 'Final206_204_Prop2SE' in file_header:
                prop = file_header.index('Final206_204_Prop2SE')
            else:
                prop = -1
                l_list.append([file_header.index('Final206_204'), file_header.index('Final206_204_Int2SE'), prop])
        else:
            l_list.append([-1, -1, -1])

        if 'Final207_204' in file_header:
            if 'Final207_204_Prop2SE' in file_header:
                prop = file_header.index('Final207_204_Prop2SE')
            else:
                prop = -1
            l_list.append([file_header.index('Final207_204'), file_header.index('Final207_204_Int2SE'), prop])
        else:
            l_list.append([-1, -1, -1])

        if 'Final208_204' in file_header:
            if 'Final208_204_Prop2SE' in file_header:
                prop = file_header.index('Final208_204_Prop2SE')
            else:
                prop = -1
            l_list.append([file_header.index('Final208_204'), file_header.index('Final208_204_Int2SE'), -1])
        else:
            l_list.append([-1, -1, -1])

        if 'Final232_204' in file_header:
            if 'Final232_204_Prop2SE' in file_header:
                prop = file_header.index('Final232_204_Prop2SE')
            else:
                prop = -1
            l_list.append([file_header.index('Final232_204'), file_header.index('Final232_204_Int2SE'), prop])
        else:
            l_list.append([-1, -1, -1])

        if 'Final238_204' in file_header:
            if 'Final238_204_Prop2SE' in file_header:
                prop = file_header.index('Final238_204_Prop2SE')
            else:
                prop = -1
            l_list.append([file_header.index('Final238_204'), file_header.index('Final238_204_Int2SE'), prop])
        else:
            l_list.append([-1, -1, -1])

        if 'Final_U_Th_Ratio' in file_header:
            if 'Final_U_Th_Ratio_Prop2SE' in file_header:
                prop = file_header.index('Final_U_Th_Ratio_Prop2SE')
            else:
                prop = -1
            l_list.append([file_header.index('Final_U_Th_Ratio'), file_header.index('Final_U_Th_Ratio_Int2SE'), prop])
        else:
            l_list.append([-1, -1, -1])

    elif imported_list[1] == "glitter":
        pass
    return l_list


def find_in_glitter(lst, predicate):
    return next((i for i, j in enumerate(lst) if predicate(j)), -1)


#returns an object of Analysis class from a given row in file
def file_to_analysis(imp_file, index):
    full_data = imp_file[0]
    pb206_u238 = []
    pb207_u235 = []
    pb208_th232 = []
    pb207_pb206 = []
    u_conc = []
    pbc = []
    pb206_pb204 = []
    pb207_pb204 = []
    pb208_pb204 = []
    th232_u238 = []
    th232_pb204 = []
    u238_pb204 = []
    final_U_Th_Ratio = []

    if imp_file[1] == 'iolite':  # iolite routine
        #replacing Iolite NaNs and no values
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

        pb206_u238.append(float(an[header[2][0]]))
        pb206_u238.append(float(an[header[2][1]]) / sigma_level)
        if header[2][2] != -1:
            pb206_u238.append(float(an[header[2][2]]) / sigma_level)
        else:
            pb206_u238.append(float(an[header[2][1]]) / sigma_level)

        pb207_u235.append(float(an[header[3][0]]))
        pb207_u235.append(float(an[header[3][1]]) / sigma_level)
        if header[3][2] != -1:
            pb207_u235.append(float(an[header[3][2]]) / sigma_level)
        else:
            pb207_u235.append(float(an[header[3][1]]) / sigma_level)

        pb208_th232.append(float(an[header[5][0]]))
        pb208_th232.append(float(an[header[5][1]]) / sigma_level)
        if header[5][2] != -1:
            pb208_th232.append(float(an[header[5][2]]) / sigma_level)
        else:
            pb208_th232.append(float(an[header[5][1]]) / sigma_level)

        pb207_pb206.append(float(an[header[6][0]]))
        pb207_pb206.append(float(an[header[6][1]]) / sigma_level)
        if header[6][2] != -1:
            pb207_pb206.append(float(an[header[6][2]]) / sigma_level)
        else:
            pb207_pb206.append(float(an[header[6][1]]) / sigma_level)

        rho = calc_rho(pb206_u238[0], pb206_u238[1], pb207_u235[0], pb207_u235[1], pb207_pb206[0], pb207_pb206[1])
        corr_coef_75_68_calculated = rho[0]
        corr_coef_86_76_calculated = rho[1]

        if header[4] != -1:
            corr_coef_75_68 = float(an[header[4]])
        else:
            corr_coef_75_68 = corr_coef_75_68_calculated

        if header[9] != -1:
            corr_coef_86_76 = float(an[header[9]])
        else:
            corr_coef_86_76 = corr_coef_86_76_calculated

        if header[7][0] != -1:
            u_conc.append(float(an[header[7][0]]))
            u_conc.append(float(an[header[7][1]]) / sigma_level)
            if header[7][2] != -1:
                u_conc.append(float(an[header[7][2]]) / sigma_level)
            else:
                u_conc.append(float(an[header[7][1]]) / sigma_level)
        else:
            u_conc.append(-1)
            u_conc.append(-1)
            u_conc.append(-1)

        if header[8][0] != -1:
            pbc.append(float(an[header[8][0]]))
            pbc.append(float(an[header[8][1]]) / sigma_level)
            if header[8][2] != -1:
                pbc.append(float(an[header[8][2]]) / sigma_level)
            else:
                pbc.append(float(an[header[8][1]]) / sigma_level)
        else:
            pbc.append(-1)
            pbc.append(-1)
            pbc.append(-1)

        if header[10][0] != -1:
            pb206_pb204.append(float(an[header[10][0]]))
            pb206_pb204.append(float(an[header[10][1]]) / sigma_level)
            if header[10][2] != -1:
                pb206_pb204.append(float(an[header[10][2]]) / sigma_level)
            else:
                pb206_pb204.append(float(an[header[10][1]]) / sigma_level)
        else:
            pb206_pb204.append(-1)
            pb206_pb204.append(-1)
            pb206_pb204.append(-1)

        if header[11][0] != -1:
            pb207_pb204.append(float(an[header[11][0]]))
            pb207_pb204.append(float(an[header[11][1]]) / sigma_level)
            if header[11][2] != -1:
                pb207_pb204.append(float(an[header[11][2]]) / sigma_level)
            else:
                pb207_pb204.append(float(an[header[11][1]]) / sigma_level)
        else:
            pb207_pb204.append(-1)
            pb207_pb204.append(-1)
            pb207_pb204.append(-1)

        if header[12][0] != -1:
            pb208_pb204.append(float(an[header[12][0]]))
            pb208_pb204.append(float(an[header[12][1]]) / sigma_level)
            if header[12][2] != -1:
                pb208_pb204.append(float(an[header[12][2]]) / sigma_level)
            else:
                pb208_pb204.append(float(an[header[12][1]]) / sigma_level)
        else:
            pb208_pb204.append(-1)
            pb208_pb204.append(-1)
            pb208_pb204.append(-1)

        if header[13][0] != -1:
            th232_pb204.append(float(an[header[13][0]]))
            th232_pb204.append(float(an[header[13][1]]) / sigma_level)
            if header[13][2] != -1:
                th232_pb204.append(float(an[header[13][2]]) / sigma_level)
            else:
                th232_pb204.append(float(an[header[13][1]]) / sigma_level)
        else:
            th232_pb204.append(-1)
            th232_pb204.append(-1)
            th232_pb204.append(-1)

        if header[14][0] != -1:
            u238_pb204.append(float(an[header[14][0]]))
            u238_pb204.append(float(an[header[14][1]]) / sigma_level)
            if header[14][2] != -1:
                u238_pb204.append(float(an[header[14][2]]) / sigma_level)
            else:
                u238_pb204.append(float(an[header[14][1]]) / sigma_level)
        else:
            u238_pb204.append(-1)
            u238_pb204.append(-1)
            u238_pb204.append(-1)

        if header[15][0] != -1:
            final_U_Th_Ratio.append(float(an[header[15][0]]))
            final_U_Th_Ratio.append(float(an[header[15][1]]) / sigma_level)
            if header[15][2] != -1:
                final_U_Th_Ratio.append(float(an[header[15][2]]) / sigma_level)
            else:
                final_U_Th_Ratio.append(float(an[header[15][1]]) / sigma_level)
        else:
            final_U_Th_Ratio.append(-1)
            final_U_Th_Ratio.append(-1)
            final_U_Th_Ratio.append(-1)

        #TEMP WORKAROUND, SINCE IOLITE 2.X does not export 232/238
        th232_u238.append(-1)
        th232_u238.append(-1 / sigma_level)
        th232_u238.append(-1 / sigma_level)

    elif imp_file[1] == 'glitter':  # glitter routine
        file_len = imp_file[2]
        sigma_level = 1.0  # check if all glitter files have this by default
        pos_isotopic_ratios_line = find_in_glitter(full_data, lambda x: '_Isotopic_ratios.' in x)
        pos_errors_line = find_in_glitter(full_data, lambda x: '_Isotopic_ratios:_1_sigma_uncertainty' in x)
        an = full_data[index + 1].split()
        an_err = full_data[index + 2 + file_len].split()

        analysis_name = an[0]
        exposure_time = 0

        u_conc = [-1, -1, -1]
        pbc = [-1, -1, -1]

        pb207_pb206.append(float(an[1]))
        pb207_pb206.append(float(an_err[1]) / sigma_level)
        pb207_pb206.append(float(an_err[1]) / sigma_level)

        pb206_u238.append(float(an[2]))
        pb206_u238.append(float(an_err[2]) / sigma_level)
        pb206_u238.append(float(an_err[2]) / sigma_level)

        pb207_u235.append(float(an[3]))
        pb207_u235.append(float(an_err[3]) / sigma_level)
        pb207_u235.append(float(an_err[3]) / sigma_level)

        pb208_th232.append(float(an[4]))
        pb208_th232.append(float(an_err[4]) / sigma_level)
        pb208_th232.append(float(an_err[4]) / sigma_level)

        if len(an) >= 6:
            th232_u238.append(float(an[5]))
            th232_u238.append(float(an_err[5]) / sigma_level)
            th232_u238.append(float(an_err[5]) / sigma_level)
        else:
            th232_u238 = [-1, -1, -1]
        pb206_pb204 = [-1, -1, -1]
        pb207_pb204 = [-1, -1, -1]
        pb208_pb204 = [-1, -1, -1]
        th232_pb204 = [-1, -1, -1]
        u238_pb204 = [-1, -1, -1]
        final_U_Th_Ratio = [-1, -1, -1]

        rho = calc_rho(pb206_u238[0], pb206_u238[1], pb207_u235[0], pb207_u235[1], pb207_pb206[0], pb207_pb206[1])
        corr_coef_75_68 = rho[0]
        corr_coef_86_76 = rho[1]


    else: #template
        sigma_level = 1
        an = full_data[index].split(",")
        analysis_name = an[0]
        exposure_time = 0

        pb208_th232.append(float(an[1]))
        pb208_th232.append(float(an[2]) / sigma_level)
        pb208_th232.append(float(an[2]) / sigma_level)

        pb206_u238.append(float(an[3]))
        pb206_u238.append(float(an[4]) / sigma_level)
        pb206_u238.append(float(an[4]) / sigma_level)

        if float(an[5]) != -1:
            pb207_u235.append(float(an[5]))
            pb207_u235.append(float(an[6]) / sigma_level)
            pb207_u235.append(float(an[6]) / sigma_level)
        else:
            temp75=float(an[7])*float(an[3])*U238_U235
            temp75err=temp75*sqrt(((float(an[8]))**2/(float(an[7]))**2)+((float(an[4]))**2/(float(an[3]))**2)+((ERR_U238_U235)**2/(U238_U235)**2))/ sigma_level
            pb207_u235.append(temp75)
            pb207_u235.append(temp75err)
            pb207_u235.append(temp75err)

        pb207_pb206.append(float(an[7]))
        pb207_pb206.append(float(an[8]) / sigma_level)
        pb207_pb206.append(float(an[8]) / sigma_level)

        rho = calc_rho(pb206_u238[0], pb206_u238[1], pb207_u235[0], pb207_u235[1], pb207_pb206[0], pb207_pb206[1])
        corr_coef_75_68 = rho[0]
        corr_coef_86_76 = rho[1]

        u_conc.append(float(an[9]))
        u_conc.append(float(an[10]) / sigma_level)
        u_conc.append(float(an[10]) / sigma_level)

        pbc.append(float(an[11]))
        pbc.append(float(an[12]) / sigma_level)
        pbc.append(float(an[12]) / sigma_level)

        pb206_pb204.append(float(an[13]))
        pb206_pb204.append(float(an[14]) / sigma_level)
        pb206_pb204.append(float(an[14]) / sigma_level)

        pb207_pb204.append(float(an[15]))
        pb207_pb204.append(float(an[16]) / sigma_level)
        pb207_pb204.append(float(an[16]) / sigma_level)

        pb208_pb204.append(float(an[17]))
        pb208_pb204.append(float(an[18]) / sigma_level)
        pb208_pb204.append(float(an[18]) / sigma_level)

        th232_pb204.append(float(an[19]))
        th232_pb204.append(float(an[20]) / sigma_level)
        th232_pb204.append(float(an[20]) / sigma_level)

        u238_pb204.append(float(an[21]))
        u238_pb204.append(float(an[22]) / sigma_level)
        u238_pb204.append(float(an[22]) / sigma_level)

        final_U_Th_Ratio.append(float(an[23]))
        final_U_Th_Ratio.append(float(an[24]) / sigma_level)
        final_U_Th_Ratio.append(float(an[24]) / sigma_level)

        th232_u238.append(float(an[25]))
        th232_u238.append(float(an[26]) / sigma_level)
        th232_u238.append(float(an[26]) / sigma_level)



    l_analysis = Analysis(analysis_name, exposure_time, pb206_u238, pb207_u235, corr_coef_75_68, corr_coef_86_76,
                          pb208_th232, pb207_pb206, u_conc, pbc, pb206_pb204, pb207_pb204, pb208_pb204,
                          th232_pb204, u238_pb204,sigma_level, final_U_Th_Ratio,  th232_u238)
    return l_analysis