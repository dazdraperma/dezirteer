# -*- coding: utf-8 -*-
from math import *
#import sys
import ctypes
try:
    from Tkinter import *
except ImportError:
    from tkinter import *

try:
    import ttk
    py3 = False
except ImportError:
    import tkinter.ttk as ttk
    py3 = True

from const import *
from math_module import *
global g_corr_type
g_corr_type = ["none", "204", "207", "208", "And."]

class NoGrainsError(Exception):
    pass

#class containing settings to be applied to graph
class GraphSettings(object):
    def __init__(self, show_multiple=False, conc_type=0, ellipses_at=1, fit_disc=False,
                 anchored=False, pdp_kde_hist=0, do_hist=False, bandwidth=50, binwidth=50, keep_previous=False):
        self.__show_multiple = show_multiple
        self.__conc_type = conc_type
        self.__ellipses_at = ellipses_at
        self.__fit_disc = fit_disc
        self.__anchored = anchored
        self.__pdp_kde_hist = pdp_kde_hist
        self.__do_hist = do_hist
        self.__bandwidth = bandwidth
        self.__binwidth = binwidth
        self.__keep_previous = keep_previous

    @property
    def show_multiple(self):
        return self.__show_multiple

    @show_multiple.setter
    def show_multiple(self, value):
        self.__show_multiple = value

    @property
    def conc_type(self):
        return self.__conc_type

    @conc_type.setter
    def conc_type(self, value):
        self.__conc_type = value

    @property
    def ellipses_at(self):
        return self.__ellipses_at

    @ellipses_at.setter
    def ellipses_at(self, value):
        self.__ellipses_at = value

    @property
    def fit_disc(self):
        return self.__fit_disc

    @fit_disc.setter
    def fit_disc(self, value):
        self.__fit_disc = value

    @property
    def anchored(self):
        return self.__anchored

    @anchored.setter
    def anchored(self, value):
        self.__anchored = value

    @property
    def pdp_kde_hist(self):
        return self.__pdp_kde_hist

    @pdp_kde_hist.setter
    def pdp_kde_hist(self, value):
        self.__pdp_kde_hist = value

    @property
    def do_hist(self):
        return self.__do_hist

    @do_hist.setter
    def do_hist(self, value):
        self.__do_hist = value

    @property
    def bandwidth(self):
        return self.__bandwidth

    @bandwidth.setter
    def bandwidth(self, value):
        self.__bandwidth = value

    @property
    def binwidth(self):
        return self.__binwidth

    @binwidth.setter
    def binwidth(self, value):
        self.__binwidth = value

    @property
    def keep_previous(self):
        return self.__keep_previous

    @keep_previous.setter
    def keep_previous(self, value):
        self.__keep_previous = value


def onChange(p_number_in_list, p_value, pars, *args, **kwargs):
    #'''p_filters, p_table, p_grainset, p_colnames''' p_table, p_grainset, p_filters, p_colnames

    if p_number_in_list == 1:
        pars[0].show_multiple = p_value
    elif p_number_in_list == 2:
        pars[0].filter_by_uconc[0] = True if p_value == 1 else False
    elif p_number_in_list == 3:
        pars[0].which_age[0] = p_value
    elif p_number_in_list == 4:
        pars[0].use_pbc = [p_value, varAgeAndersen.get()]
    elif p_number_in_list == 5:
        pars[0].filter_by_err[0] = p_value
    elif p_number_in_list == 6:
        pars[0].pos_disc_filter = p_value/100
    elif p_number_in_list == 7:
        pars[0].neg_disc_filter = p_value/100*(-1)
    elif p_number_in_list == 8:
        pars[0].disc_type[0] = p_value
    elif p_number_in_list == 9:
        pars[0].conc_type = p_value
    elif p_number_in_list == 11:
        pars[0].fit_disc = p_value
    elif p_number_in_list == 12:
        pars[0].anchored = p_value
    elif p_number_in_list == 13:
        pars[0].do_pdp = p_value
    elif p_number_in_list == 14:
        pars[0].do_kde = p_value
    elif p_number_in_list == 15:
        pars[0].do_cpdp = p_value
    elif p_number_in_list == 16:
        pars[0].do_ckde = p_value
    elif p_number_in_list == 17:
        pars[0].do_hist = p_value
    elif p_number_in_list == 18:
        pars[0].filter_by_uconc[1] = p_value
    elif p_number_in_list == 19:
        pars[0].which_age[1] = p_value
    elif p_number_in_list == 20:
        pars[0].filter_by_err[1] = p_value/100
    elif p_number_in_list == 21:
       pass
    elif p_number_in_list == 22:
        pars[0].include207235Err = p_value
    elif p_number_in_list == 23:
        pars[0].unc_type = p_value
    elif p_number_in_list == 24:
        pars[0].filter_by_commPb = p_value
    elif p_number_in_list == 25:
        pars[0].disc_type[1] = p_value
    elif p_number_in_list == 26:
        if varMinAgeCrop.get() == 0:
            pars[0].minAgeCrop = 0
            try:
                pars[0].minAgeCrop = float(p_value)
            except ValueError:
                pars[0].minAgeCrop = 0
    elif p_number_in_list == 27:
        if varMaxAgeCrop.get() == 0:
            pars[0].maxAgeCrop = EarthAge
        else:
            try:
                pars[0].maxAgeCrop = float(p_value)
            except ValueError:
                pars[0].maxAgeCrop = 0
    elif p_number_in_list == 28:
        pars[0].andersenAge = p_value
    elif p_number_in_list == 29:
        pars[0].discOrIntersect = p_value
    elif p_number_in_list == 30:
        pars[0].intersectAt = p_value+1

    sys.stdout.flush()
    fill_data_table(pars[1], pars[2], pars[0], pars[3])
    set_all_ui_elements(args[0])
#'''p_filters, p_table, p_grainset, p_colnames''' p_table, p_grainset, p_filters, p_colnames


def onGraphChange(p_graph_settings, p_number_in_list, p_value, *args, **kwargs):
    if p_number_in_list == 0:
        p_graph_settings.conc_type = p_value
    elif p_number_in_list == 1:
        p_graph_settings.conc_type = p_value
    elif p_number_in_list == 2:
        p_graph_settings.ellipses_at = p_value
    elif p_number_in_list == 7:
        p_graph_settings.pdp_kde_hist = p_value
        if p_value == 0:
            args[0].configure(state="normal")
            args[1].configure(state="disabled")
        elif p_value == 1:
            args[0].configure(state="disabled")
            args[1].configure(state="disabled")
        else:
            args[0].configure(state="disabled")
            args[1].configure(state="normal")
    elif p_number_in_list == 11:
        p_graph_settings.bandwidth = p_value
    elif p_number_in_list == 12:
        p_graph_settings.binwidth = p_value
    elif p_number_in_list == 13:
        if p_value.get() == 1:
            args[0].deselect()
    sys.stdout.flush()

def TableOnDoubleClick(event):
    #item = p_table.selection()[0]
    print("you clicked on")#, p_table.item(item, "text"))


#sets form variables
def set_Tk_var():
    global varUConc, varAgebased, varUncorrOrPbc, varErrFilter, varDiscType, varConcType, varEllipseSigma
    global varShowMultiple, varDrawKde, varPosDiscFilter, varNegDiscFilter, varFitDiscordia, varDrawPDP
    global varDrawKDE, varDrawCPDP, varDrawCKDE, varDrawHist, var_pdp_kde_hist, varAnchored, varDiscLinked2Age
    global varKeepPrev, varTypePbc, varShowCalc, varInclude207235Err, varLimitAgeSpectrum, varUncType
    global varCommPb, varMinAgeCrop, varMaxAgeCrop, varAgeCutoff, varDiscCutoff, varKDEBandwidth, varHistBinwidth
    global varAgeAndersen, varDiscPerc
    varUConc = IntVar()
    varUConc.set(1000)
    varDiscType = IntVar()
    varConcType = IntVar()
    varConcType.set(0)
    varEllipseSigma = IntVar()
    varAgebased = IntVar()
    varUncorrOrPbc = IntVar()
    varErrFilter = IntVar()
    varErrFilter.set(5)
    varShowMultiple = IntVar()
    varDrawKde = IntVar()
    varPosDiscFilter = IntVar()
    varPosDiscFilter.set(20)
    varNegDiscFilter = IntVar()
    varNegDiscFilter.set(10)
    varFitDiscordia = IntVar()
    varDrawPDP = IntVar()
    varDrawKDE = IntVar()
    varDrawCPDP = IntVar()
    varDrawCKDE = IntVar()
    varDrawHist = IntVar()
    var_pdp_kde_hist = IntVar()
    varAnchored = IntVar()
    varDiscLinked2Age = IntVar()
    varDiscLinked2Age.set(1)
    varKeepPrev = IntVar()
    varKeepPrev.set(0)
    varTypePbc = IntVar()
    varShowCalc = IntVar()
    varInclude207235Err = IntVar()
    varLimitAgeSpectrum = IntVar()
    varUncType = IntVar()
    varCommPb = IntVar()
    varMinAgeCrop = IntVar()
    varMaxAgeCrop = IntVar()
    varAgeCutoff = IntVar()
    varAgeCutoff.set(1000)
    varDiscCutoff = IntVar()
    varDiscCutoff.set(1000)
    varKDEBandwidth = IntVar()
    varKDEBandwidth.set(50)
    varHistBinwidth = IntVar()
    varHistBinwidth.set(50)
    varAgeAndersen = IntVar()
    varAgeAndersen.set(0)
    varDiscPerc = IntVar()
    varDiscPerc.set(0)


def init(pTop, pGui, *args, **kwargs):
    global w, top_level, root
    w = pGui
    top_level = pTop
    root = pTop


def destroy_window():
    # Function which closes the window.
    global top_level
    top_level.destroy()
    top_level = None


def export_table(p_grainset, p_filters, p_colnames, p_graph_settings, p_filename, p_probcum_filename):
    try:
        file = p_filename
        i = 0
        j = 0
        unc_type = int(p_filters.unc_type)
        an_list = p_grainset.analyses_list
        l_type_pbc = p_filters.use_pbc[0]
        p_grainset.good_bad_sets(p_filters)
        file.write("Analysis name"+',')
        while i < len(p_colnames):
            file.write(p_colnames[i]+',')
            i += 1
        file.write("\n")
        while j < len(p_grainset):
            pbc_corr_buff = pbc_corr(an_list[j], l_type_pbc)
            l_str = ('\n' +
                       str(an_list[j]) + ',' +

                       str(an_list[j].th232_u238[0]) + ',' +
                       str(an_list[j].th232_u238[1]) + ',' +
                       str(an_list[j].th232_u238[2]) + ',' +

                       str(an_list[j].pb208_th232[0]) + ',' +
                       str(an_list[j].pb208_th232[1]) + ',' +
                       str(an_list[j].pb208_th232[2]) + ',' +

                       str(an_list[j].pb207_pb206[0]) + ',' +
                       str(an_list[j].pb207_pb206[1]) + ',' +
                       str(an_list[j].pb207_pb206[2]) + ',' +

                       str(an_list[j].pb207_u235[0]) + ',' +
                       str(an_list[j].pb207_u235[1]) + ',' +
                       str(an_list[j].pb207_u235[2]) + ',' +

                       str(an_list[j].pb206_u238[0]) + ',' +
                       str(an_list[j].pb206_u238[1]) + ',' +
                       str(an_list[j].pb206_u238[2]) + ',' +

                       str(an_list[j].corr_coef_75_68) + ',' +
                       str(an_list[j].corr_coef_86_76) + ',' +

                       str(an_list[j].u_conc[0]) + ',' +
                       str(an_list[j].u_conc[1]) + ',' +
                       str(an_list[j].u_conc[2]) + ',' +

                       str(an_list[j].pbc[0]) + ',' +
                       str(an_list[j].pbc[1]) + ',' +
                       str(an_list[j].pbc[2]) + ',' +

                       str(an_list[j].pb206_pb204[0]) + ',' +
                       str(an_list[j].pb206_pb204[1]) + ',' +
                       str(an_list[j].pb206_pb204[2]) + ',' +

                       str(an_list[j].pb207_pb204[0]) + ',' +
                       str(an_list[j].pb207_pb204[1]) + ',' +
                       str(an_list[j].pb207_pb204[2]) + ',' +

                       str(an_list[j].pb208_pb204[0]) + ',' +
                       str(an_list[j].pb208_pb204[1]) + ',' +
                       str(an_list[j].pb208_pb204[2]) + ',' +

                       str(an_list[j].th232_pb204[0]) + ',' +
                       str(an_list[j].th232_pb204[1]) + ',' +
                       str(an_list[j].th232_pb204[2]) + ',' +

                       str(an_list[j].u238_pb204[0]) + ',' +
                       str(an_list[j].u238_pb204[1]) + ',' +
                       str(an_list[j].u238_pb204[2]) + ',' +

                       str(an_list[j].calc_age(2, p_filters.use_pbc)[0]) + ',' +
                       str(an_list[j].calc_age(2, p_filters.use_pbc)[1]) + ',' +
                       str(an_list[j].calc_age(2, p_filters.use_pbc)[2]) + ',' +

                       str(an_list[j].calc_age(3, p_filters.use_pbc)[0]) + ',' +
                       str(an_list[j].calc_age(3, p_filters.use_pbc)[1]) + ',' +
                       str(an_list[j].calc_age(3, p_filters.use_pbc)[2]) + ',' +

                       str(an_list[j].calc_age(1, p_filters.use_pbc)[0]) + ',' +
                       str(an_list[j].calc_age(1, p_filters.use_pbc)[1]) + ',' +
                       str(an_list[j].calc_age(1, p_filters.use_pbc)[2]) + ',' +

                       str(an_list[j].calc_age(0, p_filters.use_pbc)[0]) + ',' +
                       str(an_list[j].calc_age(0, p_filters.use_pbc)[1]) + ',' +
                       str(an_list[j].calc_age(0, p_filters.use_pbc)[2]) + ',' +
                       g_corr_type[l_type_pbc] + ',' +
                       str(int(pbc_corr_buff[0])) + ',' +
                       str(int(pbc_corr_buff[1])) + ',' +
                       str(int(pbc_corr_buff[2])) + ',' +
                       str(an_list[j].calc_discordance(2, p_filters.disc_type[1])) + ',' +
                       str(an_list[j].calc_discordance(3, p_filters.disc_type[1])) + ',' +

                       str(an_list[j].is_grain_good(p_filters)[0]) + ',' +
                       str(an_list[j].is_grain_good(p_filters)[1]) + ',' +

                       str(an_list[j].calc_age(an_list[j].is_grain_good(p_filters)[1], p_filters.use_pbc)[0]) + ',' +
                       str(an_list[j].calc_age(an_list[j].is_grain_good(p_filters)[1], p_filters.use_pbc)[unc_type]))
            file.write(l_str)
            j += 1
        file.write("\n" * 2)

        file.write("filter by U conc?, "
                   "Uconc filter value (if used), "
                   "correction by Pbc, "
                   "filter by measurement's error, "
                   "measurement's error value (if used), "
                   "Positive discordance filter, "
                   "Negative discordance filter, "
                   "Which age was used, "
                   "age cutoff (if used), "
                   "How was discordance calculated, "
                   "age cutoff for discordance calculation (if used) \n " +
            str(p_filters.filter_by_uconc[0]) + ',' +
                   str(p_filters.filter_by_uconc[1]) + ',' +
                   str(p_filters.use_pbc) + ',' +
                   str(p_filters.filter_by_err[0]) + ',' +
                   str(p_filters.filter_by_err[1]) + ',' +
                   str(p_filters.pos_disc_filter) + ',' +
                   str(p_filters.neg_disc_filter) + ',' +
                   str(p_filters.which_age[0]) + ',' +
                   str(p_filters.which_age[1]) + ',' +
                   str(p_filters.disc_type[0]) + ',' +
                   str(p_filters.disc_type[1]))

        file.close()

        probcum_file = open(p_probcum_filename, "w")
        age = 0
        bandwidth = p_graph_settings.bandwidth
        pdp_list = p_grainset.pdp(int(p_filters.unc_type))
        kde_list = p_grainset.kde(bandwidth)
        cpdp_list = p_grainset.cpdp(int(p_filters.unc_type))
        ckde_list = p_grainset.ckde(bandwidth)

        probcum_file.write("Age" + ',' + "PDP" + ',' + "CPDP" + ',' + "KDE" + ',' + "CKDE" + ',' +
                           'bandwidth=' + str(bandwidth))
        while age < EarthAge:
            probcum_file.write('\n' + str(age) + ',' + str(pdp_list[0][age]) + ',' + str(cpdp_list[age]) + ',' +
                               str(kde_list[0][age]) + ',' + str(ckde_list[age]))
            age += 1
        probcum_file.close()

    except AttributeError:
        pass


def fill_data_table(p_table, p_grainset, p_filters, p_colnames, *args):
    global varAgeAndersen
    if p_grainset.analyses_list == []:
        pass
    for ch in p_table.get_children():
        p_table.delete(ch)
    i = 0
    j = 0
    unc_type = int(p_filters.unc_type)
    an_list = p_grainset.analyses_list
    good_grains = p_grainset.good_bad_sets(p_filters)
    grainset = p_grainset
    filters = p_filters
    l_type_pbc = p_filters.use_pbc[0]
    p_table.heading("#0", text="Analysis name", anchor='c')
    while i < len(p_colnames):
        p_table.heading(p_colnames[i], text=p_colnames[i], anchor='c')
        p_table.column(i, width="100", anchor="c")
        while j < len(grainset):
            th232_u238 = an_list[j].th232_u238
            pb208_th232 = an_list[j].pb208_th232
            pb207_pb206 = an_list[j].pb207_pb206
            pb207_u235 = an_list[j].pb207_u235
            pb206_u238 = an_list[j].pb206_u238
            corr_coef_75_68 = an_list[j].corr_coef_75_68
            corr_coef_86_76 = an_list[j].corr_coef_86_76
            u_conc = an_list[j].u_conc
            pbc = an_list[j].pbc
            pb206_pb204 = an_list[j].pb206_pb204
            pb207_pb204 = an_list[j].pb207_pb204
            pb208_pb204 = an_list[j].pb208_pb204
            th232_pb204 = an_list[j].th232_pb204
            u238_pb204 = an_list[j].u238_pb204
            age_208_232 = an_list[j].age82
            age_207_206 = an_list[j].age76
            age_207_235 = an_list[j].age75
            age_206_238 = an_list[j].age68
            pbc204_age_206_238 = an_list[j].age68_204corr
            pbc204_age_207_235 = an_list[j].age75_204corr
            pbc204_age_208_232 = an_list[j].age82_204corr
            pbc204_age_207_206 = an_list[j].age76_204corr
            pbc207_age = an_list[j].age_207corr
            pbc208_age = an_list[j].age_208corr
            if p_filters.use_pbc[0] == 4:
                pbcAnd_age = an_list[j].calc_age(0, [4, varAgeAndersen.get()])
            else:
                pbcAnd_age = [-1, -1, -1]
            disc_76_68 = 100*an_list[j].calc_discordance(2, p_filters.disc_type[1],  p_filters.use_pbc)
            disc_75_68 = 100*an_list[j].calc_discordance(3, p_filters.disc_type[1], p_filters.use_pbc)
            is_grain_good = an_list[j].is_grain_good(filters)
            best_age = an_list[j].calc_age(is_grain_good[1], p_filters.use_pbc)
            p_table.insert('', 'end', text=(an_list[j]), values=(
                    round(th232_u238[0], 4), round(th232_u238[1], 4), round(th232_u238[2], 4),
                    round(pb208_th232[0], 4), round(pb208_th232[1], 4), round(pb208_th232[2], 4),
                    round(pb207_pb206[0], 4), round(pb207_pb206[1], 4), round(pb207_pb206[2], 4),
                    round(pb207_u235[0], 4), round(pb207_u235[1], 4), round(pb207_u235[2], 4),
                    round(pb206_u238[0], 4), round(pb206_u238[1], 4), round(pb206_u238[2], 4),
                    round(corr_coef_75_68, 2), round(corr_coef_86_76, 2),
                    round(u_conc[0], 4), round(u_conc[1], 4), round(u_conc[2], 4),
                    round(pbc[0], 4), round(pbc[1], 4), round(pbc[2], 4),
                    round(pb206_pb204[0], 1), round(pb206_pb204[1], 1), round(pb206_pb204[2], 1),
                    round(pb207_pb204[0], 1), round(pb207_pb204[1], 1), round(pb207_pb204[2], 1),
                    round(pb208_pb204[0], 1), round(pb208_pb204[1], 1), round(pb208_pb204[2], 1),
                    round(th232_pb204[0], 1), round(th232_pb204[1], 1), round(th232_pb204[2], 1),
                    round(u238_pb204[0], 1), round(u238_pb204[1], 1), round(u238_pb204[2], 1),
                    int(age_208_232[0]), int(age_208_232[1]), int(age_208_232[2]),
                    int(age_207_206[0]), int(age_207_206[1]), int(age_207_206[2]),
                    int(age_207_235[0]), int(age_207_235[1]), int(age_207_235[2]),
                    int(age_206_238[0]), int(age_206_238[1]), int(age_206_238[2]),
                    g_corr_type[l_type_pbc],
                    int(pbc204_age_206_238[0]), int(pbc204_age_206_238[1]), int(pbc204_age_206_238[2]),
                    int(pbc204_age_207_235[0]), int(pbc204_age_207_235[1]), int(pbc204_age_207_235[2]),
                    int(pbc204_age_208_232[0]), int(pbc204_age_208_232[1]), int(pbc204_age_208_232[2]),
                    int(pbc204_age_207_206[0]), int(pbc204_age_207_206[1]), int(pbc204_age_207_206[2]),
                    int(pbc207_age[0]), int(pbc207_age[1]), int(pbc207_age[2]),
                    int(pbc208_age[0]), int(pbc208_age[1]), int(pbc208_age[2]),

                    int(pbcAnd_age[0]), int(pbcAnd_age[1]), int(pbcAnd_age[2]),
                    str(varAgeAndersen.get()),
                    #int(pbc_corr(an_list[j], 4)[0]),
                    #int(pbc_corr(an_list[j], 4)[1]),
                    #int(pbc_corr(an_list[j], 4)[2]),
                    int(disc_76_68), int(disc_75_68),
                    str(is_grain_good[0]), str(is_grain_good[1]),
                    int(best_age[0]), int(best_age[unc_type])
                    ),
                    tags=str(an_list[j].is_grain_good(filters)[0]))

            j += 1
        i += 1

    p_table.tag_configure("False", background="red")
    return good_grains


def set_all_ui_elements(par):
    features_custom_state = [par.chbInclude207235Err, par.entAgeMinCrop, par.entAgeMaxCrop, par.entErrFilter,
                                par.entUconcCutoff, par.cbUConc, par.cbConcType, par.cbErrFilter,
                                par.cbEllipsesAt, par.cbWhichAge,
                                par.cbWhichConc, par.entDiscAgeFixedLim, par.cbPbc, par.entAgeCutoff,
                                par.entHistBinwidth, par.cbDensityPlotType, par.entAgeAndersen, par.cbDiscIntersect]

    for var_frame in (par.frImport, par.frAgeDisc, par.frFilter, par.frGraphSettings, par.frStatus):
        for child in var_frame.winfo_children():
            if child not in features_custom_state:
                child.configure(state=NORMAL)

    par.cbPbc.configure(state="readonly")
    par.cbUConc.configure(state="readonly")
    par.cbConcType.configure(state="readonly")
    par.cbEllipsesAt.configure(state="readonly")
    par.cbDensityPlotType.configure(state="readonly")
    par.entHistBinwidth.configure(state="disabled")
    par.cbErrFilter.configure(state="readonly")
    par.cbDiscIntersect.configure(state="readonly")
    par.cbPbc.configure(state="readonly")
    if par.cbPbc.current() in (0, 1):
        par.entPosDiscFilt.configure(state=NORMAL)
        par.entNegDiscFilt.configure(state=NORMAL)
        par.cbWhichAge.configure(state="readonly")
        par.entAgeCutoff.configure(state=NORMAL)
        par.entDiscAgeFixedLim.configure(state=NORMAL)
        par.cbWhichConc.configure(state="readonly")
        par.entAgeAndersen.configure(state=DISABLED)
    elif par.cbPbc.current() in (2, 3):
        par.entPosDiscFilt.configure(state=DISABLED)
        par.entNegDiscFilt.configure(state=DISABLED)
        par.cbWhichAge.configure(state=DISABLED)
        par.entAgeCutoff.configure(state=DISABLED)
        par.entDiscAgeFixedLim.configure(state=DISABLED)
        par.cbWhichConc.configure(state=DISABLED)
        par.entAgeAndersen.configure(state=DISABLED)
    elif par.cbPbc.current() == 4:
        par.entPosDiscFilt.configure(state=DISABLED)
        par.entNegDiscFilt.configure(state=DISABLED)
        par.cbWhichAge.configure(state=DISABLED)
        par.entAgeCutoff.configure(state=DISABLED)
        par.entDiscAgeFixedLim.configure(state=DISABLED)
        par.cbWhichConc.configure(state=DISABLED)
        par.entAgeAndersen.configure(state=NORMAL)

    if varDiscPerc.get() == 0:
        par.cbDiscIntersect.configure(state=DISABLED)
    else:
        par.cbDiscIntersect.configure(state="readonly")
    if par.cbWhichAge.current() != 1:
            par.entAgeCutoff.configure(state=DISABLED)

    if par.cbWhichConc.current() != 0:
            par.entDiscAgeFixedLim.configure(state=DISABLED)

    if par.cbErrFilter.current() == 1:
        par.entErrFilter.configure(state=NORMAL)
        par.chbInclude207235Err.configure(state=NORMAL)
    else:
        par.entErrFilter.configure(state=DISABLED)
        par.chbInclude207235Err.configure(state=DISABLED)

    if par.cbDensityPlotType.current() == 0:
        par.entKDEBandwidth.configure(state=NORMAL)
        par.entHistBinwidth.configure(state=DISABLED)
    elif par.cbDensityPlotType.current() == 1:
        par.entKDEBandwidth.configure(state=DISABLED)
        par.entHistBinwidth.configure(state=DISABLED)
    else:
        par.entKDEBandwidth.configure(state=DISABLED)
        par.entHistBinwidth.configure(state=NORMAL)

    if par.cbUConc.current() == 1:
        par.entUconcCutoff.configure(state=NORMAL)
    else:
        par.entUconcCutoff.configure(state=DISABLED)

    if varMinAgeCrop.get() == 1:
        par.entAgeMinCrop.configure(state=NORMAL)
    else:
        par.entAgeMinCrop.configure(state=DISABLED)

    if varMaxAgeCrop.get() == 1:
        par.entAgeMaxCrop.configure(state=NORMAL)
    else:
        par.entAgeMaxCrop.configure(state=DISABLED)